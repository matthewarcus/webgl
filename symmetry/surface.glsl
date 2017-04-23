// The MIT License (MIT)

// Copyright (c) 2017 Matthew Arcus

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Contains functions computing polynomial for various surfaces
// taken from Abdelaziz Nait Merzouk's Fragmentarium shaders.
// https://plus.google.com/114982179961753756261

//#define BENCHMARK
//#define FAST
//#define QUALITY

//#extension GL_EXT_frag_depth : enable
precision highp float;

uniform vec4 params1;
uniform vec4 params2;
uniform vec4 ufact;
uniform vec4 vfact;
uniform ivec4 iParams;
uniform vec4 uClock;
uniform ivec2 uWindow;
uniform sampler2D uSampler;
uniform sampler2D uNoise;
uniform samplerCube uCubeMap;
uniform mat4 uMatrix;

varying vec2 vTextureCoord;

float clock0;
float clock1;
float clock2;
float clock3;
float coffset;
int ftype;
int stype;
int ctype;

// Things that should be uniforms
const vec3 defaultColor = vec3(0.8,1.0,0.8);
//const vec3 defaultColor = vec3(1.0,0.5,0.0);
//const vec3 defaultColor = vec3(0.2,0.1,0.1);
//const vec3 light = vec3(0.0,0.707,0.707);
vec3 light;

const float ambient = 0.5;
const float diffuse = 1.0-ambient;
bool applyGamma = false;

bool colorswap = false;
bool addNoise = false;
bool doSpecular = true;
bool doDiffuse = true;
bool doHorizon = false;

const float two31 = 2147483648.0;
const float phi = 1.618033;
const float phi2 = phi*phi;
const float phi4 = phi2*phi2;

int seed = 0;

float frand() {
  // 32-bit RNG from Numerical Recipes
  seed *= 1664525;
  seed += 1013904223;
  return abs(float(seed))/two31;
}

int mix(int seed, int n) {
  seed += n;
  seed += seed*(8*1024);
  seed += seed/(128*1024);
  seed += seed*32;
  seed += seed/(256*256);
  seed += seed*2;
  seed += seed/512;
  seed += seed*2;
  return seed;
}

// Two dimensional perpendicular distance
float dist(float x0, float y0, float x1, float y1, float x2, float y2) {
   float a = abs((y2-y1)*x0 - (x2 - x1)*y0 + x2*y1 - y2*x1);
   float b = sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
   return a/b;
}

float perp(float x1, float y1, float x2, float y2) {
   float x0 = x2, y0 = 0.0;
   return dist(x0,y0,x1,y1,x2,y2);
}

// Quaternion multiplication as a matrix.
// w coordinate is real element of quaternion
mat4 qmat(vec4 q) {
  float x = q.x, y = q.y, z = q.z, t = q.w;
  return mat4( t,-z, y, x, 
               z, t,-x, y,
              -y, x, t, z,
              -x,-y,-z, t );
}

float Kummer(vec4 P) {
  float A = sqrt(2.0);
  float mu2 = 0.334 + 3.0*(1.0-cos(0.3*clock2));
  float x = P.z; float y = P.y;
  float z = P.x; float w = P.w;
  float p = w-z-A*x;
  float q = w-z+A*x;
  float r = w+z+A*y;
  float s = w+z-A*y;
  float lambda = (3.0*mu2-1.0)/(3.0-mu2);
  float k = x*x + y*y + z*z - mu2*w*w;
  return k*k-lambda*p*q*r*s;
}

float Borg(vec4 p) {
  return sin(p.x*p.y)+sin(p.y*p.z)+sin(p.x*p.z);
}

float Quadratic(vec4 P) {
  // x = y +- (y^2 - z)^0.5
  float x = P.z, y = P.y, z = P.x;
  return x*x - 2.0*x*y + z;
}

float A2(vec4 P) {
   float x = P.w, y = P.y, z = P.z, w = P.x;
  return w*w*w + x*y*z;
}

float Heart0(vec4 p) {
  float x = p.z; float y = p.y;
  float z = p.x; float w = p.w;
  float k = 2.0*x*x + 2.0*y*y + z*z -1.0;
  return k*k*k - 0.1*x*x*z*z*z - y*y*z*z*z;
}

float Heart(vec4 p) {
  float x = p.x; float y = p.z;
  float z = p.y; float w = p.w;
  float k = x*x + 9.0/4.0*y*y + z*z -1.0;
  return k*k*k - x*x*z*z*z - 9.0/80.0*y*y*z*z*z;
}

float Sphere(vec4 p) {
  // Well, a quadric of some sort, depending
  // on the projection.
  //
#if 1
  //vec4 a = vec4(1.0,cos(0.3*clock1),cos(0.2*clock1),-1.0);
  vec4 a = vec4(1.0,1.0,1.0,-1.0);
#elif 1
  mat4 a = mat4(1.0,0.0,0.0,0.0,
                0.0,cos(0.3*clock1),0.0,0.0,
                0.0,0.0,cos(0.2*clock1),0.0,
                0.0,0.0,0.0,-1.0);
#else
  float t = 0.2*clock1;
  float cost = cos(t), sint = sin(t);
  mat4 a = mat4(1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,cost,-sint,
                0.0,0.0,sint,-cost);
#endif
  return dot(a*p,p);
  //float x = p.x; float y = p.y;
  //float z = p.z; float w = p.w;
  //return a*x*x + b*y*y + c*z*z - d*w*w;
}

float Clebsch(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  float t = x+y+z+w;
  return x*x*x + y*y*y + z*z*z + w*w*w - t*t*t;
}

float Cayley(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return w*x*y + x*y*z + y*z*w + z*w*x;
}

float T9(float x) {
  float x2 = x*x;
  return x*(9.0+x2*(-120.0+x2*(432.0+x2*(-576.0+x2*256.0))));
}
float T10(float x) {
  float x2 = x*x;
  return -1.0 + x2*(50.0 + x2*(-400.0 + x2*(1120.0 + x2*(-1280.0 + x2*512.0))));
}
float Chmutov10(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return T10(x)+T10(y)+T10(z)+1.0;
}

#if 0
// These recurrences take a while to compile.
float T11(float x) {
  return 2.0*x*T10(x) - T9(x);
}
float T12(float x) {
  return 2.0*x*T11(x) - T10(x);
}
float T13(float x) {
  return 2.0*x*T12(x) - T11(x);
}
float T14(float x) {
  return 2.0*x*T13(x) - T12(x);
}
float Chmutov14(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return T14(x)+T14(y)+T14(z)+1.0;
}

float T15(float x) {
  return 2.0*x*T14(x) - T13(x);
}
float T16(float x) {
  return 2.0*x*T15(x) - T14(x);
}
float T17(float x) {
  return 2.0*x*T16(x) - T15(x);
}
float T18(float x) {
  return 2.0*x*T17(x) - T16(x);
}

float Chmutov18(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return T18(x)+T18(y)+T18(z)+1.0;
}

// T9(x) + T9(y) + T9(z) + 1 = 0 where T9(x) =
// 256x 9 − 576x 7 + 432x 5 − 120x 3 + 9x

float Chmutov9(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
   return T9(x)+T9(y)+T9(z)+1.0;
}
// T6(x) + T6(y) + T6(z) = 0 where T6(x) =
// 2x 2 (3−4x 2 ) 2 −1 = 32x^6 −48x^4 +18x^2 −1

float T6(float x) {
   return -1.0+x*x*(18.0+x*x*(-48.0+x*x*32.0));
}

float Chmutov6(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return T6(x)+T6(y)+T6(z)+1.0;
}
#endif

// The following functions for various surfaces are
// taken from Abdelaziz Nait Merzouk's Fragmentarium
// shaders.

float Barth(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  float x2 = x*x; float y2 = y*y;
  float z2 = z*z; float w2 = w*w;
  float A = 4.0*(phi2*x2-y2)*(phi2*y2-z2)*(phi2*z2-x2);
  float B = x2 + y2 + z2 - w2;
  return (1.0+2.0*phi)*w2*B*B - A;
}

float Labs7(vec4 p){
   float a = -0.140106854987125;//the real root of 7*a^3+7*a+1=0
   //Constants
   float a1= -0.0785282014969835;//(-12./7.*a-384./49.)*a-8./7.;
   float a2= -4.1583605922880200;//(-32./7.*a+24./49.)*a-4.; 
   float a3= -4.1471434889655100;//(-4.*a+24./49.)*a-4.;
   float a4= -1.1881659380714800;//(-8./7.*a+8./49.)*a-8./7.; 
   float a5= 51.9426145948147000;//(49.*a-7.)*a+50.;

   float r2= dot(p.xy,p.xy);
   vec4 p2=p*p;
   float U = (p.z+p.w)*r2+(a1*p.z+a2*p.w)*p2.z+(a3*p.z+a4*p.w)*p2.w;
   U = (p.z+a5*p.w)*U*U;
   float P = p.x*((p2.x-3.*7.*p2.y)*p2.x*p2.x+(5.*7.*p2.x-7.*p2.y)*p2.y*p2.y);
   P+= p.z*(7.*(((r2-8.*p2.z)*r2+16.*p2.z*p2.z)*r2)-64.*p2.z*p2.z*p2.z);
   return U-P;
}

float Labs(vec4 p) {
  float x = p.x;
  float y = p.y;
  float z = p.z;
  float w = p.w;
  //float t = z; z = y; y = t;
  //float a = -0.140106854987125;//the real root of 7*a^3+7*a+1=0
  //Constants
  float a1= -0.0785282014969835;//(-12./7.*a-384./49.)*a-8./7.;
  float a2= -4.1583605922880200;//(-32./7.*a+24./49.)*a-4.; 
  float a3= -4.1471434889655100;//(-4.*a+24./49.)*a-4.;
  float a4= -1.1881659380714800;//(-8./7.*a+8./49.)*a-8./7.; 
  float a5= 51.9426145948147000;//(49.*a-7.)*a+50.;
  float x2 = x*x,y2 = y*y,z2 = z*z,w2 = w*w;
  float r2= x2+y2;
  float U = (z+w)*r2+(a1*z+a2*w)*z2+(a3*z+a4*w)*w2;
  U = (z+a5*w)*U*U;
  float P = x*((x2-3.0*7.0*y2)*x2*x2+(5.0*7.0*x2-7.0*y2)*y2*y2);
  P += z*(7.0*(((r2-8.0*z2)*r2+16.0*z2*z2)*r2)-64.0*z2*z2*z2);
  return U-P;
}

float Endrass8(vec4 p){
  float x = p.x;
  float y = p.y;
  float z = p.z;
  float w = p.w;
  float x2 = x*x,y2 = y*y,z2 = z*z,w2 = w*w;
  float r2 = x2+y2;
  float U = 64.0*(x2-w2)*(y2-w2)*((x+y)*(x+y)-2.0*w2)*((x-y)*(x-y)-2.0*w2);
  float V = -4.0*(1.0+sqrt(2.0))*r2*r2+(8.0*(2.0+sqrt(2.0))*z2+2.0*(2.0+7.0*sqrt(2.0))*w2)*r2;
  V = V + z2*(-16.0*z2+8.0*(1.0-2.0*sqrt(2.0))*w2) - (1.0+12.0*sqrt(2.0))*w2*w2;
  return V*V-U;
}

float Endrass_8(vec4 p){
  vec4 p2 = p*p;
  float r2 = dot(p.xy,p.xy);
  float U = 64.*(p2.x-p2.w)*(p2.y-p2.w)*((p.x+p.y)*(p.x+p.y)-2.*p2.w)*((p.x-p.y)*(p.x-p.y)-2.*p2.w);
  float V = -4.*(1.-sqrt(2.))*r2*r2+(8.*(2.-sqrt(2.))*p2.z+2.*(2.-7.*sqrt(2.))*p2.w)*r2;
  V = V + p2.z*(-16.*p2.z+8.*(1.+2.*sqrt(2.))*p2.w) - (1.-12.*sqrt(2.))*p2.w*p2.w;
  return V*V-U;
}

float Barth10(vec4 p){//decic
  float r2 = dot(p.xyz,p.xyz);
  vec4 p2 = p*p;
  float r4 = dot(p2.xyz,p2.xyz);
  vec4 p4 = p2*p2;
  return (8.0*(p2.x-phi4*p2.y)*(p2.y-phi4*p2.z)*(p2.z-phi4*p2.x)*(r4-2.0*((p.x*p.y)*(p.x*p.y)+(p.x*p.z)*(p.x*p.z)+(p.y*p.z)*(p.y*p.z)))+(3.0+5.0*phi)*(r2-p2.w)*(r2-p2.w)*(r2-(2.0-phi)*p2.w)*(r2-(2.0-phi)*p2.w)*p2.w);
}

//   Dodecics
float Sarti12(vec4 p){
  vec4 p2 = p*p;
  vec4 p4 = p2*p2;
  float l1 = dot(p2,p2);
  float l2 = p2.x*p2.y+p2.z*p2.w;
  float l3 = p2.x*p2.z+p2.y*p2.w;
  float l4 = p2.y*p2.z+p2.x*p2.w;
  float l5 = p.x*p.y*p.z*p.w;
  float s10 = l1*(l2*l3+l2*l4+l3*l4), s11 = l1*l1*(l2+l3+l4);
  float s12=l1*(l2*l2+l3*l3+l4*l4),    s51=l5*l5*(l2+l3+l4),  s234=l2*l2*l2+l3*l3*l3+l4*l4*l4;
  float s23p=l2*(l2+l3)*l3,   s23m=l2*(l2-l3)*l3; 
  float s34p=l3*(l3+l4)*l4,       s34m=l3*(l3-l4)*l4; 
  float s42p=l4*(l4+l2)*l2,       s42m=l4*(l4-l2)*l2;
  float Q12=dot(p,p); Q12=Q12*Q12*Q12; Q12=Q12*Q12; 
  float S12=33.*sqrt(5.)*(s23m+s34m+s42m)+19.*(s23p+s34p+s42p)+10.*s234-14.*s10+2.*s11-6.*s12-352.*s51+336.*l5*l5*l1+48.*l2*l3*l4;
  return 22.*Q12-243.*S12;
}

float Roman(vec4 p) {
  float r = 2.0;
  float x = p.x, y = p.y, z = p.z;
  float eps = 0.0; //1e-4; // Fudge factor!
  // Use swizzle?
  return x*x*y*y + y*y*z*z + z*z*x*x + r*r*x*y*z - eps;
}

vec4 ctransform(vec4 p, float offset) {
  float w = cos(offset)*p.w - sin(offset)*p.z;
  float z = sin(offset)*p.w + cos(offset)*p.z;
  p.w = w; p.z = z;
  return p;
}

vec4 transform(vec4 p) {
  float k = 0.1*clock1;
  mat4 m = qmat(vec4(0.0,0.0,sin(k),cos(k)));
  p = m*p;
  return p;
}

float Fun(vec4 p) {
#if defined BENCHMARK
  //return Barth10(p);
  return Labs(p);
#else
  p = transform(p);
  //return Clebsch(p);
  //return Roman(p);
  //return Borg(p);
  //return Sphere(p);
  //return Quadratic(p);
  //return A2(p);
  //return Heart(p);
  if (stype == 0) return Labs(p);
  else if (stype == 1) return Barth(p);
  else if (stype == 2) return Endrass8(p);
  else if (stype == 3) return Barth10(p);
  else if (stype == 4) return Sarti12(p);
  else if (stype == 5) return Chmutov10(p);
  else if (stype == 6) return Endrass8(p);
  else if (stype == 7) return Kummer(p);
  //else if (stype == 8) return Endrass_8(p);
  //else if (stype == 9) return Chmutov14(p);
  //else if (stype == 10) return Roman(p);
  else return Sphere(p);
#endif  
}

bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
}

int gridpoint(float x) {
  return int(floor(x));
}

vec3 selectColor(vec4 q, vec3 eye, vec3 n) {
  float uscale = ufact[0];
  float uxfact = exp(ufact[1]);
  float uyfact = ufact[2];
  float uoffset = ufact[3];

  float vscale = vfact[0];
  float vxfact = vfact[1];
  float vyfact = vfact[2];
  float voffset = vfact[3];
  
  if (ctype == 0) {
    return defaultColor;
  }
  if (ctype == 1) {
    return textureCube(uCubeMap,reflect(eye,n)).rgb;
  }
  // Get projective coordinates
  q = transform(q);
  if (colorswap) q = ctransform(q, 1.570796);
  if (coffset != 0.0) q = ctransform(q,coffset);

  float x = q.x, y = q.y, z = q.z, w = q.w;
  if (ctype == 2) {
    // Octamap
    vec3 c = normalize(vec3(1.0,1.0,1.0));
    vec3 u = normalize(cross(c,vec3(0.0,1.0,0.0)));
    vec3 v = normalize(cross(c,u));
    // Taking abs of all coordinates reflects everything
    // into main octant.
    vec3 p = vec3(abs(x/w),abs(y/w),abs(z/w));
    //vec3 p = normalize(vec3(abs(x),abs(y),abs(z)));
    p /= 1.0/dot(p,c); // Central projection
    vec2 texCoords = uxfact*vec2(dot(p,u)+uoffset,dot(p,v)+voffset);
    vec3 col = texture2D(uSampler,texCoords).rgb;
    if (addNoise) {
      col *= texture2D(uNoise,texCoords).rgb;
    }
    return col;
  }
  if (ctype == 3) {
    return textureCube(uCubeMap,normalize(vec3(x/w,y/w,z/w))).rgb;
  }
  // Use a grid
  float R = 10.0/uxfact;
  if (ctype == 4) {
    return vec3(mod(x/w*R + uoffset,1.0),
                mod(y/w*R + uoffset,1.0),
                mod(z/w*R + uoffset,1.0));
  }
  if (ctype == 5) {
    return vec3 (abs(0.5*sin(uxfact*x/w + uoffset)+0.5),
                 abs(0.5*sin(uxfact*y/w + uoffset)+0.5),
                 abs(0.5*sin(uxfact*z/w + uoffset)+0.5));
  }
  float gridx = x/w*R + uoffset;
  float gridy = y/w*R + uoffset;
  float gridz = x/w*R + uoffset;
  seed = 12345678;
  seed = mix(seed,gridpoint(gridx));
  seed = mix(seed,gridpoint(gridy));
  seed = mix(seed,gridpoint(gridz));
  if (ctype == 6) {
    return vec3(frand(),frand(),frand());
  }
  // Position within grid cube
  float fgridx = 2.0*fract(gridx)-1.0;
  float fgridy = 2.0*fract(gridy)-1.0;
  float fgridz = 2.0*fract(gridz)-1.0;
    
  float gridr = sqrt(fgridx*fgridx +
                     fgridy*fgridy +
                     fgridz*fgridz);
  gridr = clamp(0.8*gridr,0.0,1.0);
  return gridr*vec3(.2,.2,.2) + (1.0 - gridr)*vec3(frand(),frand(),frand());
}

// Solution parameters.
#if defined FAST
const int iterations = 100;    // Maximum number of iterations
const float maxincrease = 1.1; // Largest allowed step increase.
#elif defined QUALITY
const int iterations = 300;    // Maximum number of iterations
const float maxincrease = 1.03; // Largest allowed step increase.
#else
const int iterations = 150;    // Maximum number of iterations
const float maxincrease = 1.06; // Largest allowed step increase.
#endif

const float maxstep = 1.0;     // The largest step that can be taken.
const float minstep = 0.001;  // The smallest step
const float initstep = 0.1;

const float horizon = 20.0; // limit on k

void solve(vec4 p0, vec4 r) {
  float k0 = 0.0, k1;
  float a0 = Fun(p0), a1;
  bool bracketed = false;
  bool found = false;
  float step = initstep;
  vec4 p;
  float expected = 0.0;
  for (int i = 0; i < iterations; i++) {
    if (bracketed) {
      // Once we are bracketed, just use bisection
      if (k1-k0 < minstep) {
        found = true;
        break;
      }
      float k2 = (k0 + k1)/2.0;
      //x = x0+k2*a, y = y0+k2*b, z = z0+k2*c, w = 1.0;
      p = p0+k2*r;
      float a2 = Fun(p);
      if (a0*a2 <= 0.0) {
        k1 = k2; a1 = a2;
      } else {
        k0 = k2; a0 = a2;
      }
    } else {
      k1 = k0 + step;
      if (doHorizon && k1 > horizon) break;
      //x = x0+k1*a, y = y0+k1*b, z = z0+k1*c, w = 1.0;
      p = p0 + k1*r;
      //if (x*x+y*y+z*z > radius*radius) break;
      a1 = Fun(p);
      //The idea here is to try and correct the
      // step size by seeing how close we are to
      // the curve, but it doesn't seem to work
      // very well.
      float q = abs((a1-expected)/(a1+expected));
      if (false && expected != 0.0 && q > 0.25) {
        step *= 0.5;
        expected = a0 + 0.5*(expected - a0);
      } else if (a0*a1 <= 0.0) {
        // We can hit exactly 0 - this counts as bracketed.
        bracketed = true;
      } else {
        float step0 = step;
        step = a1*step/(a0-a1);
        //step = perp(k0,a0,k1,a1); // Another nice idea...
        step = abs(step);
        step = min(step,maxstep);
        // Don't grow step by more than 10%
        // A better strategy should be possible
        // Detect overstepping & retreat.
        step = max(step,minstep);
        step = min(step,maxincrease*step0);
        if (a1 <= a0) expected = 0.0;
        else expected = a1 + step*(a1-a0)/(k1-k0);
        k0 = k1; a0 = a1;
      }
    }
  }
  if (!found) discard;

  // Compute gradient & normal
  // Should probably scale eps here
  float eps = 1e-3;
  vec2 delta = vec2(eps,0.0);
#if 0
  p = p0 + k0*r; // Ensure p corresponds to k0 and a0
  vec3 n = vec3(Fun(p + delta.xyyy), Fun(p + delta.yxyy), Fun(p + delta.yyxy)) - a0;
#else
  // Not sure how much difference this makes
  vec3 n = vec3(Fun(p + delta.xyyy) - Fun(p - delta.xyyy),
                Fun(p + delta.yxyy) - Fun(p - delta.yxyy),
                Fun(p + delta.yyxy) - Fun(p - delta.yyxy));
#endif
  float grad = abs(length(n));
  n = normalize(n);

  vec3 eye = r.xyz;
  vec3 baseColor = selectColor(p,eye,n);
  //if (abs(y) < 0.1) baseColor = vec3(1.0,1.0,0.0);
  // Point normal towards eye
  if (dot(eye,n) > 0.0) n *= -1.0;
  vec3 color = baseColor.xyz*ambient;
  float k = dot(light,n);
  if (doDiffuse && k > 0.0) {
    color += baseColor*diffuse*k;
  }
  if (doSpecular && k > 0.0) {
    float specular = pow(max(0.0,dot(reflect(light,n),eye)),4.0);
    color += 0.5*specular*vec3(1.0,1.0,1.0);
  }
  if (applyGamma) color = sqrt(color);
  //color *= 0.5+0.5*frand();
  gl_FragColor = vec4(color,1.0);
  //gl_FragDepthEXT = 0.5;
}

#if 0
  seed = int(gl_FragCoord.x*1000.0 + gl_FragCoord.y);
  //seed += int(1234.0*gl_FragCoord.x);
  //seed += int(4567.0*gl_FragCoord.y);
  seed += int(1000000.0*vTextureCoord.x);
  seed += int(1000000.0*vTextureCoord.y);
  vec4 noise = texture2D(uNoise, vTextureCoord);
  seed += int(1000000.0*noise.r);
#endif

void main(void) {
  // Set global settings
  int flags = iParams[0];
  doHorizon = nextbit(flags);
  colorswap = nextbit(flags);
  applyGamma = nextbit(flags);
  addNoise = nextbit(flags);
  doSpecular = !nextbit(flags);
  doDiffuse = !nextbit(flags);

  ftype = iParams[1];
  stype = iParams[2];
  ctype = iParams[3];

  clock0 = uClock[0];
  clock1 = uClock[1];
  clock2 = uClock[2];
  clock3 = uClock[3];
  coffset = params2[2]; //rrepeat

  light = normalize(vec3(0.5,1,-1));

  // Projection parameters
  float camera = 10.0;
  float xscale = params1[0];       // Width multiplier
  float yscale = params1[1];       // Height multiplier
  // Make sure width and height are even to keep
  // pixel centre off exact axes.
  float width = float(uWindow.x/2*2);
  float height = float(uWindow.y/2*2);
  float x = xscale*(gl_FragCoord.x - 0.5*width)/width;
  float y = yscale*(gl_FragCoord.y - 0.5*height)/height;

  vec4 p0 = vec4(0.0,0.0,-camera,1.0);
  float z = 2.0; // Fixed distance to projection plane
  vec4 r = normalize(vec4(x,y,z,0.0));
  p0 = uMatrix*p0;
  r = uMatrix*r;
  light = mat3(uMatrix)*light; // Light moves with camera

  // Reverse z direction as most of our surfaces are more
  // interesting from that direction.
  // We could roll this into our matrix.
  p0.z = -p0.z;
  r.z = -r.z;
  light.z = -light.z;
  
  // Could move ray start to radius limit.
  // Move the ray a little from the origin.
  // This avoids trouble if the camera position
  // in on the zero locus (eg. for the Roman surface).
  float eps = 0.1;
  solve(p0+eps*r,r);
}
