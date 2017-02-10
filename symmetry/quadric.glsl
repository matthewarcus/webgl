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
uniform ivec4 iParams;
uniform vec2 uClock;
uniform sampler2D uSampler;
uniform sampler2D uNoise;
uniform samplerCube uCubeMap;
varying vec2 vTextureCoord;

float clock0;
float clock1;
float coffset;
int ftype;
int stype;
int ctype;

// Things that should be uniforms
//const vec3 defaultColor = vec3(1.0,0.5,0.0);
const vec3 defaultColor = vec3(0.8,1.0,0.8);
//const vec3 defaultColor = vec3(0.2,0.1,0.1);
const vec3 light = normalize(vec3(0.0,1.0,1.0));

bool colorswap = false;
bool applyGamma = false;
bool addNoise = false;

const float two31 = pow(2.0,31.0);
int seed = 0;
float frand() {
  // 32-bit RNG from Numerical Recipes
  seed *= 1664525;
  seed += 1013904223;
  return abs(float(seed))/two31;
  //return mod(abs(float(seed)),65000.0)/65000.0; // Huh?
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

const float phi = 1.618033;
const float phi2 = phi*phi;
const float phi4 = phi2*phi2;

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
mat4 qmat(float p0, float p1, float p2, float p3) {
  // As a matrix mul:
  // p0 p1 p2 p3
  // p1 p0 -p3 p2
  // p2 p3 p0 -p1
  // p3 -p2 p1 p0
  return mat4(p0,-p1,-p2,-p3,
              p1,p0,-p3,p2,
              p2,p3,p0,-p1,
              p3,-p2,p1,p0);
}

float Sphere(vec4 p) {
  // Well, a quadric of some sort, depending
  // on the projection.
  //
  //vec4 a = vec4(1.0,cos(0.3*clock1),cos(0.2*clock1),-1.0);
#if 1
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

float Chmutov10(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return T10(x)+T10(y)+T10(z)+1.0;
}

float Chmutov14(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  return T14(x)+T14(y)+T14(z)+1.0;
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

vec4 ctransform(vec4 p, float offset) {
  float w = cos(offset)*p.w - sin(offset)*p.z;
  float z = sin(offset)*p.w + cos(offset)*p.z;
  p.w = w; p.z = z;
  return p;
}

vec4 rotate(vec4 p, float k) {
  float theta = k*0.2*clock0;
  float x = cos(theta)*p.x - sin(theta)*p.z;
  float z = sin(theta)*p.x + cos(theta)*p.z;
  p.x = x; p.z = z;
  return p;
}

vec4 transform(vec4 p) {
  // This should be done in the driver.
  float theta = 0.2*clock0;
  float x = cos(theta)*p.x - sin(theta)*p.z;
  float z = sin(theta)*p.x + cos(theta)*p.z;
  p.x = x; p.z = z;
  return p;
}

float Fun(float x, float y, float z, float w) {
#if defined BENCHMARK
  //return Barth10(x,y,z,w);
  return Labs(vec4(x,y,z,w));
#else
  vec4 p = vec4(x,y,z,w);
  p = transform(p);
  return Sphere(p);
  if (stype == 0) return Labs(p);
  else if (stype == 1) return Barth(p);
  else if (stype == 2) return Endrass8(p);
  else if (stype == 3) return Barth10(p);
  else if (stype == 4) return Sarti12(p);
  else if (stype == 5) return Chmutov10(p);
  else if (stype == 6) return Endrass_8(p);
  else if (stype == 7) return Chmutov14(p);
  else if (stype == 8) return Clebsch(p);
  else return Sphere(p);
#endif  
}

bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
}

// Not quite sure what's going on here. I think
// we have round-towards-zero. We just want a
// hash value though.
int befuddle(float x) {
  if (x < 0.0) return 2*int(-x)+2;
  else return 2*int(x)+1;
}

vec3 selectColor(vec4 q, vec3 eye, vec3 n) {
  if (ctype == 0) {
    return defaultColor;
  }
  if (ctype == 1) {
    return textureCube(uCubeMap,reflect(eye,n)).rgb;
  }
  // Get projective coordinates
  q = rotate(q,-1.0);
  if (colorswap) q = ctransform(q, 1.570796);
  else if (coffset != 0.0) q = ctransform(q,coffset);

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
    vec3 col = texture2D(uSampler,vec2(dot(p,u),dot(p,v))).rgb;
    if (addNoise) {
      col *= texture2D(uNoise,vec2(dot(p,u),dot(p,v))).rgb;
    }
    return col;
  }
  if (ctype == 3) {
    return textureCube(uCubeMap,normalize(vec3(x/w,y/w,z/w))).rgb;
  }
  // Use a grid
  float R = 10.0;
  float gridx = x/w*R;
  float gridy = y/w*R;
  float gridz = x/w*R;
  if (ctype == 4) {
    return vec3(mod(x/w*R,1.0),
                mod(y/w*R,1.0),
                mod(z/w*R,1.0));
  }
  if (ctype == 5) {
    float K = 1.0;
    return vec3 (abs(0.5*sin(K*x/w)+0.5),
                 abs(0.5*sin(K*y/w)+0.5),
                 abs(0.5*sin(K*z/w)+0.5));
  }
  seed = 12345678;
  seed = mix(seed,befuddle(gridx));
  seed = mix(seed,befuddle(gridy));
  seed = mix(seed,befuddle(gridz));
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
const float initstep = 1.0;
const float camera = 20.0;

vec3 project(vec4 p) {
  return vec3(p.x/p.w,p.y/p.w,p.z/p.w);
}

void solve(vec3 p, vec3 eye) {
  float t;
  float eps = 1e-3;
  bool found = false;
  // az^2 + b^r^2 + d^z = c
  //float a =  1.0, b = cos(clock1), c = 1.0, d = sin(clock1);
  float a =  -1.0, b = 1.0, c =  0.0, d =  0.0; // cone
  //float a =  0.0, b = 1.0, c =  0.0, d = -1.0; // paraboloid
  //float a =  1.0, b = 1.0, c =  1.0, d =  0.0; // sphere at origin
  //float a = -1.0, b = 1.0, c = -1.0, d =  0.0; // hyperboloid two-sheet
  //float a = -1.0, b = 1.0, c =  1.0, d =  0.0; // hyperboloid one-sheet
  // Revolve the surface (x) axis about y-axis
  vec3 n = normalize(rotate(vec4(1,0,0,0),1.0).xyz);
  vec3 r = eye;
  float pp = dot(p,p);
  float pn = dot(p,n);
  float pr = dot(p,r);
  float rn = dot(r,n);
  float A = (a-b)*rn*rn + b;
  float B = (a-b)*pn*rn + 0.5*d*rn + pr;
  float C = (a-b)*pn*pn + b*pp + d*pn - c;
  float disc = B*B-A*C;
  if (disc >= 0.0) {
    float t1,t2;
    disc = sqrt(disc);
    if (B > 0.0) disc = -disc;
    t1 = (-B + disc)/A;
    t2 = C/(A*t1);
    if (t1 > t2) {
      t = t1; t1 = t2; t2 = t;
    }
    float lim = 30000000000.0;
    if (t1 > eps && abs(dot(p+t1*r,n)) < lim) {
      t = t1;
      found = true;
    } else if (t2 > eps && abs(dot(p+t2*r,n)) < lim) {
      t = t2;
      found = true;
    }
  }
  if (!found /*|| x*x+y*y+z*z > radius*radius*/) {
    //gl_FragColor = vec4(textureCube(uCubeMap,eye).rgb,1.0);
    discard;
  }
  p = p+t*r; // The actual point
#if 0
  vec3 n = normalize((M*p).xyz);
#else
  // Now n is the normal to the surface
  n = normalize((2.0*a*dot(p,n)+d)*n + 2.0*b*(p-dot(p,n)*n));
#endif

  float ambient = 0.6;
  float diffuse = 1.0-ambient;

  vec3 baseColor = selectColor(vec4(p,1.0),eye,n);
  //vec3 baseColor = vec3(1.0,0.0,0.0);
  //if (p.z > 0.0) baseColor.g = 1.0;
  //if (!check) baseColor.b = 1.0;
  // Point normal towards eye
  if (dot(eye,n) > 0.0) n *= -1.0;
  vec3 color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,n));
  float specular = pow(max(0.0,dot(reflect(light,n),eye)),4.0);
  if (ctype == 1) specular = 0.0;
  color += 0.7*specular*vec3(1.0,1.0,1.0);
  if (applyGamma) color = sqrt(color);
  //color *= 0.5+0.5*frand();
  gl_FragColor = vec4(color,1.0);
  //gl_FragDepthEXT = 0.5;
}

void main(void) {
  int flags = iParams[0];
  nextbit(flags);
  colorswap = nextbit(flags);
  applyGamma = nextbit(flags);
  addNoise = nextbit(flags);

  ftype = iParams[1];
  stype = iParams[2];
  ctype = iParams[3];

  clock0 = uClock[0];
  clock1 = uClock[1];
  coffset = params2[2]; //rrepeat

  float xscale = params1[0];       // Width multiplier
  float yscale = params1[1];       // Height multiplier
  float x = 2.0*(vTextureCoord[0]-0.5)*xscale; // + 2.0*xoffset;
  float y = 2.0*(vTextureCoord[1]-0.5)*yscale; // + 2.0*yoffset;
  float z = 0.0;

#if 0
  seed = int(gl_FragCoord.x*100.0 + gl_FragCoord.y);
  //seed += int(1234.0*gl_FragCoord.x);
  //seed += int(4567.0*gl_FragCoord.y);
  seed += int(1000000.0*vTextureCoord.x);
  seed += int(1000000.0*vTextureCoord.y);
  vec4 noise = texture2D(uNoise, vTextureCoord);
  seed += int(1000000.0*noise.r);
#endif
  float x0 = 0.0, y0 = 0.0, z0 = camera;
  float a = x-x0;
  float b = y-y0;
  float c = z-z0;
  float r = sqrt(a*a + b*b + c*c);
  a /= r; b /= r; c /= r;
  // Could move ray start to radius limit.
  solve(vec3(x0,y0,z0),vec3(a,b,c));
}
