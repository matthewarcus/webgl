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

precision highp float;

//#define BENCHMARK
//#define FAST
#define QUALITY

uniform vec4 params1;
uniform vec4 params2;
uniform int uType;
uniform int uFlags;
uniform samplerCube uCubeMap;
varying vec2 vTextureCoord;

// Things that should be uniforms
//const vec3 defaultColor = vec3(1.0,0.5,0.0);
//const vec3 defaultColor = vec3(0.8,1.0,0.8);
const vec3 defaultColor = vec3(0.2,0.1,0.1);
const vec3 light = normalize(vec3(0.0,1.0,1.0));

const float phi = 1.618033;
const float phi2 = phi*phi;

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

float Sphere(float x, float y, float z, float w) {
  return x*x + y*y + z*z + w*w - 4.0;
}

float Clebsch(float x, float y, float z, float w) {
  float t = x+y+z+w;
  return x*x*x + y*y*y + z*z*z + w*w*w - t*t*t;
}

float Cayley(float x, float y, float z, float w) {
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

float Chmutov10(float x, float y, float z, float w) {
   return T10(x)+T10(y)+T10(z)+1.0;
}

float Chmutov14(float x, float y, float z, float w) {
  return T14(x)+T14(y)+T14(z)+1.0;
}

float Chmutov18(float x, float y, float z, float w) {
   return T18(x)+T18(y)+T18(z)+1.0;
}

// T9(x) + T9(y) + T9(z) + 1 = 0 where T9(x) =
// 256x 9 − 576x 7 + 432x 5 − 120x 3 + 9x

float Chmutov9(float x, float y, float z, float w) {
   return T9(x)+T9(y)+T9(z)+1.0;
}
// T6(x) + T6(y) + T6(z) = 0 where T6(x) =
// 2x 2 (3−4x 2 ) 2 −1 = 32x^6 −48x^4 +18x^2 −1

float T6(float x) {
   return -1.0+x*x*(18.0+x*x*(-48.0+x*x*32.0));
}

float Chmutov6(float x, float y, float z, float w) {
   return T6(x)+T6(y)+T6(z)+1.0;
}

// The following functions for various surfaces are
// taken from Abdelaziz Nait Merzouk's Fragmentarium
// shaders.

float Barth(float x, float y, float z, float w) {
  float x2 = x*x;
  float y2 = y*y;
  float z2 = z*z;
  float w2 = w*w;
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

   float	r2= dot(p.xy,p.xy);
   vec4 p2=p*p;
   float U = (p.z+p.w)*r2+(a1*p.z+a2*p.w)*p2.z+(a3*p.z+a4*p.w)*p2.w;
   U = (p.z+a5*p.w)*U*U;
   float P = p.x*((p2.x-3.*7.*p2.y)*p2.x*p2.x+(5.*7.*p2.x-7.*p2.y)*p2.y*p2.y);
   P+= p.z*(7.*(((r2-8.*p2.z)*r2+16.*p2.z*p2.z)*r2)-64.*p2.z*p2.z*p2.z);
   return U-P;
}

float Labs(float x, float y, float z, float w) {
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

float Endrass8(float x, float y, float z, float w){
  float x2 = x*x,y2 = y*y,z2 = z*z,w2 = w*w;
  float r2 = x2+y2;
  float U = 64.0*(x2-w2)*(y2-w2)*((x+y)*(x+y)-2.0*w2)*((x-y)*(x-y)-2.0*w2);
  float V = -4.0*(1.0+sqrt(2.0))*r2*r2+(8.0*(2.0+sqrt(2.0))*z2+2.0*(2.0+7.0*sqrt(2.0))*w2)*r2;
  V = V + z2*(-16.0*z2+8.0*(1.0-2.0*sqrt(2.0))*w2) - (1.0+12.0*sqrt(2.0))*w2*w2;
  return V*V-U;
}

float Endrass_8(float x, float y, float z, float w){
  vec4 p = vec4(x,y,z,w);
  vec4 p2 = p*p;
	float r2 = dot(p.xy,p.xy);
	float U = 64.*(p2.x-p2.w)*(p2.y-p2.w)*((p.x+p.y)*(p.x+p.y)-2.*p2.w)*((p.x-p.y)*(p.x-p.y)-2.*p2.w);
	float V = -4.*(1.-sqrt(2.))*r2*r2+(8.*(2.-sqrt(2.))*p2.z+2.*(2.-7.*sqrt(2.))*p2.w)*r2;
	V = V + p2.z*(-16.*p2.z+8.*(1.+2.*sqrt(2.))*p2.w) - (1.-12.*sqrt(2.))*p2.w*p2.w;
	return V*V-U;
}

#define PHI  1.618034
#define PHI2 2.618034
#define PHI4 6.854102

float Barth10(float x, float y, float z, float w){//decic
  vec4 P = vec4(x,y,z,w);
  float r2 = dot(P.xyz,P.xyz);
  vec4 P2 = P*P;
  float r4 = dot(P2.xyz,P2.xyz);
  vec4 P4 = P2*P2;
  return (8.0*(P2.x-PHI4*P2.y)*(P2.y-PHI4*P2.z)*(P2.z-PHI4*P2.x)*(r4-2.0*((P.x*P.y)*(P.x*P.y)+(P.x*P.z)*(P.x*P.z)+(P.y*P.z)*(P.y*P.z)))+(3.0+5.0*PHI)*(r2-P2.w)*(r2-P2.w)*(r2-(2.0-PHI)*P2.w)*(r2-(2.0-PHI)*P2.w)*P2.w);
}

//   Dodecics
float Sarti12(float x, float y, float z, float w){
  vec4 p = vec4(x,y,z,w);
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

float Fun(float x, float y, float z, float w) {
#if defined BENCHMARK
  return Barth10(x,y,z,w);
  //return Labs7(vec4(x,y,z,w));
#else
  float time = params2[3];
  float k = 0.1*time;

  // Can't decide between vectors and coordinates.
  vec4 p = vec4(x,y,z,w);
  float cosk = cos(k);
  float sink = sin(k);
  mat4 m = qmat(cosk,-sink,0.0,0.0);
  p = m*p;
  x = p.x;
  y = p.y;
  z = p.z;
  w = p.w;

  if (uType == 0) return Labs(x,y,z,w);
  else if (uType == 1) return Barth(x,y,z,w);
  else if (uType == 2) return Endrass8(x,y,z,w);
  else if (uType == 3) return Barth10(x,y,z,w);
  else if (uType == 4) return Sarti12(x,y,z,w);
  else if (uType == 5) return Chmutov10(x,y,z,w);
  else if (uType == 6) return Endrass_8(x,y,z,w);
  else if (uType == 7) return Chmutov14(x,y,z,w);
  else return Sphere(x,y,z,w);
#endif  
}

bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
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
const float camera = 10.0;
const float radius = 6.0; // Restrict view to this (unused)
const float horizon = 20.0; // limit on k
void solve(float x0, float y0, float z0,
           float a, float b, float c) {
  float k0 = 0.0, k1;
  float a0 = Fun(x0,y0,z0,1.0), a1;
  bool bracketed = false;
  bool found = false;
  float step = initstep;
  float x, y, z, w;
  float expected = 0.0;
  for (int i = 0; i < iterations; i++) {
    if (bracketed) {
      // Once we are bracketed, just use bisection
      if (k1-k0 < minstep) {
        found = true;
        break;
      }
      float k2 = (k0 + k1)/2.0;
      x = x0+k2*a, y = y0+k2*b, z = z0+k2*c, w = 1.0;
      float a2 = Fun(x,y,z,w);
      if (a0*a2 <= 0.0) {
        k1 = k2; a1 = a2;
      } else {
        k0 = k2; a0 = a2;
      }
    } else {
      k1 = k0 + step;
      if (k1 > horizon) break;
      x = x0+k1*a, y = y0+k1*b, z = z0+k1*c, w = 1.0;
      //if (x*x+y*y+z*z > radius*radius) break;
      a1 = Fun(x,y,z,w);
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
  vec3 eye = vec3(a,b,c);
  if (!found /*|| x*x+y*y+z*z > radius*radius*/) {
    //gl_FragColor = vec4(textureCube(uCubeMap,eye).rgb,1.0);
    discard;
  }
  // Compute gradient & normal
  float eps = 1e-3;
#if 1
  vec3 n = vec3(Fun(x0+k0*a+eps,y0+k0*b,z0+k0*c,1.0),
                Fun(x0+k0*a,y0+k0*b+eps,z0+k0*c,1.0),
                Fun(x0+k0*a,y0+k0*b,z0+k0*c+eps,1.0)) - a0;
#else
  vec3 n = vec3(Fun(x0+k0*a+eps,y0+k0*b,z0+k0*c,1.0) - Fun(x0+k0*a-eps,y0+k0*b,z0+k0*c,1.0),
                Fun(x0+k0*a,y0+k0*b+eps,z0+k0*c,1.0) - Fun(x0+k0*a,y0+k0*b-eps,z0+k0*c,1.0),
                Fun(x0+k0*a,y0+k0*b,z0+k0*c+eps,1.0) - Fun(x0+k0*a,y0+k0*b,z0+k0*c-eps,1.0));
#endif
  float grad = abs(length(n));
  n = normalize(n);
  float ambient = 0.6;
  float diffuse = 1.0-ambient;

  int flags = uFlags;
  bool skymap = nextbit(flags);

  vec3 baseColor;
  if (skymap) {
    baseColor = textureCube(uCubeMap,reflect(eye,n)).rgb;
  } else {
    baseColor = defaultColor;
  }
  //baseColor.r = min(5.0*grad,1.0);
  if (dot(eye,n) > 0.0) n *= -1.0;
  vec3 color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,n));
  float specular = pow(max(0.0,dot(reflect(light,n),eye)),4.0);
  color += 0.7*specular*vec3(1.0,1.0,1.0);
  gl_FragColor = vec4(color,1.0);
}

void main(void) {
  float time = params2[3];
  float xscale = params1[0];       // Width multiplier
  float yscale = params1[1];       // Height multiplier
  float x0 = 0.0, y0 = 0.0, z0 = camera;
  float x = 2.0*(vTextureCoord[0]-0.5)*xscale; // + 2.0*xoffset;
  float y = 2.0*(vTextureCoord[1]-0.5)*yscale; // + 2.0*yoffset;
  float z = 0.0;
  float a = x-x0;
  float b = y-y0;
  float c = z-z0;
  float r = sqrt(a*a + b*b + c*c);
  a /= r; b /= r; c /= r;
  // Could move ray start to radius limit.
  solve(x0,y0,z0,a,b,c);
}
