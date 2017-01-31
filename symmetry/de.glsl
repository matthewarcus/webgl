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
uniform vec4 params1;
uniform vec4 params2;
uniform int uType;
uniform int uFlags;
uniform samplerCube uCubeMap;
varying vec2 vTextureCoord;

const float phi = 1.618033;
const float phi2 = phi*phi;

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
   float r2=dot(P.xyz,P.xyz);
   vec4 P2=P*P;
   float r4=dot(P2.xyz,P2.xyz);
   vec4 P4=P2*P2;
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
  float time = params2[3];
  float k = 0.1*time;

  // Can't decide between vectors and coordinates.
  vec4 p = vec4(x,y,z,w);
  float cosk = cos(k);
  float sink = sin(k);
  if (cosk < 0.0) {
    //cosk *= -1.0; sink *= -1.0;
  }
  mat4 m = qmat(cosk,-sink,0.0,0.0);
  p = m*p;
#if 1
  x = p.x;
  y = p.y;
  z = p.z;
  w = p.w;
#else  
  x = p.x;
  y = cos(k/10.0) * p.y - sin(k/10.0)*p.z;
  z = sin(k/10.0) * p.y + cos(k/10.0)*p.z;
  w = p.w;
#endif
#if 0
  if (uType == 0) return Endrass8(x,y,z,w);
  else if (uType == 1) return Endrass_8(x,y,z,w);
  else if (uType == 2) return Barth10(x,y,z,w);
  else if (uType == 3) return Sarti12(x,y,z,w);
#else
  if (uType == 0) return Labs(x,y,z,w);
  else if (uType == 1) return Barth(x,y,z,w);
  else if (uType == 2) return Endrass8(x,y,z,w);
  else if (uType == 3) return Sarti12(x,y,z,w);
  else return Sphere(x,y,z,w);
#endif  
}

bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
}

void solve(float x0, float y0, float z0,
           float a, float b, float c) {
  float k0 = 0.0, k1;
  float a0 = Fun(x0,y0,z0,1.0), a1;
  bool bracketed = false;
  float step = 1.0;
  for (int i = 0; i < 100; i++) {
    if (bracketed) {
      float bfact = 0.5;
      // Once we are bracketed, just use bisection
      // Reduce bfact to bias towards start of interval
      if (abs(k1-k0) < 1e-3) break;
      float k2 = (1.0-bfact)*k0 + bfact*k1;
      float a2 = Fun(x0+k2*a,y0+k2*b,z0+k2*c,1.0);
      if (a0*a2 < 0.0) {
        k1 = k2; a1 = a2;
      } else {
        k0 = k2; a0 = a2;
      }
    } else {
      k1 = k0 + step;
      float x = x0+k1*a, y = y0+k1*b, z = z0+k1*c, w = 1.0;
      //float radius = 10.0;
      //if (x*x+y*y+z*z+w*w > radius*radius) break;
      a1 = Fun(x,y,z,w);
      // We can hit exactly 0 - this counts as bracketed.
      if (a0*a1 <= 0.0) {
        bracketed = true;
      } else {
        float step0 = step;
        float maxstep = 1.0;
        // minstep isn't time critical. Small makes
        // for nice cones.
        float minstep = 0.0001;
        step = a1*step/(a0-a1);
        if (step < 0.0) step = -step;
        if (step > maxstep) step = maxstep;
        // Don't grow step by more than 10%
        // A better strategy should be possible
        // Detect overstepping & retreat.
        if (step < minstep) step = minstep;
        float sfact = 1.1;
        if (step > sfact*step0) step = sfact*step0;
        k0 = k1; a0 = a1;
      }
    }
  }
  vec3 eye = vec3(a,b,c);
  vec3 light = normalize(vec3(0.0,1.0,0.5));
  if (!bracketed) {
    //gl_FragColor = vec4(textureCube(uCubeMap,eye).rgb,1.0);
    discard;
    return;
  }
  // Compute gradient & normal
  float eps = 1e-3;
#if 1
  vec3 n = normalize(vec3(Fun(x0+k0*a+eps,y0+k0*b,z0+k0*c,1.0),
                          Fun(x0+k0*a,y0+k0*b+eps,z0+k0*c,1.0),
                          Fun(x0+k0*a,y0+k0*b,z0+k0*c+eps,1.0)) - a0);
#else
  vec3 n = normalize(vec3(Fun(x0+k0*a+eps,y0+k0*b,z0+k0*c,1.0) - Fun(x0+k0*a-eps,y0+k0*b,z0+k0*c,1.0),
                          Fun(x0+k0*a,y0+k0*b+eps,z0+k0*c,1.0) - Fun(x0+k0*a,y0+k0*b-eps,z0+k0*c,1.0),
                          Fun(x0+k0*a,y0+k0*b,z0+k0*c+eps,1.0) - Fun(x0+k0*a,y0+k0*b,z0+k0*c-eps,1.0)));
#endif
  float ambient = 0.6;
  float diffuse = 1.0-ambient;

  int flags = uFlags;
  bool skymap = nextbit(flags);

  vec3 baseColor;
  if (skymap) {
    baseColor = textureCube(uCubeMap,reflect(eye,n)).rgb;
    //baseColor = textureCube(uCubeMap,n).rgb;
  } else {
    baseColor = vec3(1.0,0.5,0.0);
  }
  if (dot(eye,n) > 0.0) n *= -1.0;
  vec3 color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,n));
  float specular = pow(max(0.0,dot(reflect(light,n),eye)),2.0);
  color += 0.5*specular*vec3(1.0,1.0,1.0);
  gl_FragColor = vec4(color,1.0);
}

void main(void) {
  float xscale = params1[0];       // Width multiplier
  float yscale = params1[1];       // Height multiplier
  //float x0 = cos(time), y0 = sin(time), z0 = 6.0;
  float x0 = 0.0, y0 = 0.0, z0 = 10.0;
  float x = 2.0*(vTextureCoord[0]-0.5)*xscale; // + 2.0*xoffset;
  float y = 2.0*(vTextureCoord[1]-0.5)*yscale; // + 2.0*yoffset;
  float z = 0.0;
  float a = x-x0;
  float b = y-y0;
  float c = z-z0;
  float r = sqrt(a*a + b*b + c*c);
  a /= r; b /= r; c /= r;
  solve(x0,y0,z0,a,b,c);
}
