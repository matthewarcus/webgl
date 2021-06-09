#version 300 es
precision highp float;
out vec4 outColor;

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

precision highp float;
uniform sampler2D uSampler; // Texture sampler
uniform vec4 ufact;         // Modify u-coord
uniform vec4 vfact;         // Modify v-coord
uniform vec4 params1;
uniform vec4 params2;
uniform vec4 AB;
uniform vec4 CD;

uniform vec4 a;
uniform ivec4 mn0;
uniform ivec4 mn1;
uniform ivec4 iParams;
uniform int uType;
      
in vec2 vTextureCoord; // Could use gl_FragCoord maybe?
//varying lowp vec4 vColor;   // Not used currently

// Useful constants
const float PI = 3.14159265359;
const float TWOPI = 2.0*PI;
const float SQRT3 = 1.73205; //sqrt(3.0);

/// Complex arithmetic ///

// Normal vec2 operations work for
// addition, subtraction and
// multiplication by a scalar.
      
// Multiplication
vec2 cmul(vec2 z0, vec2 z1) {
  float x0 = z0.x; float y0 = z0.y; 
  float x1 = z1.x; float y1 = z1.y;
  return vec2(x0*x1-y0*y1,x0*y1+x1*y0);
}

// Reciprocal
vec2 cinv(vec2 z) {
  float x = z.x; float y = z.y;
  float n = 1.0/(x*x + y*y);
  return vec2(n*x,-n*y);
}

// Division
vec2 cdiv(vec2 z0, vec2 z1) {
  return cmul(z0,cinv(z1));
}

// Exponentiation
vec2 expi(float x) {
  return vec2(cos(x),sin(x));
}

vec2 E(int n, int m, float X, float Y) {
  return expi(TWOPI*(float(n)*X + float(m)*Y));
}

vec2 W(int n, int m, float X, float Y) {
  return 0.333*(E(n,m,X,Y) + E(m,-(n+m),X,Y) + E(-(n+m),n,X,Y));
}

vec2 S(int n, int m, float X, float Y) {
  return 0.25*(E(n,m,X,Y) + E(m,-n,X,Y) + E(-n,-m,X,Y) + E(-m,n,X,Y));
}
#if 0
float atanh(float r) {
  return 0.5*log((1.0+r)/(1.0-r));
}
float tanh(float x) {
  return (exp(2.0*x)-1.0)/(exp(2.0*x)+1.0);
}
#endif
// Rather tediously, ES 1.0 doesn't support bit operations
// on ints, so use these for decoding packed values.
bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
}

int next4bits(inout int n) {
  int n0 = n/16;
  int result = n-16*n0;
  n = n0;
  return result;
}

void main(void) {
  float xscale = params1[0];       // Width multiplier
  float yscale = params1[1];       // Height multiplier
  float xoffset = params1[2];
  float yoffset  = params1[3];

  float ulimit = params2[0];
  float vlimit = params2[1];
  float rrepeat = params2[2];

  float uscale = ufact[0];
  float uxfact = ufact[1]; // 2.0/X is whole width - we repeat the texture twice
  float uyfact = ufact[2];
  float uoffset = ufact[3];

  float vscale = vfact[0];
  float vxfact = vfact[1];
  float vyfact = vfact[2];
  float voffset = vfact[3];

  // x+yi is our complex number
  float x = 2.0*(vTextureCoord[0]-0.5)*xscale; // + 2.0*xoffset;
  float y = 2.0*(vTextureCoord[1]-0.5)*yscale; // + 2.0*yoffset;

  vec2 A = vec2(AB[0],AB[1]);
  vec2 B = vec2(AB[2],AB[3]);
  vec2 C = vec2(CD[0],CD[1]);
  vec2 D = vec2(CD[2],CD[3]);

  int flags = iParams[0];
  bool hexagonal = nextbit(flags);
  int ftype = uType; //next4bits(flags);
  
#if 0
  float r = sqrt(x*x+y*y);
  if (r >= 1.0) {
    x /= r*r;
    y /= r*r;
    r = 1.0/r;
  }
  float h = atanh(r);
  x *= h/r;
  y *= h/r;
#endif  
  vec2 z = vec2(x,y);
  if (ftype == 1) z = cdiv(cmul(A,z)+B,cmul(C,z)+D);
  else if (ftype == 2) z = D + cmul(z,C + cmul(z,B + cmul(z,A)));
  else if (ftype == 3) z = cinv(D + cmul(z,C + cmul(z,B + cmul(z,A))));
  x = z.x;
  y = z.y;
  // Apply chosen periodic functions
  z = vec2(0.0,0.0);
  float X,Y; // Lattice coordinates
  if (hexagonal) {
    Y = 2.0*y/SQRT3;
    X = x + Y/2.0;
    z += a[0] * W(mn0[0],mn0[1],X,Y);
    z += a[1] * W(mn0[2],mn0[3],X,Y);
    z += a[2] * W(mn1[0],mn1[1],X,Y);
    z += a[3] * W(mn1[2],mn1[3],X,Y);
#if 0
    // Mainly for framerate testing
    z += 0.5*a[0] * W(mn0[0],mn0[1],X/2.0,Y/2.0);
    z += 0.5*a[1] * W(mn0[2],mn0[3],X/2.0,Y/2.0);
    z += 0.5*a[2] * W(mn1[0],mn1[1],X/2.0,Y/2.0);
    z += 0.5*a[3] * W(mn1[2],mn1[3],X/2.0,Y/2.0);

    z += 0.3*a[0] * W(mn0[0],mn0[1],X/3.0,Y/3.0);
    z += 0.3*a[1] * W(mn0[2],mn0[3],X/3.0,Y/3.0);
    z += 0.3*a[2] * W(mn1[0],mn1[1],X/3.0,Y/3.0);
    z += 0.3*a[3] * W(mn1[2],mn1[3],X/3.0,Y/3.0);

    z += 0.1*a[0] * W(mn0[0],mn0[1],X/5.0,Y/5.0);
    z += 0.1*a[1] * W(mn0[2],mn0[3],X/5.0,Y/5.0);
    z += 0.1*a[2] * W(mn1[0],mn1[1],X/5.0,Y/5.0);
    z += 0.1*a[3] * W(mn1[2],mn1[3],X/5.0,Y/5.0);
#endif
  } else {
    X = x;
    Y = y;
    z += a[0] * S(mn0[0],mn0[1],X,Y);
    z += a[1] * S(mn0[2],mn0[3],X,Y);
    z += a[2] * S(mn1[0],mn1[1],X,Y);
    z += a[3] * S(mn1[2],mn1[3],X,Y);
  }
      
  // Rotate texture coordinates.
  // Number of complete rotations over the image width.
  float theta = PI*rrepeat*x;
  float u = cos(theta)*z.x - sin(theta)*z.y;
  float v = sin(theta)*z.x + cos(theta)*z.y;
  if ((abs(u) > ulimit || abs(v) > vlimit)) {
    discard;
  } else {
#if 0
    bool swapped = false;
    if (u < 0.0) {
      swapped = true;
      u *= -1.0;
      v *= -1.0;
    }
#endif
    // Map to actual texture coordinates
    vec2 texCoord = vec2(0.5 + 0.5*(u*uscale + uxfact*x + uyfact*y + uoffset),
                         0.5 + 0.5*(v*vscale + vxfact*x + vyfact*y + voffset));
    vec4 texColor = texture(uSampler, texCoord);
#if 0
    if (abs(u) < 0.05) {
      texColor.r = 0.0;
      texColor.g = 0.0;
      texColor.b = 0.0;
    } else if (swapped) {
      texColor.r = 1.0 - texColor.r;
      texColor.g = 1.0 - texColor.g;
      texColor.b = 1.0 - texColor.b;
    }
#endif
    outColor = texColor;
  }
}
