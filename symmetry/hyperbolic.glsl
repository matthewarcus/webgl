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

precision highp float; // Probably want maximum precision here.
precision highp int;

uniform sampler2D uSampler; // Texture sampler
uniform vec4 ufact, vfact;
uniform vec4 params1, params2;
uniform ivec4 iParams; // Misc flags & settings

varying vec2 vTextureCoord; // Could use gl_FragCoord maybe?

const bool doinvert = true; // Do I need this?

const float PI = 3.141592654;
const float TWOPI = 2.0*PI;
const float SQRT3 = 1.732050808; //sqrt(3.0);

/// Complex arithmetic ///
// File inclusion would be good (PHP?)

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

vec2 cexp(vec2 z) {
  return exp(z.x) * expi(z.y);
}

vec2 csin(vec2 z) {
  float x = z.x, y = z.y;
  return cdiv(cexp(vec2(-y,x))-cexp(vec2(y,-x)), vec2(0,2.0));
}

float atanh(float r) {
  return 0.5*log((1.0+r)/(1.0-r));
}

float tanh(float x) {
  return (exp(2.0*x)-1.0)/(exp(2.0*x)+1.0);
}

// Invert z in circle radius r, centre w = (p,0)
// z -> z - w
// z -> z*r^2/|z|^2
// z -> z + w
vec2 invert(vec2 z, vec2 w, float r2) {
  vec2 z1 = z - w;
  float k = r2/dot(z1,z1);
  return z1*k+w;
}

// Overloading for p on x-axis
vec2 invert(vec2 z, float x, float r2) {
  return invert(z,vec2(x,0),r2);
}

bool tryinvert(inout vec2 z, float p, float r2) {
  vec2 z1 = z - vec2(p,0);
  float d2 = dot(z1,z1);
  if (d2 >= r2) return false;
  z = z1*r2/d2 + vec2(p,0);
  return true;
}

// Compute the radius of the disk.
// p is the distance to the centre of the inversion
// circle for the hyperbolic triangle, r is its radius,
// so use Pythagoras to find the right angle for a tangent
// with the disk.
float diskradius(float p,float r) {
  return sqrt(p*p-r*r);
}

bool tryreflect(inout vec2 z, vec2 norm) {
  float k = dot(z,norm);
  if (k <= 0.0) {
    return false;
  } else {
    z -= 2.0*k*norm;
    return true;
  }
}

bool tryreflect(inout vec2 z, vec2 norm, vec2 w) {
  vec2 z1 = z-w;
  if (!tryreflect(z1,norm)) {
    return false;
  } else {
    z = z1 + w;
    return true;
  }
}


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

int next8bits(inout int n) {
  int n0 = n/256;
  int result = n-256*n0;
  n = n0;
  return result;
}

float square(float x) {
  return x*x;
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
  float uxfact = ufact[1];
  float uyfact = ufact[2];
  float uoffset = ufact[3];

  float vscale = vfact[0];
  float vxfact = vfact[1];
  float vyfact = vfact[2];
  float voffset = vfact[3];

  // We are using 5+8+8+4 = 25 bits for main options
  // Hopefully we will get 32 bit integers
  int flags = iParams[0]; //1 + 2*(0 + 2*(1 + 2*(0 + 2*(0 + 2*(4 + 16*(4 + 16*(1 + 16*(0))))))));

  bool hyperbolic = nextbit(flags);
  bool doubleup = nextbit(flags);
  bool mask = nextbit(flags);
  bool fundamental = nextbit(flags);
  bool chiral = nextbit(flags);
  int P = next8bits(flags) + int(hyperbolic);;
  int Q = next8bits(flags);
  int hplane = next4bits(flags);

  float theta = PI/float(P); // Central angle of triangle
  float phi = PI/float(Q); // Other angle of triangle
  // Need picture of hyperbolic region
  // Third side of hyperbolic triangle is an inversion circle.
  // ODBC are on x-axis, A is height 1 above B, so OBA is a right angle and BA = 1
  // BOA = COA = theta, OAD = phi, CAB = theta+phi
  // Maybe should scale to make radius 1 always.
  float p = cos(theta)/sin(theta) + sin(theta+phi)/cos(theta+phi);
  float r = 1.0/cos(theta+phi);
  float r2 = r*r;

  // norm2 is second radial axis, either x-axis or norm reflected in x-axis
  vec2 norm = vec2(-sin(theta),cos(theta));
  vec2 norm2 = !doubleup ? vec2(0.0,-1.0) : vec2(-sin(theta),-cos(theta));
  vec2 ci = vec2(0.0,1.0); // Complex i

  float xlim = 1.0;//ufact[1];
  float eta = 0.1*vyfact;
  float radius = diskradius(p,r);

  vec2 z = vTextureCoord;
  float x = z.x, y = z.y;

  if (hplane == 0) {
    // No half-plane
    x = -xscale*(x - 0.5);
    y = yscale*(y - 0.5);
  } else if (hplane == 4) {
    // Sides are boundaries
    x = xscale*(1.0-x);
    y = yscale*(y - 0.5);
  } else {
    x = -xscale*(x - 0.5);
    y = yscale*y;
  }

  // Scroll
  //x += xoffset*xscale;
  //y += yoffset*yscale;
  x += xoffset;
  y += yoffset;

  z = vec2(x,y);

  if (hplane == 2) {
    z = csin(PI*vec2(x,y)/xscale); // edges and bottom are boundaries
  } else if (hplane == 3) {
    z = cexp(PI*vec2(x,y)/yscale); // top and bottom are boundaries
  } else if (hplane == 4) {
    z = cexp(PI*vec2(y,x)/xscale); // swap x,y; sides are boundaries
  }

  if (hplane != 0) {
    // Then the half-plane to the disk.
    // Disk isn't the unit disk, so scale appropriately
    z = radius*cdiv(z-ci,z+ci);
  }
  
  if (doinvert) {
    if (abs(eta) < 1e-4) {
      z.x = -z.x;
    } else if (hyperbolic) {
      // Translate by reflecting in a (hyperbolic) line.
      float k = radius/sin(eta);
      z = invert(z,k,k*k-radius*radius);
    } else {
      float k = 1.0/tan(eta);
      z = invert(z,k-sin(eta),k*k);
    }
  }

  // psi is an angle here
  float psi = -uoffset;
  if (!hyperbolic) {
    psi += rrepeat * z.x;
  } else if (dot(z,z) > radius*radius) {
    // Or invert to inside the disk
    if (mask) discard;
  } else {
    // Only apply rrotation inside the disk
    psi += rrepeat*atanh(length(z)/radius);
  }

  bool found = false;
  bool odd = false;
  for (int i = 0; i < 100; i++) {
    // Fundamental region is OAB
    // OA is on x-axis, OB is at angle theta
    // AB is circle for hyperbolic case.
    // norm is normal to OB, norm2 is other radial
    // reflection - either x-axis or reflection of OA.
    // Bizarrely, adding "true &&" here avoids a compilation
    // error in Firefox _and_ Chrome on Windows.
    if (true && !tryreflect(z,norm) && !tryreflect(z,norm2)) {
      if (hyperbolic) {
        if (!tryinvert(z,p,r2)) {
          found = true;
          break;
        }
      } else {
        if (!tryreflect(z,vec2(1.0,0.0),vec2(xlim,0.0))) {
          found = true;
          break;
        }
      }
    }
    if (fundamental) discard;
    odd = !odd;
  }
  if (!found) discard;
  if (chiral && odd) z.y = -z.y;
  
  x = z.x; y = z.y;
  float x0 = cos(psi)*x - sin(psi)*y;
  float y0 = sin(psi)*x + cos(psi)*y;
  x = x0; y = y0;

  //x += uoffset;
  //y = y + 0.5 + voffset;
  // scale texture access through uxfact
  x *= 0.5 * exp(uxfact);
  y *= 0.5 * exp(uxfact);
  // and add a variable offset here?
  x += 0.5;
  y += 0.5;
  //y = y + 0.5 + voffset;
  //if (x < 0.0 || x > 1.0 || y < 0.0 || y > 1.0) discard;
  vec4 texColor = texture2D(uSampler, vec2(x,y));
  gl_FragColor = texColor;
}
