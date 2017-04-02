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
uniform vec4 ufact;
uniform vec4 vfact;
uniform ivec4 iParams;
uniform vec4 uClock;
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
const vec3 defaultColor = vec3(0.8,1.0,0.8);
//const vec3 defaultColor = vec3(1.0,0.5,0.0);
//const vec3 defaultColor = vec3(0.2,0.1,0.1);
const vec3 light = vec3(0.0,0.707,-0.707);

bool colorswap = false;
bool applyGamma = false;
bool addNoise = false;

const float two31 = 2147483648.0;
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

vec4 rotate(vec4 p, float k) {
  float theta = k*0.2*clock0;
  float x = cos(theta)*p.x - sin(theta)*p.z;
  float z = sin(theta)*p.x + cos(theta)*p.z;
  p.x = x; p.z = z;
  return p;
}

bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
}

int gridpoint(float x) {
  if (x < 0.0) return 2*int(x);
  else return 2*int(-x)+1;
}

// m is normal in model coords
// n is normal in world coords
vec3 selectColor(vec4 q, vec3 eye, vec3 mx, vec3 nx) {
  float uscale = ufact[0];
  float uxfact = ufact[1];
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
    return textureCube(uCubeMap,reflect(eye,nx)).rgb;
  }
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
    vec2 texCoords = exp(uxfact)*vec2(dot(p,u)+uoffset,dot(p,v)+voffset);
    vec3 col = texture2D(uSampler,texCoords).rgb;
    if (addNoise) {
      col *= texture2D(uNoise,texCoords).rgb;
    }
    return col;
  }
  if (ctype == 3) {
    //return textureCube(uCubeMap,normalize(vec3(x/w,y/w,z/w))).rgb;
    return textureCube(uCubeMap,mx).rgb;
  }
  // Use a grid
  float R = 10.0;
  float gridx = x/w*R;
  float gridy = y/w*R;
  float gridz = z/w*R;
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
  //seed = mix(seed,gridpoint(gridx));
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

// Solve Ax^2 + Bx + C = 0
bool quadratic(float A, float B, float C,
               out float t0, out float t1) {
  float eps = 1e-4;
  // We could probably come up with a better
  // approximation, but this will do for now.
  // Ax^2 + Bx + C = 0
  // -(b +- sqrt(b^2 - 4ac))/2a
  // Solns: x = b/a, x = c/b
#if 0
  if (abs(A) < eps) {
    if (abs(C/B) < eps) return false;
    t0 = t1 = -C/B;
    return true;
  }
#endif
  float disc = B*B-4.0*A*C;
  //disc = abs(disc);
  //if (disc < 0.0) disc = 0.0;
  if (disc < 0.0) {
    // Can allow for some rounding error around zero.
    // but needs careful calibration.
    //if (disc >= -1e-4) disc = 0.0; else
    return false;
  }
  disc = sqrt(disc);
  if (B > 0.0) disc = -disc;
  t0 = (-B + disc)/(2.0*A);
  t1 = C/(A*t0);
  if (t0 > t1) {
    float t = t0; t0 = t1; t1 = t;
  }
  return true;
}

vec4 qmul(vec4 p, vec4 q) {
  return vec4(cross(p.xyz,q.xyz),0.0) +
    vec4(0.0,0.0,0.0,p.w*q.w-dot(p.xyz,q.xyz)) +
    p.w*q + q.w*p;
}
vec4 conj(vec4 p) {
  return vec4(-p.xyz,p.w);
}
vec3 qrot(vec4 p, vec3 q) {
  return qmul(conj(p), qmul(vec4(q,0.0), p)).xyz;
}
vec3 iqrot(vec4 p, vec3 q) {
  return qmul(p, qmul(vec4(q,0.0), conj(p))).xyz;
}
vec4 exp(vec4 p, float t) {
  return vec4(sin(t)*p.xyz,cos(t));
}

void solve(vec3 p, vec3 r) {
  //vec4 rot = normalize(vec4(1.0,2.0,3.0,4.0));
  // Pretty arbitrary rotation axis
  vec4 rot = normalize(vec4(1.0,1.0,1.0,0.0));
  rot = exp(rot,0.1*clock0);
  // The quadric parameters
  vec3 trans = vec3(0.0,0.0,0.0);
  vec3 P,Q;
  float R;
  if (stype == 0) {
    //vec3 P = vec3(0.0,cos(clock1),sin(clock1));
    // Hyperbolic paraboloid & hyperboloid
    P = vec3(1.0,-1.0,sin(clock1));
    Q = vec3(0.0,0.0,-1.0);
    R = 0.0;
  } else if (stype == 1) {
    P = vec3(sin(clock1),2.0,3.0);
    Q = vec3(1.0,0.0,0.0);
    R = -1.0;
  } else if (stype == 2) {
    //vec3 P = vec3(0.0,cos(clock1),sin(clock1));
    // Transition between hyperbolic cylinder and plane pairs
    // Impressive rounding error for R == 0
    P = vec3(1.0,-1.0,0.0);
    Q = vec3(0.0,0.0,0.0);
    R = sin(clock1);
  } else {
    // Transition between cones & hyperboloids
    P = vec3(1.0,-2.0,3.0);
    Q = vec3(0.0,0.0,0.0);
    R = sin(0.309*clock1);//0.0;
  }
  // Turn into rotation matrix4 for production
  p = qrot(rot,p);
  r = qrot(rot,r);
  p += trans;
  // Gradient is P*p+Q so centre is where P*p = -Q
  // If the surface passes through this point,
  // then we have a singularity.

  // We "rebase" the projection point
  // to the closest point to the origin
  // on the ray. This avoids rounding errors
  // eg. for cones with vertex at origin.
  // We could use the "centre" -Q/P, but
  // have problems with zero division etc.
  // On the other hand, if the intersection point
  // is back at the ray origin, we just move the
  // error back there. Possible solution: solve once,
  // then if the solution is dubious (close to origin
  // or singularity, move the projection point to
  // rough solution & redo).
  // P is typically eg. (1,-1,0) & we evaluate (x^2 - y^2) the
  // hard way. Pr.r = ax^2 + by^2 + cy^2.
  // (1,-1,0)*(x,y).(x,y)
  // (a,b,c)*(x,y,z).(x,y,z) = ax*x+b*y*y+c*z*z

  //P(p+r).(p+r) = Pp.r + 2Pp.r + Pr.r
  //P(p-r).(p-r) = Pp.r - 2Pp.r + Pr.r
  float t = -dot(p,r);
  vec3 p0 = p;
  p += t*r;
  float A = dot(P*r,r);
  float B = dot(2.0*P*p+Q,r);
  float C = dot(P*p + Q,p) + R;

  float eps = 1e-2;
  float t0,t1;
  bool found = false;
  if (quadratic(A,B,C,t0,t1)) {
    if (t+t0 > eps) {
      vec3 p0 = p+t0*r;
      p = p0;
      found = true;
    }
    if (!found && t+t1 > eps) {
      vec3 p1 = p+t1*r;
      p = p1;
      found = true;
    }
  }
  if (!found) discard;

  // Now p is the actual point, set n to the normal to the surface
  vec3 n = normalize(2.0*P*p + Q);
  p -= trans; // Just undo translation, which is in object frame
  vec3 m = n;
  n = iqrot(rot,m);
  r = iqrot(rot,r);
  //p = iqrot(rot,p);
  vec3 baseColor = selectColor(vec4(p,1.0),r,m,n);
  //if (p.z < -10.0) baseColor = vec3(1.0,0.0,0.0);
  float ambient = 0.6;
  float diffuse = 1.0-ambient;

  //vec3 baseColor = vec3(1.0,0.0,0.0);
  //if (p.z > 0.0) baseColor.g = 1.0;
  //if (!check) baseColor.b = 1.0;
  // Point normal towards eye
  if (dot(r,n) > 0.0) n *= -1.0;
  vec3 color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,n));
  float specular = pow(max(0.0,dot(reflect(light,n),r)),4.0);
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

  float camera = -10.0; // Use OpenGL coords - -z is towards viewer

#if 0
  seed = int(gl_FragCoord.x*100.0 + gl_FragCoord.y);
  //seed += int(1234.0*gl_FragCoord.x);
  //seed += int(4567.0*gl_FragCoord.y);
  seed += int(1000000.0*vTextureCoord.x);
  seed += int(1000000.0*vTextureCoord.y);
  vec4 noise = texture2D(uNoise, vTextureCoord);
  seed += int(1000000.0*noise.r);
#endif
  float x0 = 0.0-0.0, y0 = 0.0+0.0, z0 = camera;
  float a = x-x0;
  float b = y-y0;
  float c = z-z0;
  float r = sqrt(a*a + b*b + c*c);
  a /= r; b /= r; c /= r;
  // Could move ray start to radius limit.
  solve(vec3(x0,y0,z0),vec3(a,b,c));
}

// (1,1,1), (0,0,0), -1  // sphere
// (1,1,-1), (0,0,0), -1 // hyperboloid I
// (1,1,-1), (0,0,0), 1  // hyperboloid II (via cone)
// (1,1,0), (0,0,1), 0
// (1,-1,0), (0,0,1), 0
// (1,1,0), (0,0,-1), 0
