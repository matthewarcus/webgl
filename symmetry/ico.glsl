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
uniform vec4 uClock;
uniform sampler2D uSampler;
uniform sampler2D uNoise;
uniform samplerCube uCubeMap;
uniform ivec2 uWindow;
uniform mat4 uMatrix;

varying vec2 vTextureCoord;

float iGlobalTime;
vec2 iResolution;
int stype = 0;
int ctype = 0;
bool oneface = false;

// Utilities
bool nextbit(inout int n) {
  int n0 = n/2;
  bool result = bool(n-2*n0);
  n = n0;
  return result;
}

struct Plane {
  int c;   // Color class
  vec3 n;  // Normal
  float d; // Distance from origin (all the same, in fact)
  vec3 u;  // Main face axis
  vec3 v;  // Secondary face axis
};

Plane planes[20];
vec3 colors[5];

// Set up the planes. plane 0 is "reference" plane, planes 1-10 intersect with
// reference sector of plane 0, planes 11 & 12 also intersect but outside the
// bounded part of the face so don't appear in the classic stellations.
void init() {
  planes[0] = Plane(0,vec3(0.3568, 0, 0.9342), 1.5115, vec3(0, -1, 0), vec3(0.9342, 0, -0.3568));
  planes[1] = Plane(1,vec3(0.5774, 0.5774, 0.5774), 1.5115, vec3(0.809, -0.5, -0.309), vec3(0.1103, 0.6455, -0.7558)); // P
  planes[2] = Plane(2,vec3(0.5774, -0.5774, 0.5774), 1.5115, vec3(0.5, -0.309, -0.809), vec3(0.6455, 0.7558, 0.1103)); // P
  planes[3] = Plane(3,vec3(0.9342, -0.3568, 0), 1.5115, vec3(0, 0, 1), vec3(-0.3568, -0.9342, 0));
  planes[4] = Plane(4,vec3(0.9342, 0.3568, 0), 1.5115, vec3(0, 0, -1), vec3(-0.3568, 0.9342, 0));
  planes[5] = Plane(2,vec3(0.5774, 0.5774, -0.5774), 1.5115, vec3(0.5, 0.309, 0.809), vec3(0.6455, -0.7558, -0.1103));
  planes[6] = Plane(1,vec3(0, -0.9342, 0.3568), 1.5115, vec3(-0.5, 0.309, 0.809), vec3(-0.866, -0.1784, -0.4671));
  planes[7] = Plane(1,vec3(0.3568, 0, -0.9342), 1.5115, vec3(0, 1, 0), vec3(0.9342, 0, 0.3568));
  planes[8] = Plane(0,vec3(0.5774, -0.5774, -0.5774), 1.5115, vec3(0.809, 0.5, 0.309), vec3(0.1103, -0.6455, 0.7558));
  planes[9] = Plane(4,vec3(0, -0.9342, -0.3568), 1.5115, vec3(0.5, -0.309, 0.809), vec3(-0.866, -0.1784, 0.4671));
  planes[10] = Plane(2,vec3(-0.5774, -0.5774, -0.5774), 1.5115, vec3(-0.5, -0.309, 0.809), vec3(-0.6455, 0.7558, -0.1103));
  planes[11] = Plane(0,vec3(0, 0.9342, -0.3568), 1.5115, vec3(-0.5, -0.309, -0.809), vec3(-0.866, 0.1784, 0.4671));
  planes[12] = Plane(3,vec3(-0.5774, -0.5774, 0.5774), 1.5115, vec3(0.309, -0.809, -0.5), vec3(0.7558, -0.1103, 0.6455));
  planes[13] = Plane(3,vec3(0, 0.9342, 0.3568), 1.5115, vec3(0.5, 0.309, -0.809), vec3(-0.866, 0.1784, -0.4671));
  planes[14] = Plane(2,vec3(-0.5774, 0.5774, 0.5774), 1.5115, vec3(-0.5, 0.309, -0.809), vec3(-0.6455, -0.7558, 0.1103));
  planes[15] = Plane(4,vec3(-0.3568, 0, 0.9342), 1.5115, vec3(-0.809, -0.5, -0.309), vec3(0.4671, -0.866, 0.1784));
  planes[16] = Plane(4,vec3(-0.5774, 0.5774, -0.5774), 1.5115, vec3(0.309, 0.809, 0.5), vec3(0.7558, 0.1103, -0.6455));
  planes[17] = Plane(1,vec3(-0.9342, 0.3568, 0), 1.5115, vec3(-0.309, -0.809, -0.5), vec3(-0.1784, -0.4671, 0.866));
  planes[18] = Plane(0,vec3(-0.9342, -0.3568, 0), 1.5115, vec3(-0.309, 0.809, 0.5), vec3(-0.1784, 0.4671, -0.866));
  planes[19] = Plane(3,vec3(-0.3568, 0, -0.9342), 1.5115, vec3(-0.809, 0.5, 0.309), vec3(0.4671, 0.866, -0.1784));

  colors[0] = vec3(1.0,0.0,0.0);
  colors[1] = vec3(1.0,1.0,0.0);
  colors[2] = vec3(0.0,1.0,0.0);
  colors[3] = vec3(0.0,1.0,1.0);
  colors[4] = vec3(0.0,0.0,1.0);
}

const vec3 defaultColor = vec3(0.8,1.0,0.8);
const vec3 light0 = vec3(0.0,0.707,-0.707);
vec3 light; // Rotated light goes here
const float ambient = 0.3;
const float diffuse = 1.0-ambient;
bool applyGamma = false;

vec3 getColor(int face) {
  for (int i = 0; i < 5; i++) {
    if (i == face) return colors[i];
  }
  return defaultColor;
}

// Fold point into the reference sector/fundamental region.
// Return "parity" ie. whether original point is part of
// a mirrored region.
bool fold(inout vec2 p) {
  const vec2 n = vec2(0.5,-0.866);
  for (int i = 0; i < 4; i++) {
    if (dot(p,n) > 0.0) p = p - 2.0*(dot(p,n))*n;
    else if (p.x < 0.0) p.x = -p.x;
    else return (i - i/2*2) == 1; // i%2, WebGL style
  }
  return true;
}

// Boolean array indicating which side of plane i the
// point is. We only use a[1-12] in fact.
bool a[20];

#if 0
// Inlined and unrolled.
void classify(vec3 p) {
  const float D = 1.5115;
  a[0] = dot(p,vec3(0.3568, 0, 0.9342)) <= D;
  a[1] = dot(p,vec3(0.5774, 0.5774, 0.5774)) <= D;
  a[2] = dot(p,vec3(0.5774, -0.5774, 0.5774)) <= D;
  a[3] = dot(p,vec3(0.9342, -0.3568, 0)) <= D;
  a[4] = dot(p,vec3(0.9342, 0.3568, 0)) <= D;
  a[5] = dot(p,vec3(0.5774, 0.5774, -0.5774)) <= D;
  a[6] = dot(p,vec3(0, -0.9342, 0.3568)) <= D;
  a[7] = dot(p,vec3(0.3568, 0, -0.9342)) <= D;
  a[8] = dot(p,vec3(0.5774, -0.5774, -0.5774)) <= D;
  a[9] = dot(p,vec3(0, -0.9342, -0.3568)) <= D;
  a[10] = dot(p,vec3(-0.5774, -0.5774, -0.5774)) <= D;
  a[11] = dot(p,vec3(0, 0.9342, -0.3568)) <= D;
  a[12] = dot(p,vec3(-0.5774, -0.5774, 0.5774)) <= D;
}
#endif

// Convenience functions for Coxeter's faces
bool face0() { return a[2]; }
bool face1() { return a[3]; }
bool face2() { return a[4] && a[6] && !a[3]; }
bool face3() { return a[8] && !a[6]; }
bool face4() { return a[1] && a[8] && !a[4]; }
bool face5() { return !a[8] && a[1] && a[6]; }
bool face6() { return !a[8] && a[4] && a[9]; }
bool face7() { return !a[1] && a[8]; }
bool face8() { return a[4] && !a[9]; }
bool face9() { return !a[8] && !a[1] && a[5]; }
bool face10() { return !a[6] && !a[4] && a[9]; }
bool face11() { return a[7] && !a[4] && !a[9]; }
bool face12() { return a[7] && !a[5]; }
bool face13() { return !a[7] && (a[9] || a[5] && a[10]); }

// Coxeter's f1 region. This is the enantiomorphous part.
bool f1(bool parity) {
  if (parity) {
    return face5() || face6();
  } else {
    return face9() || face10();
  }
}

// Find first intersection for ray p,r
vec3 solve(vec3 p, vec3 r) {
  float t0 = 1e10;
  bool found = false;
  vec3 n;
  vec2 uv0;
  int type;
  for (int i = 0; i < 20; i++) {
    if (oneface && i != 15) continue;
    Plane s = planes[i];
    // (p + tr).n = d = p.n + tr.n
    float t = (s.d - dot(p,s.n))/dot(r,s.n);
    // If no better that already found, continue
    if (t <= 0.0 || t >= t0) continue;
    vec3 q = p + t*r; // Point on plane
    vec3 q0 = q - s.d*s.n; // Move in to plane
    vec2 uv = vec2(dot(q0,s.u),dot(q0,s.v)); // And get plane coordinates
    bool parity = fold(uv); // Fold into fundamental region
    q0 = uv[0]*planes[0].u + uv[1]*planes[0].v; //
    q0 = q0 + planes[0].d*planes[0].n; // And move to plane 0
    // 1: P, continuation to inner base of outer spines
    // 2: P, boundary of original face, shell A
    // 3: Q, low 3-sided pyramids on original face, shell B
    // 4: R, high 5-sided pyramids, part of shell F
    // 5: T, side of single outer spine, shell H
    // 6: Q, side of 2, shell C, 5, shell E
    // 7: line U
    // 8: line T
    // 9: line T'
    // 10:
    // Just need the first 12 planes
    // My Firefox seems to unroll up to 10 iterations,
    // we want this unrolled for speed
#if 0
    // Not clear this is any faster really
    classify(q0);
#else
    const float D = 1.5115;
    for (int j = 1; j <= 10; j++) {
      a[j] = dot(q0,planes[j].n) <= D;
    }
    a[11] = dot(q0,planes[11].n) <= D;
    a[12] = dot(q0,planes[12].n) <= D;
#endif    
    bool accept = false;
    if (stype == 0) {
      accept = accept || (parity? face3() || face6() : face3() || face5() || face9() || face10());
    } else if (stype == 1) {
      accept = face4() || face7() || (parity ? face5() : face9() || face10() || face6());
    } else if (stype == 2) {
      accept = face4() || face8() || face11() || (parity ? face5() || face10() : face6() || face9()); // De1f1g2
    } else if (stype == 3) {
      accept = f1(parity); // f1
    } else if (stype == 4) {
      accept = face8() || face9() || face10();  // Shell F
    } else if (stype == 5) {
      accept = face13(); // Shell H
    } else if (stype == 6) {
      accept = (a[12] && a[5]) || (a[6] && a[11]); // Infinity and beyond!
    } else {
      accept = (!a[5] && !a[10]) || !a[11] || face0();
    }
    //accept = face0(); // A
    //accept = face1(); // B
    //accept = face2(); // C
    //accept = face3() || face4(); // D

    //accept = face5() || face6() || face7(); // E
    //accept = face11() || face12(); // G

    //accept = face3() || face5(); // e1
    //accept = a[5] && a[9] && !a[8]; // f1
    //accept = face10() || face12(); // g1

    //accept = face3() || face6() || face9() || face10(); // e1f1
    //accept = face3() || face6() || face9() || face12(); // e1f1g1
    //accept = face5() || face6() || face9() || face12(); // f1g1

    //accept = face4() || face6() || face7(); // e2
    //accept = face7() || face8(); // f2
    //accept = face8() || face9() || face11(); // g2
    
    //accept = face4() || face6() || face8(); // e2f2
    //accept = face4() || face6() || face9() || face11(); // e2f2g2
    //accept = face7() || face9() || face11(); // f2g2

    //accept = face4() || face5(); // De1
    //accept = face7() || face9() || face10(); // Ef1
    //accept = face8() || face9() || face12(); // Fg1

    //accept = face4() || face6() || face9() || face10(); // De1f1
    //accept = face4() || face6() || face9() || face12(); // De1f1g1
    //accept = face7() || face9() || face12(); // Ef1g1
    // USW
    //accept = a[8] || (parity ? a[9] && a[5] : false);
    //accept = a[4]; // De2f2 (exterior)
    //accept = face8() || face3() || face4() || face2() || face7() || face1() || face0();
    //accept = face8() || face10() || face12() || face0() || face7();
    //accept = face0() || face7() || face9() || face12();
    if (!accept) continue;
    found = true;
    t0 = t; n = -s.n; uv0 = uv; type = s.c;
  }
  if (!found) discard;
  vec3 baseColor;
  if (ctype == 0) {
    baseColor = texture2D(uSampler,uv0).rgb;
  } else if (ctype == 1) {
    baseColor = getColor(type);
  } else {
    baseColor = defaultColor;
  }
  if (dot(r,n) > 0.0) {
    n *= -1.0;
  }
        
  vec3 color = baseColor.xyz*(ambient+diffuse*max(0.0,dot(light,n)));
  float specular = pow(max(0.0,dot(reflect(light,n),r)),4.0);
  color += 0.5*specular*vec3(1.0,1.0,1.0);
  //color += 0.5*specular*baseColor.xyz;
  if (applyGamma) color = sqrt(color);
  return color;
}

void main()
{
  stype = iParams[2];
  ctype = iParams[3];
  int flags = iParams[0];
  nextbit(flags);
  oneface = nextbit(flags);
  //colorswap = nextbit(flags);
  applyGamma = nextbit(flags);
  //addNoise = nextbit(flags);
  bool rotatescene = nextbit(flags);
  float clock0 = uClock[0];
  init();

  float xscale = params1[0];       // Width multiplier
  //float yscale = params1[1];       // Height multiplier
  float t = clock0*0.2;
  float cost = cos(t);
  float sint = sin(t);
  float width = float(uWindow[0]);
  float height = float(uWindow[1]);
  float camera = 10.0+5.0*params1[2];
  float x = xscale*(gl_FragCoord.x - 0.5*width)/width;
  float y = xscale*(gl_FragCoord.y - 0.5*height)/width;
  vec3 p = vec3(0.0, 0.0, -camera); // Eye position
  vec3 r = normalize(vec3(x, y, 1.0)); // Ray

  mat3 rmat = mat3(cost,0,-sint,
                   0,1,0,
                   sint,0,cost);
  light = light0;
  p = mat3(uMatrix) * p;
  r = mat3(uMatrix) * r;
  light = mat3(uMatrix) * light;

  p = rmat * p;
  r = rmat * r;

  if (rotatescene) {
    light *= rmat;
  }
  vec3 color = solve(p,r);
  //vec3 color = vec3(1,1,0);
  gl_FragColor = vec4(color,1.0);
}
