#version 300 es
precision mediump float;
uniform float iTime;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform sampler2D uSampler;
#define iChannel0 uSampler
#define iChannel1 uSampler
#define iChannel2 uSampler
#define iChannel3 uSampler

out vec4 outColor;

void mainImage( out vec4 fragColor, in vec2 fragCoord );

void main() {
  mainImage(outColor, gl_FragCoord.xy);
}

vec4 iChannelResolution[4] = vec4[](vec4(0),vec4(1),vec4(2),vec4(3));

#define texelFetch(a,b,c) (vec4(0))

////////////////////////////////////////////////////////////////////////////////
//
// Coloured Wallpaper. Generate the 12 non-hexagonal wallpaper groups, and
// use symmetries to generate a consistent colouring.
// Matthew Arcus (mla), 2020
//
// <mouse>: change offset into texture
// <up>/<down>: zoom in/out
// <page up>/<page down>: number of colours, 1-8
// <left>/<right>: select groups, default cycles through them all
// c: colouring of regions
// e: enlarge texture
// r: use equilateral rhombus for basic tile
// o: offset into texture
// 1-6: draw various grid lines
// 
////////////////////////////////////////////////////////////////////////////////

bool docolor = true;
bool dolines = true;
float scale = 2.2;
int ncolors = 5;
float lod = 1.0;

float PI = 3.141592654;

// Map to lattice separately as it's useful to have the
// resulting intermediate point.
vec2 lattice(vec2 z, out ivec2 index) {
  vec2 z1 = round(0.5*z);
  index = ivec2(z1);
  z -= 2.0*z1;
  return z;
}

void try(inout vec2 z, inout int index, bool test, vec2 z1) {
  index = (index << 1) | int(test);
  if (test) z = z1;
}

void try(inout float x, inout int index, bool test, float x1) {
  index = (index << 1) | int(test);
  if (test) x = x1;
}

bool alert = false;

// Wallpaper groups - also Conway's orbiform notation.
int p1 = 1; // o
int p2 = 2; // 2222
int pm = 3; // **
int pg = 4; // xx
int cm = 5; // *x
int pmm = 6; // *2222
int pmg = 7; // 22*
int pgg = 8; // 22x
int cmm = 9; // 2*22
int p4 = 10; // 442
int p4m = 11; //*442
int p4g = 12; // 4*2

// How many bits for the index.
int nbits[12] = int[](0,1,1,1,1,2,2,2,2,2,3,3);

vec2 wallpaper(vec2 z, int type, out int index) {
  index = 0;
  if (type == p1) {
    // Nothing to do
  } else if (type == p2 || type == pgg) {
    // 2-fold rotation
    try(z,index,z.x < 0.0,-z); // Rotate 180 about 0
    if (type == pgg) {
      // (viii) reflect up, glide, then reflect back - pgg
      bool p = z.y < 0.0;
      if (p) z.y = -z.y;
      try(z,index,z.x+z.y > 1.0,vec2(1.0-z.x,z.y-1.0));
      if (p) z.y = -z.y;
    }
  } else if (type == pm || type == pmg) {
    try(z.x,index,z.x < 0.0,-z.x); // Reflect in x=0 - pm
    if (type == pmg) {
      // (vii) reflect in x, then rotate about (0.5,0) - pmg
      try(z,index,z.y < 0.0, vec2(1.0-z.x,-z.y));
    }
  } else if (type == pg) {
    try(z,index,z.x < 0.0, vec2(z.x+1.0,-z.y)); // Glide along y=0
  } else if (type == cm || type == cmm) {
    try(z,index,z.x+z.y < 0.0,-z.yx); // reflect in x+y=0
    if (type == cmm) {
      // (ix) also reflect in x = y - cmm
      try(z,index,z.y > z.x,z.yx);
    }
  } else if (type == pmm || type == p4m) {
    index |= 2*int(z.x < 0.0) + int(z.y < 0.0); // 2 bits
    z = abs(z); // (vi) reflect quadrants - pmm
    if (type == p4m) {
      // (xi) 4 reflection axes - p4m
      try(z,index,z.x + z.y > 1.0, 1.0-z.yx); // Reflect in x+y=1
    }
  } else if (type == p4 || type == p4g) {
    // 4-fold rotation
    // (x) rotate quadrants - p4
    index |= 2*int(z.x < 0.0) + int(z.y < 0.0);
    if (z.x < 0.0) { z = -z; } // Rotate 180 about 0
    if (z.y < 0.0) z = vec2(-z.y,z.x); // Rotate 90 about 0
    if (type == p4g) {
      // (xii) additionally: reflect in x+y=1 - p4g
      try(z,index,z.x+z.y > 1.0,1.0-z.yx); // Reflect in x+y=1
    }
  } else {
    alert = true;
  }
  return z;
}

float linedist(vec2 pos, vec2 a, vec2 b) {
  vec2 pa = pos - a;
  vec2 ba = b - a;
  float h = dot(pa, ba) / dot(ba, ba);
  float d = length(pa - ba * h);
  return d;
}

int mymod(int n, int m) {
  // glsl % undefined for -ve arguments
  if (n < 0) return m-1-(-n-1)%m;
  else return n%m;
}

vec3 indexcolor(int index) {
  index = mymod(index, ncolors);
  vec3 col = vec3(0);
  if (index == 0) col = vec3(1,0,0);
  if (index == 1) col = vec3(1,1,0);
  if (index == 2) col = vec3(0,1,0);
  if (index == 3) col = vec3(0,0,1);
  if (index == 4) col = vec3(0,1,1);
  if (index == 5) col = vec3(1,0,1);
  if (index == 6) col = vec3(1);
  if (index == 7) col = vec3(0.5);
  return 0.5+0.5*col;
}

bool dooffset = true;
bool doenlarge = true;

vec3 getcolor(vec2 z, int type, ivec4 ivector) {
  vec2 offset = vec2(0);
  if (dooffset) offset += vec2(sin(0.1*iTime),sin(0.123*iTime));
  if (iMouse.x > 0.0) {
    offset += 2.0*(2.0*iMouse.xy-iResolution.xy)/iResolution.y;
  }
  vec2 uv = 0.5*z-0.5+offset;
  if (doenlarge) uv *= 0.5;
  vec4 tt = textureLod(iChannel0,uv,lod);
  if (!docolor) return tt.xyz;
  // Generate a colour from the index. Lots of ways this could
  // be done, eg. just use z and w so whole tile gets same colour.
  int index = 0;
  int n = nbits[type];
  if (n == 3) {
    index = ivector.x;
    // xor the top bit (bit 3) with the other bits
    index ^= -((index>>2)&1);
    index &= 3;
    // offset for the tile
    index += ivector.z + ivector.w;
  } else {
    index = ((ivector.z + ivector.w)<<n) - ivector.x;
  }
  vec3 color = indexcolor(index);
  color *= length(tt.xyz);
  return color;
}

const int CHAR_0 = 48;
const int CHAR_C = 67;
const int CHAR_E = 69;
const int CHAR_F = 70;
const int CHAR_H = 72;
const int CHAR_M = 77;
const int CHAR_O = 79;
const int CHAR_R = 82;
const int CHAR_S = 83;
const int CHAR_T = 84;
const int CHAR_U = 85;
const int KEY_PAGE_UP = 33;
const int KEY_PAGE_DOWN = 34;
const int KEY_LEFT = 37;
const int KEY_UP = 38;
const int KEY_RIGHT = 39;
const int KEY_DOWN = 40;

vec4 store(int i,int j) {
  return texelFetch(iChannel3, ivec2(i,j),0);
}

int keycount(int key) {
  return int(store(0,key).x);
}

bool key(int key) {
  return texelFetch(iChannel2,ivec2(key,2),0).x != 0.0;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  int type = mymod(keycount(KEY_RIGHT)-keycount(KEY_LEFT),13);
  if (type == 0) type = 1+int(0.2*iTime)%12;
  ncolors = clamp(5+keycount(KEY_PAGE_DOWN)-keycount(KEY_PAGE_UP),1,8);
  docolor = !key(CHAR_C);
  dooffset = !key(CHAR_O);
  doenlarge = key(CHAR_E);
  scale *= exp(0.1*float(keycount(KEY_DOWN)-keycount(KEY_UP)));
  lod = log2(scale*iChannelResolution[0].y/iResolution.y);
  if (doenlarge) lod -= 1.0;
  vec2 z = (2.0*fragCoord-iResolution.xy)/iResolution.y;
  z *= scale;
  mat2 m = mat2(1);
  //m = mat2(-0.8,0.5,0.8,0.5);
  // Use rhombus for basic tile - not all groups work with this.
  if (key(CHAR_R)) m = mat2(1,0,0.5,0.5*sqrt(3.0));
  z = inverse(m)*z;
  ivec4 index;
  vec2 z0 = z;
  z = lattice(z,index.zw);
  vec2 z1 = m*z;
  z = wallpaper(z,type,index.x);
  z = m*z;
  //if (z != z0) alert = true; // Show region
  vec3 col = getcolor(z,type,index);
  if (dolines) {
    // Draw various grid lines.
    float d = 1e8;
    vec2 p = m*vec2(1,1);
    vec2 q = m*vec2(1,-1);
    vec2 r = m*vec2(-1,-1);
    vec2 s = m*vec2(-1,1);
    vec2 pq = 0.5*(p+q);
    vec2 rs = 0.5*(r+s);
    vec2 qr = 0.5*(q+r);
    vec2 sp = 0.5*(s+p);
    if (!key(CHAR_0+1)) {
      d = min(d,linedist(z1,p,q));
      d = min(d,linedist(z1,q,r));
      d = min(d,linedist(z1,r,s));
      d = min(d,linedist(z1,s,p));
    }
    if (key(CHAR_0+2)) d = min(d,linedist(z1,pq,rs));
    if (key(CHAR_0+3)) d = min(d,linedist(z1,qr,sp));
    if (key(CHAR_0+4)) d = min(d,linedist(z1,p,r));
    if (key(CHAR_0+5)) d = min(d,linedist(z1,q,s));

    if (key(CHAR_0+6)) {
      d = min(d,linedist(z1,pq,qr));
      d = min(d,linedist(z1,qr,rs));
      d = min(d,linedist(z1,rs,sp));
      d = min(d,linedist(z1,sp,pq));
    }

    col *= smoothstep(0.01,0.02,d);
  }
  if (alert) col *= 0.6;
  fragColor = vec4(col,1);
}
