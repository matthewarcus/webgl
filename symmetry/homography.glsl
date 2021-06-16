#version 300 es
precision highp float;
uniform float iTime;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform sampler2D uSampler;
#define iChannel0 uSampler
#define LOCAL

uniform uvec4 uKeysPressed;
uniform uvec4 uKeysToggled;

out vec4 outColor;

void mainImage( out vec4 fragColor, in vec2 fragCoord );

void main() {
  mainImage(outColor, gl_FragCoord.xy);
}

#define texelFetch(a,b,c) (vec4(0))
//#define key(code) (texelFetch(iChannel3, ivec2((code),2),0).x != 0.0)
#define key(code) ((uKeysToggled[code>>5]&(1U<<(code&31))) != 0U)

////////////////////////////////////////////////////////////////////////////////

const float PI = 3.14159265359;

bool alert = false;
void assert(bool b) {
  if (!b) alert = true;
}

const int CHAR_X = 88;
const int CHAR_Y = 89;
const int CHAR_Z = 90;

// Project from (0,1) to y = 0 plane
// Return homogeneous coordinates
vec2 stereographic(vec2 z) {
  if (z.y == 1.0) return vec2(1,0);
  else return vec2(z.x,1.0-z.y);
}

// Inverse projection from y = 0 (as a homogeneous coordinate) to unit sphere
// with (0,1) as projection point.
vec2 istereographic(vec2 z) {
  if (z.y == 0.0) return vec2(0,1);
  // |k(x,0)+(1-k)(0,1)| = |kx,1-k| = 1 => k²x²+1-2k+k² = 1 => k(x²+1) = 2
  float x = z.x/z.y;
  float k = 2.0/(x*x+1.0);
  return vec2(k*x,1.0-k);
}

// Find a projective mapping taking p0,p1,p2,p4 to
// triangle of reference and unit point, ie:
// p0 -> (1,0,0), p1 -> (0,1,0), p2 -> (0,0,1), p3 -> (1,1,1)
// No three points collinear.
mat3 rproject(vec3 p0, vec3 p1, vec3 p2, vec3 p3) {
  // Just an inverse for the first three points
  // (the triangle of reference). No inverse if collinear.
  mat3 m = inverse(mat3(p0,p1,p2)); // column major!
  vec3 p3a = m*p3;
  // Then scale each row so the unit point (1,1,1) is correct
  m = transpose(m);
  // zero components here only if not collinear
  m[0] /= p3a[0];
  m[1] /= p3a[1];
  m[2] /= p3a[2];
  m = transpose(m);
  return m;
}

// Return matrix that maps p0,p1,p2 to (1,0),(0,1) and (1,1)
mat2 rproject(vec2 p0, vec2 p1, vec2 p2) {
  mat2 m = inverse(mat2(p0,p1)); // column major!
  vec2 p2a = m*p2;
  // Then scale each row so the unit point (1,1) is correct
  m = transpose(m);
  // zero components here only if not collinear
  m[0] /= p2a[0];
  m[1] /= p2a[1];
  m = transpose(m);
  return m;
}

vec2 circle(float t) {
  return vec2(sin(t),cos(t));
}

vec2 icircle0(vec2 z) {
  return vec2(atan(z.y,z.x),length(z));
}

vec2 norm(vec3 p) {
  return p.xy/p.z;
}

// Screen coords to P2 coords
vec2 map(vec2 p) {
  return (2.0*p - iResolution.xy) / iResolution.y;
}

vec2 cmul(vec2 z1, vec2 z2) {
  return mat2(z1.x,z1.y,-z1.y,z1.x)*z2;
}

vec2 cinv(vec2 z) {
  return vec2(z.x,-z.y)/dot(z,z);
}

vec2 cdiv(vec2 z1, vec2 z2) {
  return cmul(z1,cinv(z2));
}

// Invert unit circle to upper half plane
// Centre of inversion (0,-1), r² = 2
vec2 cayley(vec2 z) {
  z -= vec2(0,1);
  z *= 2.0/dot(z,z);
  z += vec2(0,1);
  return z;
}

float segment(vec2 p, vec2 a, vec2 b) {
  p -= a; b -= a;
  float h = dot(p,b)/dot(b,b);
  h = clamp(h, 0.0, 1.0);
  return length(p-b*h);
}

float line(vec2 p, vec2 a, vec2 b) {
  if (a == b) return 1e8;
  p -= a; b -= a;
  float h = dot(p,b)/dot(b,b);
  return length(p-b*h);
}

// color = point(color,z0,p,vec3(1,0,0));
vec3 point(vec3 color, vec2 z0, vec2 p, vec3 pcolor) {
  color *= 0.4+0.6*smoothstep(0.03,0.04,abs(distance(z0,p)));
  color = mix(pcolor,color,0.4+0.6*smoothstep(0.02,0.03,abs(distance(z0,p))));
  return color;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 z0 = 1.1*map(fragCoord.xy);
  vec2 mouse = vec2(1,0.5); // Identity
  if (iMouse.x > 0.0) mouse = map(iMouse.xy);
  float a = 0.0;
  a = 0.2*iTime;
  float b = PI*mouse.x;
  float c = PI*mouse.y;
  vec2 p = circle(a);
  vec2 q = circle(b);
  vec2 r = circle(c);
  // p1 etc. are 1d homogeneous
  vec2 p1 = stereographic(p);
  vec2 q1 = stereographic(q);
  vec2 r1 = stereographic(r);
  // h is homography on (real) line: p1,q1,r1,=> (1,0),(0,1),(1,1)
  mat2 h = rproject(p1,q1,r1);
  // The fourth point, parameter -1
  vec2 s1 = inverse(h)*vec2(1,-1); // On the line
  vec2 s = istereographic(s1);     // On the circle

  // Now p,q,r,s are base points for collineation
  // m takes p,q,r,s to unit points
  mat3 m = rproject(vec3(p,1),vec3(q,1),vec3(r,1),vec3(s,1));

  // n takes compass points to unit points
  vec3 p0 = vec3(0,1,1);  // = istereographic(vec2(1,0));
  vec3 q0 = vec3(0,-1,1); // = istereographic(vec2(0,1));
  vec3 r0 = vec3(1,0,1);  // = istereographic(vec2(1,1));
  vec3 s0 = vec3(-1,0,1); // = istereographic(vec2(1,-1));
  mat3 n = rproject(p0,q0,r0,s0);

  // combine
  m = inverse(n)*m;
  mat3 minv = inverse(m);

  vec3 color = vec3(0);
  int k = 0;
  if (!key(CHAR_X)) {
    vec3 p = vec3(z0,1);
    p = m*p;
    color += texture(iChannel0,norm(p)).xyz;
    k++;
  }
  if (!key(CHAR_Y)) {
    vec2 z = z0;
    // Map to half plane
    z = cayley(z);
    //assert(z.y > 0.0);
    {
      float a = h[0][0];
      float b = h[1][0];
      float c = h[0][1];
      float d = h[1][1];
      // Negative discriminant means inside/outside are flipped.
      //assert(a*d-b*c > 0.0);
      z = cdiv(a*z+vec2(b,0),c*z+vec2(d,0));
    }
    z = cayley(z);
    color += 0.5*texture(iChannel0,z).xyz;
    k++;
  }
  if (k == 0) color = vec3(1,1,0.5);
  else if (k > 1) color /= float(k);
  // z0 is untransformed point
  color *= 0.4+0.6*smoothstep(0.005,0.02,abs(length(z0)-1.0));
  color *= 0.4+0.6*smoothstep(0.005,0.02,line(z0,vec2(0,1),p));
  color *= 0.4+0.6*smoothstep(0.005,0.02,line(z0,vec2(0,1),q));
  color *= 0.4+0.6*smoothstep(0.005,0.02,line(z0,vec2(0,1),r));
  color *= 0.4+0.6*smoothstep(0.005,0.02,line(z0,vec2(0,1),s));
  color *= 0.4+0.6*smoothstep(0.005,0.02,abs(z0.y));

  color = point(color,z0,vec2(0,1),vec3(1));
  color = point(color,z0,p,vec3(1,0,0));
  color = point(color,z0,vec2(p1.x/p1.y,0),vec3(1,0,0));
  color = point(color,z0,q,vec3(0,1,0));
  color = point(color,z0,vec2(q1.x/q1.y,0),vec3(0,1,0));
  color = point(color,z0,r,vec3(0,0,1));
  color = point(color,z0,vec2(r1.x/r1.y,0),vec3(0,0,1));
  color = point(color,z0,s,vec3(1,1,0));
  color = point(color,z0,vec2(s1.x/s1.y,0),vec3(1,1,0));

  if (alert) color.r = 1.0;
  fragColor = vec4(pow(1.0*color,vec3(0.4545)),1);
}
