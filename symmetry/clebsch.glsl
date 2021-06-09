#version 300 es
precision mediump float;
uniform float iTime;
uniform vec2 iResolution;
uniform vec4 iMouse;

out vec4 outColor;

void mainImage( out vec4 fragColor, in vec2 fragCoord );

void main() {
  mainImage(outColor, gl_FragCoord.xy);
}

////////////////////////////////////////////////////////////////////////////////

bool highlight = false;

void assert(bool t) {
  if (!t) highlight = true;
}

const float PI =  3.141592654;
const float TWOPI = 2.0 * PI;
const float PHI = 1.618033989;
const vec4 I = vec4(1,1,1,1);
const vec3 AXIS = vec3(1,1,1); // axis for rotation
const int NLINES = 27;
vec4 lines[NLINES*2];
int colors[NLINES];

// Lighting params
vec3 light;
float ambient;
float diffuse;

// Quaternion multiplication as a matrix.
// w coordinate is real element of quaternion
// This gives a Clifford translation in R4
mat4 qmat(vec4 q) {
  float x = q.x, y = q.y, z = q.z, t = q.w;
  return mat4( t,-z, y, x, 
               z, t,-x, y,
              -y, x, t, z,
              -x,-y,-z, t );
}

#if 0
int quadratic(float A, float B, float C, out vec3 x) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return 0;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x[0] = (-B-D)/(2.0*A);
   x[1] = C/(A*x[0]);
   return 2;
}

// Numerical Recipes algorithm for solving cubic equation
int cubic0(vec4 coeffs, out vec3 x) {
  float a = coeffs[0], b = coeffs[1];
  float c = coeffs[2], d = coeffs[3];
  if (a == 0.0) return quadratic(b,c,d,x);
  //if (d == 0.0) return quadratic(a,b,c,x); // Need 0 too.
  float tmp = a; a = b/tmp; b = c/tmp; c = d/tmp;
  // solve x^3 + ax^2 + bx + c = 0
  float Q = (a*a-3.0*b)/9.0;
  float R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  float R2 = R*R, Q3 = Q*Q*Q;
  if (R2 < Q3) {
    float X = clamp(R/sqrt(Q3),-1.0,1.0);
    float theta = acos(X);
    float S = sqrt(Q); // Q must be positive since 0 <= R2 < Q3
    x[0] = -2.0*S*cos(theta/3.0) - a/3.0;
    x[1] = -2.0*S*cos((theta+2.0*PI)/3.0) - a/3.0;
    x[2] = -2.0*S*cos((theta+4.0*PI)/3.0) - a/3.0;
    return 3;
  } else {
    float A = -sign(R)*pow(abs(R)+sqrt(R2-Q3),0.3333);
    float B = A == 0.0 ? 0.0 : Q/A;
    x[0] = (A+B) - a/3.0;
    return 1;
  }
}

int cubic(vec4 coeffs, out vec3 x) {
  // Some ill-conditioned coeffs can cause problems
  // The worst is fixed by solving for reciprocal
  bool flip = false;
  if (abs(coeffs[0]) < abs(coeffs[3])) {
    flip = true;
    coeffs.xyzw = coeffs.wzyx;
  }
  int nroots = cubic0(coeffs,x);
  if (flip) x = 1.0/x;
  return nroots;
}
#else
// The Lanczos method
float evalcubic(float x, float a, float b, float c, float d) {
  return ((x*a+b)*x+c)*x+d;
}

float evalquad(float x, float a, float b, float c) {
  return (x*a+b)*x+c;
}

// Solve a*x**2 + b*x + c == 0
bool quadratic(float a, float b, float c, out vec2 t) {
   float d = b*b - 4.0*a*c;
   if (d < 0.0) return false;
   d = sqrt(d);
   if (b < 0.0) d = -d;
   t.x = 0.5*(-b-d)/a;
   t.y = c/(a*t.x);
   if (t.x > t.y) t.xy = t.yx; // Sort results
   return true;
}

// Find real root of x**3 + a*x**2 + b*x + c
float cubic1(float a, float b, float c) {
  if (c == 0.0) return 0.0;
  float k = 1.0; // Multiplication factor for answer
  // First, ensure c is negative
  if (c > 0.0) {
    a = -a; c = -c; k = -k;
  }
  // Now c is negative, but may be very small,
  // in which case we return an approximation.
  // Note that we expect a positive root.
  if (c > -1e-6) {
    if (b > 1e-10) return -k*c/b;
    // If abs(b) is small then the x^2 factor becomes
    // significant and we can just return 0.
    if (abs(b) < 1e-4) return 0.0;
    // Otherwise look elsewhere
  }
  // Now substitute to make c = -1
  float r = pow(-c,1.0/3.0); // cube root
  // Is is worth polishing r?
  a /= r; b /= r*r; c = -1.0; k *= r;
  // Now bracket a root between 0 and 1
  // We may need to solve for 1/x
  bool reciprocate = evalcubic(1.0,1.0,a,b,c) < 0.0;
  if (reciprocate) {
    // Swap a and b and negate.
    float t = a; a = -b; b = -t;
  }
  vec2 res;
  float x = 0.0;
  // Chebyshev polynomial: |32x**3 - 48x**2 + 18x - 1| < 1
  if (quadratic(a+1.5, b-0.5625, c+0.03125,res)) {
    // Find root closest to unit interval,
    x = (abs(res[0]-0.5) < abs(res[1]-0.5))? res[0]: res[1];
  }
  float eps = 1e-6;
  for (int i = 0; i < 6; i++) {
    // 6 iterations should be enough for anything
    // 2 usually suffices
    float t = evalcubic(x,1.0,a,b,c)/evalquad(x,3.0,2.0*a,b);
    x -= t;
    if (abs(t) < eps) break;
  }
  // and we are done.
  if (reciprocate) return k/x;
  return k*x;
}

int cubic(vec4 coeffs, out vec3 res) {
  float a = coeffs[1]/coeffs[0];
  float b = coeffs[2]/coeffs[0];
  float c = coeffs[3]/coeffs[0];
  float x = cubic1(a,b,c);
  res.x = x;
  vec2 tmp;
  float A = 1.0;
  float B = a+x;
  //float C = b+x*(a+x); // Less accurate it seems
  float C = (x == 0.0) ? b: -c/x;
  if (!quadratic(A,B,C,tmp)) return 1;
  res.yz = tmp.xy;
  // Sort the roots, yz are already sorted
  if (res.x > res.y) res.xy = res.yx;
  if (res.y > res.z) res.yz = res.zy;
  return 3;
}
#endif

// Initialize array of lines & colors
// I could probably do this procedurally or just use
// a static array for GLSL 3.0
void initlines() {
  lines[0] =  vec4(1,-1,0,0); lines[1] =  vec4(0,0,1,0);
  lines[2] =  vec4(1,-1,0,0); lines[3] =  vec4(0,0,0,1);
  lines[4] =  vec4(1,0,-1,0); lines[5] =  vec4(0,1,0,0);
  lines[6] =  vec4(1,0,-1,0); lines[7] =  vec4(0,0,0,1);
  lines[8] =  vec4(1,0,0,-1); lines[9] =  vec4(0,1,0,0);
  lines[10] = vec4(1,0,0,-1); lines[11] = vec4(0,0,1,0);

  lines[12] = vec4(0,1,-1,0); lines[13] = vec4(1,0,0,0);
  lines[14] = vec4(0,1,-1,0); lines[15] = vec4(0,0,0,1);
  lines[16] = vec4(0,1,0,-1); lines[17] = vec4(1,0,0,0);
  lines[18] = vec4(0,1,0,-1); lines[19] = vec4(0,0,1,0);
  lines[20] = vec4(0,0,1,-1); lines[21] = vec4(1,0,0,0);
  lines[22] = vec4(0,0,1,-1); lines[23] = vec4(0,1,0,0);

  lines[24] = vec4(1,-1,0,0); lines[25] = vec4(0,0,1,-1);
  lines[26] = vec4(1,0,-1,0); lines[27] = vec4(0,1,0,-1);
  lines[28] = vec4(1,0,0,-1); lines[29] = vec4(0,1,-1,0);

  lines[30] = vec4(1,PHI,-1,0); lines[31] = vec4(PHI,1,0,-1);
  lines[32] = vec4(1,PHI,0,-1); lines[33] = vec4(PHI,1,-1,0);
  lines[34] = vec4(1,0,PHI,-1); lines[35] = vec4(PHI,-1,1,0);
  lines[36] = vec4(1,-1,PHI,0); lines[37] = vec4(PHI,0,1,-1);
  lines[38] = vec4(1,-1,0,PHI); lines[39] = vec4(PHI,0,-1,1);
  lines[40] = vec4(1,0,-1,PHI); lines[41] = vec4(PHI,-1,0,1);

  lines[42] = vec4(-1,1,PHI,0); lines[43] = vec4(0,PHI,1,-1);
  lines[44] = vec4(0,1,PHI,-1); lines[45] = vec4(-1,PHI,1,0);
  lines[48] = vec4(-1,1,0,PHI); lines[49] = vec4(0,PHI,-1,1);
  lines[46] = vec4(0,1,-1,PHI); lines[47] = vec4(-1,PHI,0,1);
  lines[50] = vec4(-1,0,1,PHI); lines[51] = vec4(0,-1,PHI,1);
  lines[52] = vec4(0,-1,1,PHI); lines[53] = vec4(-1,0,PHI,1);

  colors[0] = colors[5] = colors[9] = 5;
  colors[1] = colors[3] = colors[7] = 2;
  colors[2] = colors[4] = colors[11] = 7;
  colors[6] = colors[8] = colors[10] = 4;
  colors[12] = colors[13] = colors[14] = 6;
  // This colouring shows a "double six" for 12 lines.
  colors[15] = colors[17] = colors[19] = colors[21] = colors[23] = colors[25] = 1;
  colors[16] = colors[18] = colors[20] = colors[22] = colors[24] = colors[26] = 3;
}

vec4 getColor(int i) {
  if (i == 0) return vec4(0.5,0.5,1,1);
  if (i == 1) return vec4(1,0,0,1);
  if (i == 2) return vec4(1,1,0,1);
  if (i == 3) return vec4(0,0,1,1);
  if (i == 4) return vec4(0,1,0,1);
  if (i == 5) return vec4(0,1,1,1);
  if (i == 6) return vec4(1,0,1,1);
  if (i == 7) return vec4(0.5,0.5,0.1,1);
  return vec4(0,0,0,1);
}

// p is in world coordinates
float eval(vec4 p,mat4 m) {
  p = m*p;
  float pI = dot(p,I);
  return dot(p*p,p) - pI*pI*pI;
}

// Equation is: pp.p - (p.I)^3 = 0
// A = rr.r - (r.I)^3
// B = 3(rr.p -(p.I)(r.I)^2)
// C = 3(pp.r -(p.I)^2(r.I))
// D = pp.p - (p.I)^3
float surface(vec4 p0, vec4 r0) {
  vec4 p = p0;
  vec4 r = r0;
  float pI = dot(p,I);
  float rI = dot(r,I);

  float A = dot(r*r - rI*rI*I, r);
  float B = 3.0*dot(r*r - rI*rI*I, p);
  float C = 3.0*dot(p*p - pI*pI*I, r);
  float D = dot(p*p - pI*pI*I, p);

  vec3 res;
  int nroots = cubic(vec4(A,B,C,D),res);
  // return smallest positive root
  for (int i = 0; i < nroots; i++) {
    if (res[i] > 0.0) return res[i];
  }
  return 1e8;
}
    
// Find the (squared) distance to the line in R3 but
// using homogeneous coordinates.
float tryline(vec4 q, vec4 p, vec4 r) {
  // |r.w| < |p.w| and if both are small we have a
  // line at infinity so we can ignore it
  if (abs(p.w) < 1e-2) return 1e8;
   r = p.w*r - r.w*p; // r.w = 0
   q = p.w*q-p;
   float qr = dot(q,r);
   float d2 = (dot(q,q) - (qr*qr/dot(r,r)))/(p.w*p.w);
   return d2;
}

// Multiply 2 quaternions
vec4 qmul(vec4 p, vec4 q) {
  // xyz is ijk, w is real part
  vec3 P = p.xyz, Q = q.xyz;
  return vec4(cross(P,Q)+p.w*Q+q.w*P,p.w*q.w-dot(P,Q));
}
  
vec4 solve(vec4 p, vec4 r) {
  vec3 normal = vec3(1,0,0);
  float theta = iTime*0.1;
  // M rotates surface in R4
  vec3 a = normalize(AXIS);
  mat4 M = qmat(vec4(sin(theta)*a,cos(theta)));
  mat4 Minv = qmat(vec4(-sin(theta)*a,cos(theta)));
  initlines();

  const float u = 0.025;
  const float u2 = u*u;
  int colorindex = -1; // Selected color
  float t = surface(M*p,M*r);
  if (t < 1e8) {
    // Compute normal by finding gradient
    vec4 q = p+t*r;
    float x0 = eval(q,M);
    vec2 eps = vec2(1e-2,0);
    normal = vec3 (eval(q+eps.xyyy,M)-x0,
                   eval(q+eps.yxyy,M)-x0,
                   eval(q+eps.yyxy,M)-x0);
    normal = normalize(normal);
    colorindex = 0;

    // Now find the nearest line
    float d = u2;
    for (int i = 0; i < NLINES; i++) {
      vec4 p4 = Minv*lines[2*i];
      vec4 r4 = Minv*lines[2*i+1];
      if (abs(p4.w) < abs(r4.w)) {
        vec4 t = p4; p4 = r4; r4 = t;
      }
      float d0 = tryline(q,p4,r4);
      if (d0 < d){
        d = d0;
        colorindex = colors[i];
      }
    }
  }
  if (colorindex == -1) return vec4(0,0,0,1);
  if (dot(r.xyz,normal) > 0.0) {
    normal = -normal; // Face forwards
  }
  vec4 baseColor = getColor(colorindex);
  vec3 color;
  color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,normal));
  float specular = pow(max(0.0,dot(reflect(light,normal),r.xyz)),4.0);
  color += 0.7*specular*vec3(1.0,1.0,1.0);
  return vec4(sqrt(color),1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  light = normalize(vec3(0.0,1.0,-1.0));
  ambient = 0.6;
  diffuse = 1.0-ambient;

  vec2 uv = 2.0*fragCoord.xy/iResolution.xy - 1.0;
  vec3 p = vec3(0.0, 0.0, 6.0);
  vec3 r = normalize(vec3(iResolution.x/iResolution.y * uv.x, uv.y, -3.0));
  // Rotate camera according to mouse position
  float xrot = 0.0; // About x-axis
  float yrot = 0.0; // About y-axis
  if (iMouse.x > 0.0) {
    yrot = TWOPI*(iMouse.x-0.5*iResolution.x)/iResolution.x;
    xrot = TWOPI*(iMouse.y-0.5*iResolution.y)/iResolution.y;
  }
  mat3 mx = mat3(1,0,0,
                 0,cos(xrot),-sin(xrot),
                 0,sin(xrot),cos(xrot));
  mat3 my = mat3(cos(yrot),0,-sin(yrot),
                 0,1,0,
                 sin(yrot),0,cos(yrot));
  mat3 m = mx*my;
  light *= m;
  p *= m;
  r *= m;
  fragColor = solve(vec4(p,1),vec4(r,0));
  if (highlight) fragColor = vec4(1,0,0,1);
}
