#version 300 es
//#extension GL_EXT_shader_texture_lod : enable
//#extension GL_OES_standard_derivatives : enable

precision highp float;

out vec4 outColor;

uniform vec4 uClock;
uniform vec4 params1;
uniform vec4 params2;
uniform vec4 params3;
uniform mat4 uMatrix;
uniform vec2 iResolution;
uniform float K0;
uniform int iTime;

const int NLINES = 27;
uniform int nlines;
uniform vec4 uLines[2*NLINES];
uniform int uColors[NLINES];
uniform samplerCube uCubeMap;

////////////////////////////////////////////////////////////////////////////////

const float PI =  3.141592654;
const float PHI = 1.618033989;
const vec4 I = vec4(1,1,1,1);
const float eps = 1e-3;
//const float K0 = 0.27;
//const float K0 = 2.0;
//const float K0 = 10.23;
//const float K0 = 0.5;
//const float K0 = 0.3;
//const vec4 K = vec4(1,1.0/4.0,1.0/9.0,1.0/25.0);
//const float K0 = eps+1.0/16.0;
//const float K0 = 0.5; //eps+0.25;
vec4 K;

float clock0;
float clock1;

bool highlight = false; // Debugging

mat4 m;
mat4 minv;

int colorindex = -1;

vec3 light;
float ambient;
float diffuse;

// Quaternion multiplication as a matrix.
// w coordinate is real element of quaternion
mat4 qmat(vec4 q) {
  float x = q.x, y = q.y, z = q.z, t = q.w;
  return mat4( t,-z, y, x, 
               z, t,-x, y,
              -y, x, t, z,
              -x,-y,-z, t );
}

void swap(inout float x1, inout float x2) {
    float t = x1; x1 = x2; x2 = t;
}

void sort(inout float x1, inout float x2) {
  if (x1 > x2) swap(x1,x2);
}

void sort(inout float x1, inout float x2, inout float x3) {
  sort(x1,x2); sort(x2,x3); sort(x1,x2);
}

// Use vec2 return?
bool quadratic0(float A, float B, float C, out float x1, out float x2) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return false;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x1 = (-B-D)/(2.0*A);
   x2 = C/(A*x1);
   sort(x1,x2);
   return true;
}

int quadratic(float A, float B, float C, out vec3 x) {
  float x0,x1;
  if (!quadratic0(A,B,C,x0,x1)) return 0;
  x[0] = x0; x[1] = x1;
  return 2;
}

int cubic0(float a, float b, float c, float d, out vec3 x) {
#if 0
  if (d >= 0.0 && abs(d/a) < 1e-4) d = 1e-4*abs(a);
  else if (d < 0.0 && abs(d/a) < 1e-4) d = -1e-4*abs(a);
  if (abs(d/a) < 1e-4) {
    //highlight = true;
  }
  if (a == 0.0) {
    //highlight = true;
    return quadratic(b,c,d,x);
  }
  if (d == 0.0) {
    //highlight = true;
    return quadratic(a,b,c,x);
  }
#endif
#if 0
  //if (abs(c) < 1e-3) highlight = true;
  if (abs(a) < 1e-3) {
    //highlight = true;
    // This helps, but isn't perfect.
    float A = 1.0;
    float C = d/b;
    float B = (c-a*C)/b;
    return quadratic(A,B,C,x);
  }
#endif
#if 0
  int nroots = 1;
  {
    // Find out how many roots to expect by looking at
    // the derivative. We have three roots just when
    // there are 2 different points where the derivative
    // goes to zero.
    float x0, x1;
    if (quadratic0(3.0*a,2.0*b,c,x0,x1) && x0 != x1) {
      float y0 = evalcubic(x0,a,b,c,d);
      float y1 = evalcubic(x1,a,b,c,d);
      if (y0*y1 < 0.0) nroots = 3;
    }
  }
#endif
  float tmp = a; a = b/tmp; b = c/tmp; c = d/tmp;
  // solve x^3 + ax^2 + bx + c = 0
  float Q = (a*a-3.0*b)/9.0;
  float R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  float R2 = R*R, Q3 = Q*Q*Q;
  if (R2 < Q3) {
    //if (nroots != 3) highlight = true;
    float X = clamp(R/sqrt(Q3),-1.0,1.0);
    float theta = acos(X);
    float S = sqrt(Q); // Q must be positive since 0 <= R2 < Q3
    x[0] = -2.0*S*cos(theta/3.0) - a/3.0;
    x[1] = -2.0*S*cos((theta+2.0*PI)/3.0) - a/3.0;
    x[2] = -2.0*S*cos((theta+4.0*PI)/3.0) - a/3.0;
    return 3;
  } else {
    //if (nroots != 1) highlight = true;
    float A = -sign(R)*pow(abs(R)+sqrt(R2-Q3),0.3333);
    float B = A == 0.0 ? 0.0 : Q/A;
    x[0] = (A+B) - a/3.0;
    return 1;
  }
}

int cubic(float A, float B, float C, float D, out vec3 x) {
  int nroots;
  // Some ill-conditioned coeffs can cause problems
  // The worst is fixed by solving for reciprocal.
  // Jim Blinn suggests a more sophisticated method
  // solving both ways & selecting best roots from each.
  if (abs(A) > abs(D)) {
    nroots = cubic0(A,B,C,D,x);
  } else {
    nroots = cubic0(D,C,B,A,x);
    for (int i = 0; i < 3; i++) {
      if (i == nroots) break;
      x[i] = 1.0/x[i];
    }
  }
  if (nroots == 2) sort(x[0],x[1]);
  else if (nroots == 3) sort(x[0],x[1],x[2]);
  return nroots;
}

vec4 getColor(int i) {
  if (i == 0) return vec4(0.5,0.5,1,1);
  if (i == 1) return vec4(1,0,0,1);
  if (i == 2) return vec4(0,0,1,1);
  if (i == 3) return vec4(1,1,0,1);
  if (i == 4) return vec4(0,1,0,1);
  if (i == 5) return vec4(0,1,1,1);
  if (i == 6) return vec4(1,0,1,1);
  if (i == 7) return vec4(0.5,0.5,0.1,1);
  return vec4(0,0,0,1);
}

// p is in world coordinates
float eval(vec4 p) {
  p = m*p;
#if 1
  return p.x*p.w - p.y*p.z;
#else
  float pI = dot(p,I);
  return dot(K*p*p,p) - pI*pI*pI;
#endif
}

// Equation is: Kpp.p - (p.I)^3 = 0
// A = rr.r - (r.I)^3
// B = 3(rr.p -(p.I)(r.I)^2)
// C = 3(pp.r -(p.I)^2(r.I))
// D = pp.p - (p.I)^3
int surface(vec4 p0, vec4 r0, inout vec3 x) {
  vec4 p = m*p0;
  vec4 r = m*r0;
#if 0
  float A = r.x*r.w - r.y*r.z;
  float B = p.x*r.w + r.x*p.w - p.y*r.z - r.y*p.z;
  float C = p.x*p.w - p.y*p.z;
  int nroots = quadratic(A,B,C,x);
#else
  float pI = dot(p,I);
  float rI = dot(r,I);
#if 0
  float A = dot(K*r*r,r) - rI*rI*rI;
  float B = 3.0*(dot(K*r*r,p) - pI*rI*rI);
  float C = 3.0*(dot(K*p*p,r) - pI*pI*rI);
  float D = dot(K*p*p,p)-pI*pI*pI;
#else
  // Possibly this is a better way.
  float A = dot(K*r*r - rI*rI*I, r);
  float B = 3.0*dot(K*r*r - rI*rI*I, p);
  float C = 3.0*dot(K*p*p - pI*pI*I, r);
  float D = dot(K*p*p - pI*pI*I, p);
#endif
  int nroots = cubic(A,B,C,D,x);
#endif
  return nroots;
}
    
// Find the distance to the line in R3
void tryline(vec4 q, vec4 p, vec4 r, int color, inout float d2best) {
  const float u = 0.025;
  const float u2 = u*u;
  if (abs(p.w) < 1e-2) return;
   //p /= p.w; // p.w = 1
   r = p.w*r - r.w*p; // r.w = 0
   //r = normalize(r); // r.r = 1
   q = p.w*q-p;
   float qr = dot(q,r); // = dot(p.w*q0-p,r)
   float d2 = (dot(q,q) - (qr*qr/dot(r,r)))/(p.w*p.w);
   if (d2 < min(u2,d2best)) {
      d2best = d2;
      colorindex = color;
   }
}

vec4 solve(vec4 p, vec4 r, float tmin) {
  float theta = clock1*0.1; 
#if 1
  mat4 minv0 = 0.25*mat4(-1,1,1,1,
                         1,-1,1,1,
                         1,1,-1,1,
                         -1,-1,-1,1);
  mat4 m0 = 0.25*mat4(-1,1,1,-1,
                      1,-1,1,-1,
                      1,1,-1,-1,
                      1,1,1,1);
#else
  mat4 minv0 = mat4(0.5,4,-0.8660253882408142,-3,0.8660253882408142,4,0.5,-3,-1,4,0,-3,0,0,0,3);
  mat4 m0 = mat4(0.21132487058639526,0.3660253882408142,-0.5773502588272095,-0,0.052831217646598816,0.09150634706020355,0.10566243529319763,0.2499999850988388,-0.7886751294136047,0.6339746117591858,0.15470051765441895,-0,-0,-0,-0,0.3333333432674408);
#endif

#if 0
  m = mat4(1.0, 0.0, 0.0, 0.0,
                0.0, cos(theta), 0.0, -sin(theta),
                0.0, 0.0, 1.0, 0.0,
                0.0, sin(theta), 0.0, cos(theta));
  minv = mat4(1.0, 0.0, 0.0, 0.0,
                0.0, cos(theta), 0.0, sin(theta),
                0.0, 0.0, 1.0, 0.0,
                0.0, -sin(theta), 0.0, cos(theta));
#elif 0
  m = mat4(cos(theta), 0.0, -sin(theta), 0.0,
           0.0, 1.0, 0.0, 0.0,
           sin(theta), 0.0, cos(theta), 0.0,
                0.0, 0.0, 0.0, 1.0);
  minv = mat4(cos(theta), 0.0, sin(theta), 0.0,
                0.0, 1.0, 0.0, 0.0,
                -sin(theta), 0.0, cos(theta), 0.0,
                0.0, 0.0, 0.0, 1.0);
#else
  vec3 a = normalize(vec3(1,1,1));
  //vec3 a = normalize(vec3(0,0,1));
  m = qmat(vec4(sin(theta)*a,cos(theta)));
  minv = qmat(vec4(sin(-theta)*a,cos(-theta)));
#endif
  //m = m0*m; minv = minv*minv0;

  //initlines();

  bool dolines = false;
  dolines = true;
  vec3 x = vec3(0);
  int nroots = surface(p,r,x);
  //sort(x[0],x[1],x[2]);
  vec3 color = vec3(0);
  color = vec3(0);
  //color = textureCube(uCubeMap,r.xyz).rgb;
  for (int root = 0; root < 3; root++) {
    if (root == nroots) break;
    float t = x[root];
    //if (t < 0.0) highlight = true;
    if (t > tmin) {
      vec4 p1 = p+t*r;
      float x0 = eval(p1);
      vec2 eps = vec2(1e-2,0);
      vec3 normal = vec3 (eval(p1+eps.xyyy)-x0,
                     eval(p1+eps.yxyy)-x0,
                     eval(p1+eps.yyxy)-x0);
      normal = normalize(normal);
      colorindex = 0;
      if (dolines) {
        vec4 q = p+t*r;
        float d2best = 1e10;
        for (int i = 0; i < 27; i++) {
          if (i == nlines) break;
          vec4 p = minv*uLines[2*i];
          vec4 r = minv*uLines[2*i+1];
          if (abs(p.w) < abs(r.w)) {
            tryline(q,r,p,uColors[i],d2best);
          } else {
            tryline(q,p,r,uColors[i],d2best);
          }
        }
      }
      if (dot(vec3(r),normal) > 0.0) {
        normal = -normal; // Face forwards
      }
      float mix = 0.5;
      //float fog = min(1.0/t,1.0);
      color = clamp(color,0.0,1.0);
      //color = 0.4*color + 0.1*fog*vec3(1,1,1);
      color = 0.4*color;
      vec3 baseColor = getColor(colorindex).xyz;
      color += mix*baseColor*ambient;
      float k = dot(light,normal);
      if (k > 0.0) {
        color += mix*baseColor*diffuse*k;
        float specular = pow(max(0.0,dot(reflect(light,normal),vec3(r))),10.0);
        color += 0.7*specular*vec3(1.0,1.0,1.0);
      }
      break;
    }
  }
  if (highlight) color.r = 1.0;
  return vec4(sqrt(color),1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  clock0 = uClock[0];
  clock1 = uClock[1];
  //K0 = 1.0 + sin(iTime*0.2);
  K = vec4(1,1,1,K0);
  float scale = params1[0];       // Width multiplier

  light = normalize(vec3(0.0,1.0,-1.0));
  ambient = 0.3;
  diffuse = 1.0-ambient;

  vec2 uv = scale*(fragCoord.xy/iResolution.xy - 0.5);
  float camera = 4.0;
  vec4 p = vec4(0.0, 0.0, -camera,1.0);
  vec4 r = normalize(vec4(iResolution.x/iResolution.y * uv.x, uv.y, 2.0, 0.0));
  float k = -dot(p,r);
  p += k*r;
  mat4 rmat = uMatrix;
  light = mat3(rmat) * light;
  p = rmat * p;
  r = rmat * r;
  light.z = -light.z; p.z = -p.z; r.z = -r.z;
  fragColor = solve(p,r,-k);
}

////////////////////////////////////////////////////////////////////////////////

void main() {
  mainImage(outColor, gl_FragCoord.xy);
}
