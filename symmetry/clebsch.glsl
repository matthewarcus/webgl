precision highp float;
uniform float iGlobalTime;
uniform vec2 iResolution;

////////////////////////////////////////////////////////////////////////////////

const float PI =  3.141592654;
const float PHI = 1.618033989;
const vec4 I = vec4(1,1,1,1);
const vec3 AXIS = vec3(1,1,1); // axis for rotation
//const int NLINES = 27;
const int NLINES = 12;
vec4 lines[NLINES*2];
int colors[NLINES];

// Rotation & inverse in R4
mat4 m4;
mat4 m4inv;

// Lighting params
vec3 light;
float ambient;
float diffuse;

float d2best = 1e10; // Smallest distance to line
int colorindex = -1; // Selected color

// Quaternion multiplication as a matrix.
// w coordinate is real element of quaternion
mat4 qmat(vec4 q) {
  float x = q.x, y = q.y, z = q.z, t = q.w;
  return mat4( t,-z, y, x, 
               z, t,-x, y,
              -y, x, t, z,
              -x,-y,-z, t );
}

int quadratic(float A, float B, float C, out vec3 x) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return 0;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x[0] = (-B-D)/(2.0*A);
   x[1] = C/(A*x[0]);
   return 2;
}

// NR algorithm for solving cubic equation
int cubic0(float a, float b, float c, float d, out vec3 x) {
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

int cubic(float A, float B, float C, float D, out vec3 x) {
  int nroots;
  // Some ill-conditioned coeffs can cause problems
  // The worst is fixed by solving for reciprocal
  if (abs(A) > abs(D)) {
    nroots = cubic0(A,B,C,D,x);
  } else {
    nroots = cubic0(D,C,B,A,x);
    for (int i = 0; i < 3; i++) {
      x[i] = 1.0/x[i];
    }
  }
  return nroots;
}

// Initialize array of lines & colors
void initlines() {
  lines[0] =  vec4(1,-1,0,0); lines[1] =  vec4(0,0,1,0);
  lines[2] =  vec4(1,-1,0,0); lines[3] =  vec4(0,0,0,1);
  lines[4] =  vec4(1,0,-1,0); lines[5] =  vec4(0,1,0,0);
  lines[6] =  vec4(1,0,-1,0); lines[7] =  vec4(0,0,0,1);
  lines[8] =  vec4(1,0,0,-1); lines[9] =  vec4(0,1,0,0);
  lines[10] = vec4(1,0,0,-1); lines[11] = vec4(0,0,1,0);

#if 0
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
#endif

  colors[0] = colors[5] = colors[9] = 5;
  colors[1] = colors[3] = colors[7] = 2;
  colors[2] = colors[4] = colors[11] = 7;
  colors[6] = colors[8] = colors[10] = 4;
#if 0
  colors[12] = colors[13] = colors[14] = 6;
  // This colouring shows a "double six" for 12 lines.
  colors[15] = colors[17] = colors[19] = colors[21] = colors[23] = colors[25] = 1;
  colors[16] = colors[18] = colors[20] = colors[22] = colors[24] = colors[26] = 3;
#endif
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
float eval(vec4 p) {
  p = m4*p;
  float pI = dot(p,I);
  return dot(p*p,p) - pI*pI*pI;
}

// Equation is: pp.p - (p.I)^3 = 0
// A = rr.r - (r.I)^3
// B = 3(rr.p -(p.I)(r.I)^2)
// C = 3(pp.r -(p.I)^2(r.I))
// D = pp.p - (p.I)^3
bool surface(vec4 p0, vec4 r0, float min, inout float t0, out vec3 normal) {
  vec4 p = m4*p0;
  vec4 r = m4*r0;
  float pI = dot(p,I);
  float rI = dot(r,I);
#if 0
  float A = dot(r*r,r) - rI*rI*rI;
  float B = 3.0*(dot(r*r,p) - pI*rI*rI);
  float C = 3.0*(dot(p*p,r) - pI*pI*rI);
  float D = dot(p*p,p)-pI*pI*pI;
#else
  // Possibly this is a better way.
  float A = dot(r*r - rI*rI*I, r);
  float B = 3.0*dot(r*r - rI*rI*I, p);
  float C = 3.0*dot(p*p - pI*pI*I, r);
  float D = dot(p*p - pI*pI*I, p);
#endif
  vec3 x;
  int nroots = cubic(A,B,C,D,x);
  bool found = false;
  for (int i = 0; i < 3; i++) {
    if (i == nroots) break;
    if (min < x[i] && x[i] < t0) {
      t0 = x[i];
      found = true;
    }
  }
  // Compute normal by finding gradient
  vec4 p1 = p0+t0*r0;
  float x0 = eval(p1);
  vec2 eps = vec2(1e-2,0);
  normal = vec3 (eval(p1+eps.xyyy)-x0,
                 eval(p1+eps.yxyy)-x0,
                 eval(p1+eps.yyxy)-x0);
  normal = normalize(normal);
  return found;
}
    
// Find the distance to the line in R3
void tryline(vec4 q, vec4 p, vec4 r, int color) {
  const float u = 0.1; //0.025;
  const float u2 = u*u;
  if (abs(p.w) < 1e-2) return;
   r = p.w*r - r.w*p; // r.w = 0
   q = p.w*q-p;
   float qr = dot(q,r);
   float d2 = (dot(q,q) - (qr*qr/dot(r,r)))/(p.w*p.w);
   if (d2 < min(u2,d2best)) {
      d2best = d2;
      colorindex = color;
   }
}

vec4 solve(vec3 p, vec3 r, float min) {
  vec3 normal = vec3(1,0,0);
  float theta = iGlobalTime*0.1; 
  vec3 a = normalize(AXIS);
  m4 = qmat(vec4(sin(theta)*a,cos(theta)));
  m4inv = qmat(vec4(sin(-theta)*a,cos(-theta)));
  initlines();

  float t = 1e10;
  if (surface(vec4(p,1),vec4(r,0),min,t,normal)) {
    colorindex = 0;
    vec4 q = vec4(p+t*r,1);
    for (int i = 0; i < NLINES; i++) {
      vec4 p = m4inv*lines[2*i];
      vec4 r = m4inv*lines[2*i+1];
      if (abs(p.w) < abs(r.w)) {
        tryline(q,r,p,colors[i]);
      } else {
        tryline(q,p,r,colors[i]);
      }
    }
  }
  if (colorindex == -1) return vec4(0,0,0,1);
  if (dot(r,normal) > 0.0) {
    normal = -normal; // Face forwards
  }
  vec4 baseColor = getColor(colorindex);
  vec3 color;
  color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,normal));
  float specular = pow(max(0.0,dot(reflect(light,normal),r)),4.0);
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
  float k = -dot(p,r);
  p += k*r;
  fragColor = solve(p,r,-k);
}

////////////////////////////////////////////////////////////////////////////////

void main() {
  mainImage(gl_FragColor, gl_FragCoord.xy);
}
