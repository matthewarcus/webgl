//#extension GL_EXT_shader_texture_lod : enable
//#extension GL_OES_standard_derivatives : enable

precision highp float;
uniform vec4 uClock;
uniform vec2 iResolution;

const int NLINES = 27;
uniform float K0;
uniform int nlines;
uniform vec4 uLines[2*NLINES];
uniform int uColors[NLINES];

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

#if 0
vec4 lines[NLINES*2];
int colors[NLINES];
int nlines = 3;
#endif

float d2best = 1e10;
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
  if (a == 0.0) {
    //highlight = true;
    return quadratic(b,c,d,x);
  }
  if (d == 0.0) {
    //highlight = true;
    return quadratic(a,b,c,x);
  }
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
  // The worst is fixed by solving for reciprocal
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

#if 0
void initlines() {
  if (K0 != 1.0) {
    // Lines on the base
    lines[0] = minv*vec4(1,0,0,0); lines[1] = minv*vec4(0,1,-1,0);
    lines[2] = minv*vec4(0,1,0,0); lines[3] = minv*vec4(1,0,-1,0);
    lines[4] = minv*vec4(0,0,1,0); lines[5] = minv*vec4(1,-1,0,0);
    if (false && K0 == 0.25) {
      nlines = 9;
      // Lines from "central" vertex of tetrahedron to "outer" vertices
      lines[6] = minv*vec4(-1,-1,-1,2); lines[7] = minv*vec4(1,-1,-1,2);
      lines[8] = minv*vec4(-1,-1,-1,2); lines[9] = minv*vec4(-1,1,-1,2);
      lines[10] = minv*vec4(-1,-1,-1,2); lines[11] = minv*vec4(-1,-1,1,2);

      // Lines between "outer" vertices
      lines[12] = minv*vec4(1,-1,-1,2); lines[13] = minv*vec4(-1,1,-1,2);
      lines[14] = minv*vec4(-1,1,-1,2); lines[15] = minv*vec4(-1,-1,1,2);
      lines[16] = minv*vec4(-1,-1,1,2); lines[17] = minv*vec4(1,-1,-1,2);
    } else {
      // For K0 = 0.25 get z0 = z1 = -0.5
      // For K0 = 0.5, get C = 0.5, (-3 +- sqrt(9 - 6))/6 = -1/2(1 +- 1/sqrt(3));
      // For K0 = 1, get z0 = 0, z1 = -1
      float A = 3.0, B = 3.0, C = 1.0-K0, z0, z1;
      if (quadratic0(A,B,C,z0,z1)) {
        nlines += 12;
        lines[6] = minv*vec4(0,0,z0,1); lines[7] = minv*vec4(1,-1,0,0);
        lines[8] = minv*vec4(0,z0,0,1); lines[9] = minv*vec4(1,0,-1,0);
        lines[10] = minv*vec4(z0,0,0,1); lines[11] = minv*vec4(0,1,-1,0);
        lines[12] = minv*vec4(0,0,z1,1); lines[13] = minv*vec4(1,-1,0,0);
        lines[14] = minv*vec4(0,z1,0,1); lines[15] = minv*vec4(1,0,-1,0);
        lines[16] = minv*vec4(z1,0,0,1); lines[17] = minv*vec4(0,1,-1,0);

        {
          const int base = 18;//2*nlines;
          lines[base+0] = minv*vec4(z0,z1,0,1); lines[base+1] = minv*vec4(0,0,1,0);
          lines[base+2] = minv*vec4(z0,0,z1,1); lines[base+3] = minv*vec4(0,1,0,0);
          lines[base+4] = minv*vec4(0,z0,z1,1); lines[base+5] = minv*vec4(1,0,0,0);

          lines[base+6+0] = minv*vec4(z1,z0,0,1); lines[base+6+1] = minv*vec4(0,0,1,0);
          lines[base+6+2] = minv*vec4(z1,0,z0,1); lines[base+6+3] = minv*vec4(0,1,0,0);
          lines[base+6+4] = minv*vec4(0,z1,z0,1); lines[base+6+5] = minv*vec4(1,0,0,0);
        }
        float b = 0.607644, a = 0.444341;
        //float a = 0.5125,  b = 0.5;
        //float a = 1.85802,  b = 5.50585;
        //float a = 0.844, b = 2.371; // K0 = 2
        //float a = 0.475, b = 1.052; // K0 = 0.5
        //float a = 0.434, b = 0.692;
        if (true) {
          nlines += 6;
          const int base = 30;//2*nlines;
          lines[base+0] = minv*vec4(a, -a,z0,1); lines[base+1] = minv*vec4(b, z1,-b,1);
          lines[base+2] = minv*vec4(a, z0,-a,1); lines[base+3] = minv*vec4(b, -b,z1,1);
          lines[base+4] = minv*vec4(-a, a,z0,1); lines[base+5] = minv*vec4(z1, b,-b,1);
          lines[base+6] = minv*vec4(-a,z0, a,1); lines[base+7] = minv*vec4(z1,-b, b,1);
          lines[base+8] = minv*vec4(z0, a,-a,1); lines[base+9] = minv*vec4(-b, b,z1,1);
          lines[base+10] = minv*vec4(z0,-a, a,1); lines[base+11] = minv*vec4(-b,z1, b,1);
        }
        {
          nlines += 6;
          swap(a,b);
          a = -a; b = -b;
          const int base = 42;//2*nlines;
          lines[base+0] = minv*vec4(a, -a,z0,1); lines[base+1] = minv*vec4(b, z1,-b,1);
          lines[base+2] = minv*vec4(a, z0,-a,1); lines[base+3] = minv*vec4(b, -b,z1,1);
          lines[base+4] = minv*vec4(-a, a,z0,1); lines[base+5] = minv*vec4(z1, b,-b,1);
          lines[base+6] = minv*vec4(-a,z0, a,1); lines[base+7] = minv*vec4(z1,-b, b,1);
          lines[base+8] = minv*vec4(z0, a,-a,1); lines[base+9] = minv*vec4(-b, b,z1,1);
          lines[base+10] = minv*vec4(z0,-a, a,1); lines[base+11] = minv*vec4(-b,z1, b,1);
        }
      }
    }
    colors[0] = colors[1] = colors[2] = 0;
    colors[3] = colors[4] = colors[5] = 1;
    colors[6] = colors[7] = colors[8] = 2;
    colors[9] = colors[10] = colors[11] = 3;
    colors[12] = colors[13] = colors[14] = 4;
    colors[15] = colors[16] = colors[17] = 5;
    colors[18] = colors[19] = colors[20] = 5;
    colors[21] = colors[22] = colors[23] = 6;
    colors[24] = colors[25] = colors[26] = 6;
  } else {
    nlines = 27;
    lines[0] = minv*vec4(1,-1,0,0); lines[1] = minv*vec4(0,0,1,0);
    lines[2] = minv*vec4(1,-1,0,0); lines[3] = minv*vec4(0,0,0,1);
    lines[4] = minv*vec4(1,0,-1,0); lines[5] = minv*vec4(0,1,0,0);
    lines[6] = minv*vec4(1,0,-1,0); lines[7] = minv*vec4(0,0,0,1);
    lines[8] = minv*vec4(1,0,0,-1); lines[9] = minv*vec4(0,1,0,0);
    lines[10] = minv*vec4(1,0,0,-1); lines[11] = minv*vec4(0,0,1,0);

    lines[12] = minv*vec4(0,1,-1,0); lines[13] = minv*vec4(1,0,0,0);
    lines[14] = minv*vec4(0,1,-1,0); lines[15] = minv*vec4(0,0,0,1);
    lines[16] = minv*vec4(0,1,0,-1); lines[17] = minv*vec4(1,0,0,0);
    lines[18] = minv*vec4(0,1,0,-1); lines[19] = minv*vec4(0,0,1,0);
    lines[20] = minv*vec4(0,0,1,-1); lines[21] = minv*vec4(1,0,0,0);
    lines[22] = minv*vec4(0,0,1,-1); lines[23] = minv*vec4(0,1,0,0);

    lines[24] = minv*vec4(1,-1,0,0); lines[25] = minv*vec4(0,0,1,-1);
    lines[26] = minv*vec4(1,0,-1,0); lines[27] = minv*vec4(0,1,0,-1);
    lines[28] = minv*vec4(1,0,0,-1); lines[29] = minv*vec4(0,1,-1,0);

    lines[30] = minv*vec4(1,PHI,-1,0); lines[31] = minv*vec4(PHI,1,0,-1);
    lines[32] = minv*vec4(1,PHI,0,-1); lines[33] = minv*vec4(PHI,1,-1,0);
    lines[34] = minv*vec4(1,-1,PHI,0); lines[35] = minv*vec4(PHI,0,1,-1);
    lines[36] = minv*vec4(1,0,PHI,-1); lines[37] = minv*vec4(PHI,-1,1,0);
    lines[38] = minv*vec4(1,-1,0,PHI); lines[39] = minv*vec4(PHI,0,-1,1);
    lines[40] = minv*vec4(1,0,-1,PHI); lines[41] = minv*vec4(PHI,-1,0,1);

    lines[42] = minv*vec4(-1,1,PHI,0); lines[43] = minv*vec4(0,PHI,1,-1);
    lines[44] = minv*vec4(0,1,PHI,-1); lines[45] = minv*vec4(-1,PHI,1,0);
    lines[46] = minv*vec4(-1,1,0,PHI); lines[47] = minv*vec4(0,PHI,-1,1);
    lines[48] = minv*vec4(0,1,-1,PHI); lines[49] = minv*vec4(-1,PHI,0,1);
    lines[50] = minv*vec4(-1,0,1,PHI); lines[51] = minv*vec4(0,-1,PHI,1);
    lines[52] = minv*vec4(0,-1,1,PHI); lines[53] = minv*vec4(-1,0,PHI,1);

    colors[0] = colors[2] = colors[4] = colors[6] = colors[8] = colors[10] = 0;
    colors[1] = colors[3] = colors[5] = colors[7] = colors[9] = colors[11] = 1;
    colors[12] = colors[13] = colors[14] = 3;
    colors[15] = colors[17] = colors[19] = colors[21] = colors[23] = colors[25] = 2;
    colors[16] = colors[18] = colors[20] = colors[22] = colors[24] = colors[26] = 4;
  }
}
#endif

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
  float pI = dot(p,I);
  return dot(K*p*p,p) - pI*pI*pI;
}

// Equation is: Kpp.p - (p.I)^3 = 0
// A = rr.r - (r.I)^3
// B = 3(rr.p -(p.I)(r.I)^2)
// C = 3(pp.r -(p.I)^2(r.I))
// D = pp.p - (p.I)^3
bool surface(vec4 p0, vec4 r0, float min, inout float t0, out vec3 normal) {
  vec4 p = m*p0;
  vec4 r = m*r0;
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
  float t = t0;
  vec3 x;
  int nroots = cubic(A,B,C,D,x);
  bool found = false;
  for (int i = 0; i < 3; i++) {
    if (i == nroots) break;
    float t = x[i];
    if (min < t && t < t0) {
#if 0
      // Insert test for point inclusion here
      vec4 p1 = 10.0*(p+t*r);
      float k = 0.5;
      if (p1.x-floor(p1.x) < k &&
          p1.z-floor(p1.z) < k) {
        t0 = t;
        found = true;
      }
      float r2 = 0.0;
      vec4 A = vec4(1,1,1,0);
      if (dot(A*p1,p1) < r2) {
        t0 = t;
        found = true;
      }
#else
      t0 = t;
      found = true;
#endif
    }
  }
  if (!found) return false;
  {
    vec4 p1 = p0+t0*r0;
    float x0 = eval(p1);
    vec2 eps = vec2(1e-2,0);
    normal = vec3 (eval(p1+eps.xyyy)-x0,
                   eval(p1+eps.yxyy)-x0,
                   eval(p1+eps.yyxy)-x0);
    normal = normalize(normal);
  }
  return true;
}
    
// Find the distance to the line in R3
void tryline(vec4 q, vec4 p, vec4 r, int color) {
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

vec4 solve(vec3 p, vec3 r, float min) {
  float t = 1e10;
  vec3 normal = vec3(1,0,0);
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
  //vec3 a = normalize(vec3(0,1,0));
  m = qmat(vec4(sin(theta)*a,cos(theta)));
  minv = qmat(vec4(sin(-theta)*a,cos(-theta)));
#endif
  //m = m0*m; minv = minv*minv0;

  //initlines();

  bool dosurface = true;
  bool dolines = false;
  dolines = true;
  if (dosurface && surface(vec4(p,1),vec4(r,0),min,t,normal)) {
    colorindex = 0;
    vec4 q = vec4(p+t*r,1);
    for (int i = 0; i < 27; i++) {
      if (i == nlines) break;
      vec4 p = minv*uLines[2*i];
      vec4 r = minv*uLines[2*i+1];
      if (abs(p.w) < abs(r.w)) {
        tryline(q,r,p,uColors[i]);
      } else {
        tryline(q,p,r,uColors[i]);
      }
    }
  }
  if (highlight) colorindex = 1;
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
  clock0 = uClock[0];
  clock1 = uClock[1];
  K = vec4(1,1,1,K0);

  light = normalize(vec3(0.0,1.0,-1.0));
  ambient = 0.6;
  diffuse = 1.0-ambient;

  vec2 uv = 2.0*fragCoord.xy/iResolution.xy - 1.0;
  vec3 p = vec3(0.0, 0.0, 6.0);
  vec3 r = normalize(vec3(iResolution.x/iResolution.y * uv.x, uv.y, -3.0));
  float k = -dot(p,r);
  //k = 0.0;
  p += k*r;
  float phi = clock0*0.123; 
  mat3 rmat = mat3(cos(phi), 0.0, -sin(phi),
                   0.0, 1.0, 0.0,
                   sin(phi), 0.0, cos(phi));
  light *= rmat;
  p *= rmat;
  r *= rmat;
  fragColor = solve(p,r,-k);
}

////////////////////////////////////////////////////////////////////////////////

void main() {
  mainImage(gl_FragColor, gl_FragCoord.xy);
}
