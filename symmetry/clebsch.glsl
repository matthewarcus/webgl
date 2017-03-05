precision highp float;
uniform float iGlobalTime;
uniform vec2 iResolution;

////////////////////////////////////////////////////////////////////////////////

const float PI =  3.141592654;
const float PHI = 1.618033989;
const vec4 I = vec4(1,1,1,1);

bool highlight = false; // Debugging

vec3 project(vec4 p) {
  return p.xyz/p.w;
}

bool isnan(float x) {
  return x != x;
}

// Quaternion multiplication as a matrix.
// w coordinate is real element of quaternion
mat4 qmat(vec4 q) {
  float x = q.x, y = q.y, z = q.z, t = q.w;
  return mat4( t,-z, y, x, 
               z, t,-x, y,
              -y, x, t, z,
              -x,-y,-z, t );
}

// Use vec2 return?
bool quadratic0(float A, float B, float C, out float x0, out float x1) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return false;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x0 = (-B-D)/(2.0*A);
   //x1 = (-B+D)/(2.0*A);
   x1 = C/(A*x0);
   if (x1 < x0) {
     float tmp = x0; x0 = x1; x1 = tmp;
   }
   return true;
}

bool quadratic(float A, float B, float C, float min, inout float x) {
  float x0,x1;
  if (!quadratic0(A,B,C,x0,x1)) return false;
  bool result = true;
  if (min <= x0 && x0 < x) x = x0;
  else if (min <= x1 && x1 < x) x = x0;
  else result = false;
  return result;
}

float evalcubic(float x, float a, float b, float c, float d) {
  return d+x*(c+x*(b+x*a));
}

bool cubic(float a, float b, float c, float d, float min, inout float x) {
  // Need a good way of handling ill-conditioned coeffs.
  if (a == 0.0) return quadratic(b,c,d,min,x);
  //if (abs(c) < 1e-3) highlight = true;
  if (abs(a) < 5e-3) {
    //highlight = true;
#if 1
    // This helps, but isn't perfect.
    float A = 1.0;
    float C = d/b;
    float B = (c-a*C)/b;
    return quadratic(A,B,C,min,x);
#endif
  }
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
  // return smallest +ve solution
  float Q = (a*a-3.0*b)/9.0;
  float R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  float R2 = R*R, Q3 = Q*Q*Q;
  bool found = false;
  if (R2 < Q3) {
    //if (nroots != 3) highlight = true;
    float X = clamp(R/sqrt(Q3),-1.0,1.0);
    float theta = acos(X);
    float S = sqrt(Q); // Q must be positive since 0 <= R2 < Q3
    float t = -2.0*S*cos(theta/3.0) - a/3.0;
    if (t >= min && (t < x || !found)) { x = t; found = true; }
    t = -2.0*S*cos((theta+2.0*PI)/3.0) - a/3.0;
    if (t >= min && (t < x || !found)) { x = t; found = true; }
    t = -2.0*S*cos((theta+4.0*PI)/3.0) - a/3.0;
    if (t >= min && (t < x || !found)) { x = t; found = true; }
  } else {
    //if (nroots != 1) highlight = true;
    float A = -sign(R)*pow(abs(R)+sqrt(R2-Q3),0.3333);
    float B = A == 0.0 ? 0.0 : Q/A;
    float t = (A+B) - a/3.0;
    if (t >= min && t < x) { x = t; found = true; }
  }
  return found;
}

const int NLINES = 27;
vec4 lines[NLINES*2];
int colors[NLINES];

void init() {
  lines[0] = vec4(1,-1,0,0); lines[1] = vec4(0,0,1,0);
  lines[2] = vec4(1,-1,0,0); lines[3] = vec4(0,0,0,1);
  lines[4] = vec4(1,0,-1,0); lines[5] = vec4(0,1,0,0);
  lines[6] = vec4(1,0,-1,0); lines[7] = vec4(0,0,0,1);
  lines[8] = vec4(1,0,0,-1); lines[9] = vec4(0,1,0,0);
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
  lines[34] = vec4(1,-1,PHI,0); lines[35] = vec4(PHI,0,1,-1);
  lines[36] = vec4(1,0,PHI,-1); lines[37] = vec4(PHI,-1,1,0);
  lines[38] = vec4(1,-1,0,PHI); lines[39] = vec4(PHI,0,-1,1);
  lines[40] = vec4(1,0,-1,PHI); lines[41] = vec4(PHI,-1,0,1);

  lines[42] = vec4(-1,1,PHI,0); lines[43] = vec4(0,PHI,1,-1);
  lines[44] = vec4(0,1,PHI,-1); lines[45] = vec4(-1,PHI,1,0);
  lines[46] = vec4(-1,1,0,PHI); lines[47] = vec4(0,PHI,-1,1);
  lines[48] = vec4(0,1,-1,PHI); lines[49] = vec4(-1,PHI,0,1);
  lines[50] = vec4(-1,0,1,PHI); lines[51] = vec4(0,-1,PHI,1);
  lines[52] = vec4(0,-1,1,PHI); lines[53] = vec4(-1,0,PHI,1);

  colors[0] = colors[2] = colors[4] = colors[6] = colors[8] = colors[10] = 0;
  colors[1] = colors[3] = colors[5] = colors[7] = colors[9] = colors[11] = 1;
  colors[12] = colors[13] = colors[14] = 2;
  colors[15] = colors[17] = colors[19] = colors[21] = colors[23] = colors[25] = 3;
  colors[16] = colors[18] = colors[20] = colors[22] = colors[24] = colors[26] = 4;
}

vec4 getColor(int i) {
  if (i == 0) return vec4(0.5,0.5,1,1);
  if (i == 1) return vec4(1,0,0,1);
  if (i == 2) return vec4(1,1,0,1);
  if (i == 3) return vec4(0,1,0,1);
  if (i == 4) return vec4(0,0,1,1);
  if (i == 5) return vec4(0,1,1,1);
  if (i == 6) return vec4(1,0,1,1);
  return vec4(0.5,0.5,0.5,1);
}

vec4 gradient(vec4 p) {
  float x = p.x; float y = p.y; float z = p.z; float w = p.w;
  float A = dot(p,I);
  return -3.0*vec4((A+x)*(A-x),(A+y)*(A-y),(A+z)*(A-z),(A+w)*(A-w));
}

mat4 m;
mat4 minv;

// 2x^2 - (x+y+z+w)^3
// x^3 + 3x^2(y+z+w) + 3x(yy+yz+yw+zz+zw+ww)
// Equation is: pp.p - (p.I)^3 = 0
bool surface(vec4 p, vec4 r, float min, inout float t, out vec3 normal) {
  p = m*p;
  r = m*r;
  float pI = dot(p,I);
  float rI = dot(r,I);
#if 0
  float A = dot(p*p,p)-pI*pI*pI;
  float B = 3.0*(dot(p*p,r) - pI*pI*rI);
  float C = 3.0*(dot(r*r,p) - pI*rI*rI);
  float D = dot(r*r,r) - rI*rI*rI;
#else
  // Possibly this is a better way.
  float A = dot(p*p - pI*pI*I, p);
  float B = 3.0*dot(p*p - pI*pI*I, r);
  float C = 3.0*dot(r*r - rI*rI*I, p);
  float D = dot(r*r - rI*rI*I, r);
#endif
  if (!cubic(D,C,B,A,min,t)) return false;
  vec4 tmp = p+t*r; // Actual surface point
  vec4 n4 = minv*gradient(tmp);  // Rotate n4 here
  normal = normalize(n4.xyz); // E3 component of P4 normal
  return true;
}
    
// A = pp.p - (p.I)^3
// B = 3(pp.r -(p.I)^2(r.I))
// C = 3(rr.p -(p.I)(r.I)^2)
// D = rr.r - (r.I)^3

const vec3 light = normalize(vec3(0.0,1.0,-1.0));
const float ambient = 0.6;
const float diffuse = 1.0-ambient;

vec4 solve(vec3 p, vec3 r, float min) {
  float t = 1e10;
  float u = 0.015;
  float u2 = u*u;
  vec3 q0, s0;
  int colorindex = -1;
  vec3 normal = vec3(1,0,0);
  float theta = iGlobalTime*0.1; 
  mat4 minv0 = 0.25*mat4(-1,1,1,1, 1,-1,1,1, 1,1,-1,1, -1,-1,-1,1);
  mat4 m0 = 0.25*mat4(-1,1,1,-1, 1,-1,1,-1, 1,1,-1,-1, 1,1,1,1);
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
  m = qmat(vec4(0.0,0.0,sin(theta),cos(theta)));
  minv = qmat(vec4(0.0,0.0,sin(theta),-cos(theta)));
#endif
  //m = m0*m; minv = minv*minv0;

  bool dosurface = true;
  bool dolines = true;
  if (dosurface && surface(vec4(p,1),vec4(r,0),min,t,normal)) {
    colorindex = 0;
  }
  if (dolines) {
    for (int i = 0; i < NLINES; i++) {
      vec3 q,s;
      {
        // Convert 2 projective points into an E3
        // point+direction line. I'm sure there is
        // a cleaner way to do this, but this seems
        // to work.
        vec4 p0 = minv*lines[2*i];
        vec4 p1 = minv*lines[2*i+1];
        float eps = 1e-1;
        if (abs(p0.w) < eps) p0 += p1;
        if (abs(p1.w) < eps) p1 += p0;
        if (abs(p0.w) < eps*eps) continue;
        if (abs(p1.w) < eps*eps) continue;
        q = project(p0);
        s = project(p1);
        s = normalize(s-q); // Line direction
      }
      //q = vec3(0,0,0);
      //s = vec3(0,1,0);
      vec3 p1 = p-q; // Move p to line space
      float rs = dot(r,s);
      float ps = dot(p1,s);
      float pr = dot(p1,r);
      float pp = dot(p1,p1);
      float A = 1.0 - rs*rs;
      float B = 2.0*(pr - ps*rs);
      float C = pp - ps*ps - u2;
      if (quadratic(A,B,C,min,t)) {
        q0 = q; s0 = s; colorindex = 1+colors[i];
        p1 += t*r; // Final point in line space
        p1 -= dot(p1,s)*s; // Subtract parallel part
        normal = normalize(p1);
      }
    }
  }
  if (highlight) colorindex = 1;
  if (colorindex == -1) return vec4(0,0,0,1);
  vec4 baseColor = getColor(colorindex);
  vec3 color;
  color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,normal));
  float specular = pow(max(0.0,dot(reflect(light,normal),r)),4.0);
  color += 0.7*specular*vec3(1.0,1.0,1.0);
  return vec4(sqrt(color),1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  init();
  vec2 uv = 2.0*fragCoord.xy/iResolution.xy - 1.0;
  vec3 p = vec3(0.0, 0.0, 10.0);
  vec3 r = normalize(vec3(iResolution.x/iResolution.y * uv.x, uv.y, -3.0));
  float k = -dot(p,r);
  //k = 0.0;
  p += k*r;
  fragColor = solve(p,r,-k);
}

////////////////////////////////////////////////////////////////////////////////

void main() {
  mainImage(gl_FragColor, gl_FragCoord.xy);
}
