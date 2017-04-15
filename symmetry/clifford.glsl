//#extension GL_OES_standard_derivatives : enable

precision highp float;
uniform sampler2D uSampler;
uniform samplerCube uCubeMap;
uniform vec2 iResolution;
uniform vec4 uClock;
uniform vec4 params1;
uniform vec4 params2;
uniform mat4 uMatrix;

float clock0;
float clock1;
float clock2;
float clock3;

////////////////////////////////////////////////////////////////////////////////

#define CUBEMAP

float C = 4.0;  // z coord of 3D camera
float K = 2.0;  // w coord of 4D camera
float K2 = 4.0; // K*K

const int NSTEPS = 16;
const float infinity = 1e10; //1.0/0.0;

const float GRID = 0.0;

bool highlight = false; // Debug
void assert(bool test) {
  if (!test) highlight = true;
}

const float PI =  3.141592654;
const float TWOPI =  2.0*PI;

void sort(inout float x1, inout float x2) {
  if (x1 > x2) {
    float t = x1; x1 = x2; x2 = t;
  }
}

vec4 getColor(int i) {
  if (i == 0) return vec4(0.5,0.5,1,1);
  if (i == 1) return vec4(0.5,1,0.5,1);
  if (i == 2) return vec4(1,0,0,1);
  if (i == 3) return vec4(1,1,0,1);
  if (i == 4) return vec4(0,1,0,1);
  if (i == 5) return vec4(0,0,1,1);
  if (i == 6) return vec4(0,1,1,1);
  if (i == 7) return vec4(1,0,1,1);
  return vec4(0.5,0.5,0.5,1);
}

bool quadratic(float A, float B, float C, out float x1, out float x2) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return false;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x1 = (-B-D)/(2.0*A);
   x2 = C/(A*x1);
   if (x1 > x2) {
     float x0 = x1; x1 = x2; x2 = x0;
   }
   return true;
}

// Quadratic but clamp discriminant to non-negative.
void cquadratic(float A, float B, float C, out float x1, out float x2) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) D = 0.0; // Clamp discriminant
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x1 = (-B-D)/(2.0*A);
   x2 = C/(A*x1);
   sort(x1,x2);
}

float F(vec4 p) {
  const vec4 A = vec4(1,1,0,0);
  return dot(A*p,p) - 0.5;
}

// We are projecting from w = K to w=0
// Assume the point is on the hypersphere
// so |unproject(p,n)| = 1 and
// project(unproject(p,n)) = p
vec4 unproject(vec3 p, int parity) {
  float A = dot(p,p)+K2;
  float B = -2.0*K2;
  float C = K2 - 1.0;
  float t1, t2;
  cquadratic(A,B,C,t1,t2);
  if (parity == 1) t1 = t2;
  return vec4(t1*p,K*(1.0-t1));
}

vec3 project(vec4 p) {
  float a = K/(K-p.w);
  return a*p.xyz;
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

mat4 m,minv;

float round(float x) {
  return floor(x+0.5);
}

bool gridpoint(vec3 p, int parity) {
  if (GRID == 0.0) return false;
  float N = GRID;
  const float halfwidth = 0.02;
  vec4 p4 = m*unproject(p,parity);
  float u = atan(p4.y,p4.x)/PI;
  if (abs(N*u - round(N*u)) < halfwidth) return true;
  float v = atan(p4.z,p4.w)/PI;
  if (abs(N*v - round(N*v)) < halfwidth) return true;
  return false;
}

// Taken from NR, simplified by using a fixed number of
// iterations and removing negative modulus case.
// Modulus is passed in as k^2 (_not_ 1-k^2 as in NR).
void sncndn(float u, float k2,
            out float sn, out float cn, out float dn) {
  float emc = 1.0-k2;
  float a,b,c;
  const int N = 4;
  float em[N],en[N];
  a = 1.0;
  dn = 1.0;
  for (int i = 0; i < N; i++) {
    em[i] = a;
    emc = sqrt(emc);
    en[i] = emc;
    c = 0.5*(a+emc);
    emc = a*emc;
    a = c;
  }
  // Nothing up to here depends on u, so
  // could be precalculated.
  u = c*u; sn = sin(u); cn = cos(u);
  if (sn != 0.0) {
    a = cn/sn; c = a*c;
    for(int i = N-1; i >= 0; i--) {
      b = em[i];
      a = c*a;
      c = dn*c;
      dn = (en[i]+a)/(b+a);
      a = c/b;
    }
    a = 1.0/sqrt(c*c + 1.0);
    if (sn < 0.0) sn = -a;
    else sn = a;
    cn = c*sn;
  }
}

// Complex sn. uv are coordinates in a rectangle, map to
// the upper half plane with a Jacobi elliptic function.
// Note: uses k^2 as parameter.
vec2 sn(vec2 z, float k2) {
  float snu,cnu,dnu,snv,cnv,dnv;
  sncndn(z.x,k2,snu,cnu,dnu);
  sncndn(z.y,1.0-k2,snv,cnv,dnv);
  float a = 1.0/(1.0-dnu*dnu*snv*snv);
  return a*vec2(snu*dnv, cnu*dnu*snv*cnv);
}

vec2 cn(vec2 z, float k2) {
  float snu,cnu,dnu,snv,cnv,dnv;
  sncndn(z.x,k2,snu,cnu,dnu);
  sncndn(z.y,1.0-k2,snv,cnv,dnv);
  float a = 1.0/(1.0-dnu*dnu*snv*snv);
  return a*vec2(cnu*cnv,-snu*dnu*snv*dnv);
}

vec2 getuv(float u, float v, out bool flip) {
  u += 0.05*clock2;
  v += 0.05*clock3;
  u = u - floor(u);
  v = v - floor(v);
#if defined CUBEMAP
  u = 2.0*u-1.0;
  v = 2.0*v-1.0;
  // -1 <= u,v <= 1
  if (u < -v) { u = -u; v = -v; }
  flip = true;
  if (u < 0.0) u = -u;
  else if (v < 0.0) v = -v;
  else if (u > 1.0) u = 2.0 - u;
  else if (v > 1.0) v = 2.0 - v;
  else flip = false;
  u = 1.0-u; v = 1.0-v;
  // Now rotate square by 5PI/4
  // and enlarge by 1.854;
  vec2 uv = cn(1.854*vec2(-u-v,u-v),0.5);
  //uv *= 0.97;
  uv = 0.5*(uv+1.0);
  return uv;
#else
  // 0 <= u,v < 1
  float u0 = v+u-1.0, v0 = v-u;
  u = u0; v = v0;
  if (u < -1.0 || u > 1.0) highlight = true;
  if (v < -1.0 || v > 1.0) highlight = true;
  // -1 < u,v < 1
  if (u < 0.0) { u = -u; v = -v; }
  if (v > 0.5) { v = 1.0-v; u = 1.0-u; }
  if (v < -0.5) { v = -1.0-v; u = 1.0-u; }
  v += 0.5;
  return vec2(u,v);
#endif
}

vec3 spherical(vec2 latlon) {
  float lat = latlon.x;
  float lon = latlon.y;
  return vec3(cos(lon)*cos(lat),sin(lat),sin(lon)*cos(lat));
}

float eval(vec3 p, vec3 r, float t, int parity) {
  // q is point on ray in ball
  vec3 q = p+t*r;
  vec4 q4 = unproject(q,parity);
  //if (abs(length(q4)-1.0) > 0.01) highlight = true;
  q4 = m*q4;;
  float x = F(q4);
  //if (x == 0.0) highlight = true;
  //if (x != x) highlight = true;
  return x;
}

float bisect(vec3 p, vec3 r, float t0, float t1, float x0, int parity) {
  //if (t0 >= t1) highlight = true;
  //if (eval(p,r,t0,parity) != x0) highlight = true;
  //if (eval(p,r,t0,parity) * eval(p,r,t1,parity) > 0.0) highlight = true;
  const int N = NSTEPS;
  for (int i = 0; i < N; i++) {
    if (x0 == 0.0) highlight = true;
    float t = 0.5*(t0+t1);
    float x = eval(p,r,t,parity);
    if (x*x0 <= 0.0) {
      t1 = t;
    } else {
      t0 = t; x0 = x;
    }
  }
  t0 = 0.5*(t0+t1);
  //if (abs(eval(p,r,t0,parity)) > 0.001) highlight = true;
  return t0;
}

bool solve(vec3 p, vec3 r, int parity, inout float tres) {  
  // p is camera, r is ray direction
  // (p + tr).(p + tr) = a2
  // p.p - a2 + 2tp.r + t^2 = 0
  float a2 = K2/(K2-1.0);
  float A = 1.0;
  float B = 2.0*dot(p,r);
  float C = dot(p,p) - a2;
  float t1,t2;
  if (!quadratic(A,B,C,t1,t2)) {
    return false;
  }
  //assert(t1 <= t2);
  //assert(t2 - t1 > 1e-2);
  // t1 and t2 are near and far R3 sphere intersections
  const int N = NSTEPS; // subdivide interval
  float t0 = t1;
  float x0 = eval(p,r,t0,parity);

  // If we have already hit the surface, return.
  // Allowing negative x0 avoids problems with the join
  if (abs(x0) < 1e-3 && t0 < tres) {
    tres = t0;
    //highlight = true;
    return true;
  }
  for (int i = 1; i <= N; i++) {
    float t = (float(N-i)*t1 + float(i)*t2)/float(N);
    float x = eval(p,r,t,parity);
    if (x0*x <= 0.0) {
      float t = bisect(p,r,t0,t,x0,parity);
      if (t < tres) {
        tres = t;
        return true;
      } else {
        return false;
      }
    } else {
      x0 = x; t0 = t;
    }
  }
  return false;
}

vec4 getMapColor(vec3 q, int color) {
  bool flip;
  vec4 p4 = m*unproject(q,color);
  float u = 0.5+atan(p4.y,p4.x)/TWOPI;
  float v = 0.5+atan(p4.w,p4.z)/TWOPI;
  vec2 uv = getuv(u,v,flip);
#if !defined CUBEMAP
  return texture2D(uSampler,uv);
#else
  // Convert uv to spherical coords
  // We should be able to do this more directly
  u = uv.x; v = uv.y;
  u = 2.0*u-1.0;
  v = 2.0*v-1.0;
  float r = sqrt(u*u + v*v);
  float lon = atan(u,v);
  float lat = 0.5*PI-2.0*atan(r);
  if (flip) lat = -lat;
  //float coslat = 2*r/(1+r*r);
  //float sinlat = (1-r*r)/(1+r*r);
  vec3 cubedir = spherical(vec2(lat,lon));
  return textureCube(uCubeMap,cubedir);
#endif
}

void testImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 uv = fragCoord.xy/iResolution.xy;
  bool flip;
  uv = getuv(uv.x,uv.y,flip);
#if !defined CUBEMAP
  fragColor = texture2D(uSampler,uv);
#else
  // Convert uv to spherical coords
  // We should be able to do this more directly
  float u = 2.0*uv.x-1.0, v = 2.0*uv.y-1.0;
  float r = sqrt(u*u + v*v);
  float lon = atan(u,v);
  float lat = 0.5*PI-2.0*atan(r);
  if (flip) lat = -lat;
  //float coslat = 2*r/(1+r*r);
  //float sinlat = (1-r*r)/(1+r*r);
  vec3 cubedir = spherical(vec2(lat,lon));
  fragColor = textureCube(uCubeMap,cubedir);
#endif
}

vec3 light;
float ambient;
float diffuse;

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  light = normalize(vec3(0.0,1.0,-1.0)); // -ve z is lit from front
  ambient = 0.5;
  diffuse = 0.5;

  float ar = iResolution.x/iResolution.y; // aspect ratio
  vec2 uv = 2.0*(fragCoord.xy/iResolution.xy - 0.5);
  // -1 <= u,v <= 1
  float scale = params1[1]; 
  vec2 screenpoint = scale*vec2(ar*uv.x, uv.y);
  vec3 p = vec3(0.0, 0.0, -C);                  // point of projection
  vec3 r = normalize(vec3(screenpoint,0) - p); // ray

  // Try to shift projection point into scene
  // There are better ways to do this
  if (C > 4.0) p += (C-4.0)*r;
  if (C < -4.0) p += (C+4.0)*r;

  // Set matrices
  float theta = 0.1*clock1;
  float phi = 0.1234*clock0;
  {
    // Quaternion rotation about dir
    vec3 dir = normalize(vec3(0,1,0));
    m = qmat(vec4(-sin(theta)*dir,cos(theta)));
    minv = qmat(vec4(sin(theta)*dir,cos(theta)));
  }
  // model rotation about y axis
  mat3 rmat = mat3(cos(phi),0,-sin(phi),
                   0,1,0,
                   sin(phi),0,cos(phi));
  float t = infinity;
  int color = -1;
  // Move the camera & ray into model coordinates
  p = mat3(uMatrix) * p;
  r = mat3(uMatrix) * r;
  p = rmat*p;
  r = rmat*r;
  if (solve(p,r,0,t)) color = 0;
  if (solve(p,r,1,t)) color = 1;
  vec3 q = p+t*r;
  if (color == -1) {
    fragColor = vec4(0,0,0,1);
  } else {
    vec3 normal;
    {
      vec4 q4 = m*unproject(q,color);
      // Reproject from q4 for consistency
      vec3 q0 = project(minv*q4);
      //assert(F(q4) < 1e-6);
      float eps = 1e-2;
      float x = q4.x, y = q4.y, z = q4.z, w = q4.w;
      // Find two other nearby points on the surface
      vec4 q4a = q4 + eps*vec4(y,-x,0,0);
      vec4 q4b = q4 + eps*vec4(0,0,w,-z);
      //assert(abs(F(q4a)) < 1e-4);
      //assert(abs(F(q4b)) < 1e-4);
      // Transform back to model space
      vec3 qa = project(minv*q4a);
      vec3 qb = project(minv*q4b);
      //assert(length(qa-q) > 1e-3);
      //assert(length(qb-q) > 1e-3);
      //assert(F(q4b) < 1e-2);
      // And find the normal to the triangle.
      normal = normalize(cross(qa-q0,qb-q0));
    }
    //normal = normalize(cross(dFdx(q), dFdy(q)));

    // with RH rule, n should be facing in right direction,
    if (dot(r,normal) > 0.0) {
      normal = -normal;
    }
    //if (gridpoint(q,color)) color = 2;
    vec4 baseColor;
    if (color > 1) {
      baseColor = getColor(color);
    } else {
      baseColor = getMapColor(q,color);
    }
    //light = mat3(uMatrix) * light;
    light = rmat*light;
    vec3 color = baseColor.xyz*(ambient+diffuse*max(0.0,dot(light,normal)));
    float specular = pow(max(0.0,dot(reflect(light,normal),r)),4.0);
    color += 0.5*specular*vec3(1.0,1.0,1.0);
    fragColor = vec4(color,1.0);
    //fragColor = vec4(sqrt(color),1.0);
  }
  if (highlight) {
    fragColor.r = 1.0;
  }
}

////////////////////////////////////////////////////////////////////////////////

void main() {
  clock0 = uClock[0];
  clock1 = uClock[1];
  clock2 = uClock[2];
  clock3 = uClock[3];
  C = 4.0 + 1.0*params1[2];
  K = 2.0 + 1.0*params1[3];
  K2 = K*K;
  mainImage(gl_FragColor, gl_FragCoord.xy);
  //testImage(gl_FragColor, gl_FragCoord.xy);
}
