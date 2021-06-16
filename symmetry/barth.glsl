#version 300 es
precision highp float;
uniform float iTime;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform sampler2D uSampler;
uniform samplerCube uCubeMap;
#define iChannel0 uCubeMap

uniform uvec4 uKeysPressed;
uniform uvec4 uKeysToggled;

out vec4 outColor;

void mainImage( out vec4 fragColor, in vec2 fragCoord );

void main() {
  mainImage(outColor, gl_FragCoord.xy);
}

#define texelFetch(x,y,z) (vec4(0))

#define key(code) ((uKeysToggled[code>>5]&(1U<<(code&31))) != 0U)

////////////////////////////////////////////////////////////////////////////////

bool alert = false;
void assert(bool b) {
  if (!b) alert = true;
}

#if !defined(key)
#define key(code) (texelFetch(iChannel3, ivec2((code),2),0).x != 0.0)
#endif

const int CHAR_A = 65;
const int CHAR_B = 66;
const int CHAR_C = 67;
const int CHAR_D = 68;
const int CHAR_E = 69;

const int CHAR_S = 83;
const int CHAR_T = 84;
const int CHAR_X = 88;
const int CHAR_Y = 89;
const int CHAR_Z = 90;

bool doMirror = false;
bool doGamma = true;
bool doSpecular = true;
bool doDiffuse = true;

const float ambient = 0.5;
const float diffuse = 0.5;
float specpower = 4.0;

const vec3 defaultColor = vec3(0.7,1.0,1.0);

const float PI =  3.141592654;
const float TWOPI = 2.0*PI;
const float phi = 1.618033;
const float phi2 = phi*phi;

vec3 light;

vec4 qmul(vec4 p, vec4 q) {
  vec3 P = p.xyz, Q = q.xyz;
  return vec4(p.w*Q+q.w*P+cross(P,Q),p.w*q.w-dot(P,Q));
}

vec2 rotate(vec2 p, float t) {
  return cos(t)*p + sin(t)*vec2(-p.y,p.x);
}

vec4 transform4(vec4 p) {
  if (!key(CHAR_S)) {
    p *= 2.0/dot(p,p);
    p.w -= 1.0;
  }
  float k = 0.1*PI*iTime;
  vec3 axis = vec3(0,0,1);
  axis = normalize(axis);
  p = qmul(p,vec4(sin(k)*axis,cos(k)));
  if (key(CHAR_C)) p.xz = rotate(p.xz,0.55); // Pentagon centred
  if (key(CHAR_T)) p.yz = rotate(p.yz,0.36);   // Triangle centred
  return p;
}

float Barth(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  float x2 = x*x; float y2 = y*y;
  float z2 = z*z; float w2 = w*w;
  float A = 4.0*(phi2*x2-y2)*(phi2*y2-z2)*(phi2*z2-x2);
  float B = x2 + y2 + z2 - w2;
  float K = 0.5+0.5*sin(-0.1*PI*iTime); // K is mixed parameter, K=0.5 is "the" Barth Sextic
  K = max(K,0.01);
  return (1.0-K)*(1.0+2.0*phi)*w2*B*B - K*A;
}

float Clebsch(vec4 p) {
  float t = dot(p,vec4(1));
  return dot(p*p,p)-t*t*t;
}

float Kummer(vec4 P) {
  float A = sqrt(2.0);
  float mu2 = 0.334 + 3.0*(1.0-sin(0.2*iTime));
  float x = P.x; float y = P.y;
  float z = P.z; float w = P.w;
  float p = w-z-A*x;
  float q = w-z+A*x;
  float r = w+z+A*y;
  float s = w+z-A*y;
  float lambda = (3.0*mu2-1.0)/(3.0-mu2);
  float k = x*x + y*y + z*z - mu2*w*w;
  return k*k-lambda*p*q*r*s;
}

float Sarti12(vec4 p){
  vec4 p2 = p*p;
  vec4 p4 = p2*p2;
  float l1 = dot(p2,p2);
  float l2 = p2.x*p2.y+p2.z*p2.w;
  float l3 = p2.x*p2.z+p2.y*p2.w;
  float l4 = p2.y*p2.z+p2.x*p2.w;
  float l5 = p.x*p.y*p.z*p.w;
  float s10 = l1*(l2*l3+l2*l4+l3*l4), s11 = l1*l1*(l2+l3+l4);
  float s12=l1*(l2*l2+l3*l3+l4*l4),    s51=l5*l5*(l2+l3+l4),  s234=l2*l2*l2+l3*l3*l3+l4*l4*l4;
  float s23p=l2*(l2+l3)*l3,   s23m=l2*(l2-l3)*l3; 
  float s34p=l3*(l3+l4)*l4,       s34m=l3*(l3-l4)*l4; 
  float s42p=l4*(l4+l2)*l2,       s42m=l4*(l4-l2)*l2;
  float Q12=dot(p,p); Q12=Q12*Q12*Q12; Q12=Q12*Q12; 
  float S12=33.*sqrt(5.)*(s23m+s34m+s42m)+19.*(s23p+s34p+s42p)+10.*s234-14.*s10+2.*s11-6.*s12-352.*s51+336.*l5*l5*l1+48.*l2*l3*l4;
  return 22.*Q12-243.*S12;
}

float Fun(vec3 p0) {
  vec4 p = vec4(p0,1);
  p = transform4(p);
  return Barth(p);
  //return Clebsch(p);
  //return Sarti12(p);
  //return Kummer(p);
}

vec3 skybox(vec3 dir) {
  return pow(textureLod(iChannel0,dir,0.0).rgb,vec3(2.2));
}

vec3 getColor(vec3 r, vec3 n) {
  if (!doMirror) {
    //return defaultColor;
    return abs(n);
  } else {
    vec3 color = skybox(reflect(r,n));
    return color;
    return min(vec3(0.75),sqrt(color));
  }
}

// Solution parameters.
const int iterations = 200;    // Maximum number of iterations
const float maxincrease = 1.03; // Largest allowed step increase.
const float maxstep = 1.0;     // The largest step that can be taken.
const float minstep = 0.001;  // The smallest step
const float initstep = 0.1;

vec3 solve(vec3 p0, vec3 r) {
  float k0 = 0.0, k1;
  vec3 p = p0;
  float a0 = Fun(p), a1;
  bool bracketed = false;
  bool found = false;
  float step = initstep;
  vec3 color = skybox(r.xyz);
  for (int i = 0; i < iterations; i++) {
    k1 = k0 + step;
    p = p0 + k1*r;
    a1 = Fun(p);
    if (a0*a1 <= 0.0) {
      // We can hit exactly 0 - this counts as bracketed.
      bracketed = true;
      break;
    } else {
      float step0 = step;
      step = a1*step/(a0-a1);
      step = abs(step);
      step = min(step,maxstep);
      // Don't grow step by more than a certain amount
      // Maybe detect overstepping & retreat.
      step = max(step,minstep);
      step = min(step,maxincrease*step0);
      k0 = k1; a0 = a1;
    }
  }
  if (!bracketed) return color;
  for (int i = 0; i < 10; i++) {
    // Once we are bracketed, just use bisection
    if ((k1-k0)/k0 < minstep) break;
    float k2 = (k0 + k1)/2.0;
    //x = x0+k2*a, y = y0+k2*b, z = z0+k2*c, w = 1.0;
    p = p0+k2*r;
    float a2 = Fun(p);
    if (a0*a2 <= 0.0) {
      k1 = k2; a1 = a2;
    } else {
      k0 = k2; a0 = a2;
    }
  }

  // Compute gradient & normal
  // Should probably scale eps here
  float eps = 1e-4;
  vec2 delta = k0*vec2(eps,0.0);
#if 1
  p = p0 + k0*r; // Ensure p corresponds to k0 and a0
  vec3 n = vec3(Fun(p + delta.xyy), Fun(p + delta.yxy), Fun(p + delta.yyx)) - a0;
#else
  // Not sure how much difference this makes
  vec3 n = vec3(Fun(p + delta.xyy) - Fun(p - delta.xyy),
                Fun(p + delta.yxy) - Fun(p - delta.yxy),
                Fun(p + delta.yyx) - Fun(p - delta.yyx));
#endif
  vec4 p4 = transform4(vec4(p,1));
  float grad = abs(length(n));
  n = normalize(n);

  // Point normal towards eye
  //assert(dot(r,n) < 0.0);
  if (dot(r,n) > 0.0) n *= -1.0;
  vec3 baseColor = getColor(r,n);
  color = baseColor;
  if (!doMirror) {
    color *= ambient;
    float k = dot(light,n);
    if (doDiffuse && k > 0.0) {
      color += baseColor*diffuse*k;
    }
    if (doSpecular && k > 0.0) {
      float specular = pow(max(0.0,dot(reflect(light,n),r)),specpower);
      color += 0.8*specular*vec3(1);
    }
  }
  //color = mix(skybox(r.xyz),color,max(0.0,(50.0-k0)/50.0));
  return color;
}

vec3 transform(vec3 p) {
  if (iMouse.x > 0.0) {
    float phi = (2.0*iMouse.x-iResolution.x)/iResolution.x*PI;
    float theta = (2.0*iMouse.y-iResolution.y)/iResolution.y*PI;
    p.yz = rotate(p.yz,theta);
    p.zx = rotate(p.zx,-phi);
  }
  //p.zx = rotate(p.zx,iTime * 0.2);
  return p;
}

void mainImage(out vec4 fragColor, vec2 fragCoord) {
  doMirror = !key(CHAR_X);
  light = vec3(0.5,1,1);
  // Projection parameters
  float camera = 6.0;
  vec3 p = vec3(0,0,camera);
  p = transform(p);
  light = transform(light); // Light moves with camera
  light = normalize(light);
  vec3 color = vec3(0);
  int AA = key(CHAR_A) ? 2 : 1;
  for (int i = 0; i < AA; i++) {
    for (int j = 0; j < AA; j++) {
      vec2 uv = (2.0*fragCoord+vec2(i,j) - iResolution.xy)/iResolution.y;
      vec3 r = normalize(vec3(uv,-2));
      r = transform(r);
      r = normalize(r);
      color += solve(p,r);
    }
  }
  color /= float(AA*AA);
  color = pow(color,vec3(0.4545));
  if (alert) color.r = 1.0;
  fragColor = vec4(color,1);
}
