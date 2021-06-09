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

const int CHAR_X = 88;
const int CHAR_Y = 89;
const int CHAR_Z = 90;

bool doMirror = true;
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

vec4 transform(vec4 p) {
  //return p;
  float k = 0.123*iTime;
  vec3 axis = vec3(0,0,-1);
  axis = normalize(axis);
  p = qmul(p,vec4(sin(k)*axis,cos(k)));
  p.xz = rotate(p.xz,0.55); // Pentagon centred
  //p.yz = rotate(p.yz,0.36);   // Triangle centred
  return p;
}

float Barth(vec4 p) {
  float x = p.x; float y = p.y;
  float z = p.z; float w = p.w;
  float x2 = x*x; float y2 = y*y;
  float z2 = z*z; float w2 = w*w;
  float A = 4.0*(phi2*x2-y2)*(phi2*y2-z2)*(phi2*z2-x2);
  float B = x2 + y2 + z2 - w2;
  float K = 0.5+0.5*sin(0.2*iTime); // K is mixed parameter, K=0.5 is "the" Barth Sextic
  K = max(K,0.01);
  return (1.0-K)*(1.0+2.0*phi)*w2*B*B - K*A;
}

float Fun(vec4 p) {
  p = transform(p);
  return Barth(p);
}

vec3 skybox(vec3 dir) {
  return pow(textureLod(iChannel0,dir,0.0).rgb,vec3(2.2));
}

vec3 getColor(vec4 q, vec3 eye, vec3 n) {
  q = transform(q);
#if 0
  return abs(n);
#elif 0
  return defaultColor;
#else
  vec3 color = skybox(reflect(eye,n));
  return color;
  return min(vec3(0.75),sqrt(color));
#endif
}

// Solution parameters.
const int iterations = 200;    // Maximum number of iterations
const float maxincrease = 1.03; // Largest allowed step increase.
const float maxstep = 1.0;     // The largest step that can be taken.
const float minstep = 0.001;  // The smallest step
const float initstep = 0.1;

vec3 solve(vec4 p0, const vec4 r) {
  float k0 = 0.0, k1;
  float a0 = Fun(p0), a1;
  bool bracketed = false;
  bool found = false;
  float step = initstep;
  vec4 p;
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
  float eps = 1e-3;
  vec2 delta = vec2(eps,0.0);
  assert(length(p) < 200.0);
#if 0
  p = p0 + k0*r; // Ensure p corresponds to k0 and a0
  vec3 n = vec3(Fun(p + delta.xyyy), Fun(p + delta.yxyy), Fun(p + delta.yyxy)) - a0;
#else
  // Not sure how much difference this makes
  vec3 n = vec3(Fun(p + delta.xyyy) - Fun(p - delta.xyyy),
                Fun(p + delta.yxyy) - Fun(p - delta.yxyy),
                Fun(p + delta.yyxy) - Fun(p - delta.yyxy));
#endif
  float grad = abs(length(n));
  n = normalize(n);

  vec3 eye = r.xyz;
  // Point normal towards eye
  if (dot(eye,n) > 0.0) n *= -1.0;
  vec3 baseColor = getColor(p,eye,n);
  color = baseColor;
  if (!doMirror) {
    color *= ambient;
    float k = dot(light,n);
    if (doDiffuse && k > 0.0) {
      color += baseColor*diffuse*k;
    }
    if (doSpecular && k > 0.0) {
      float specular = pow(max(0.0,dot(reflect(light,n),eye)),specpower);
      color += 0.8*specular*vec3(1);
    }
  }
  color = mix(skybox(r.xyz),color,max(0.0,(50.0-k0)/50.0));
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
  doMirror = key(CHAR_X);
  light = vec3(0.5,1,1);
  // Projection parameters
  float camera = 5.0;

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
      color += solve(vec4(p,1),vec4(r,0));
    }
  }
  color /= float(AA*AA);
  color = pow(color,vec3(0.4545));
  if (alert) color.r = 1.0;
  fragColor = vec4(color,1);
}
