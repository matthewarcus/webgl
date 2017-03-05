precision highp float;
uniform float iGlobalTime;
uniform vec2 iResolution;

////////////////////////////////////////////////////////////////////////////////

struct Ray {
  vec3 q;               // origin
  vec3 d;               // direction
};

struct Hit {
  float t;      // solution to p=q+t*d
  vec3 n;       // normal
  vec3 color;
};

struct Sphere {
  float r;      // radius
  vec3 p;       // centre
  vec3 color;
};

bool invert = true;

bool intersectSphere(Sphere s, Ray ray, out Hit hit) {
  vec3 p = s.p;
  float c = length(p);
  float r = s.r;
  // This inverts the sphere (in the origin).
  if (invert) {
    float k = 1.0/((c+r)*(c-r));
    r *= k; c *= k; p *= k;
  }
  vec3 q = ray.q, d = ray.d;
  float B = dot(q-p,d);
  float C = dot(q,q)-2.0*dot(q,p)+(c+r)*(c-r);
  float D = B*B - C;
  if (D < 0.0) return false;
  D = sqrt(D);
  float t,t1;
  if (B >= 0.0) {
    t = -B-D; t1 = C/t;
  } else {
    t1 = -B+D; t = C/t1;
  }
  if (t < 0.0) t = t1;
  if (t < 0.0) return false;
  hit = Hit(t, (q+t*d-p)/r, s.color);
  return true;
}

vec3 vertices[12];
int colors[12];
void setVertices() {
  float phi = 0.80902;
  // Three golden rectangle, oriented to
  // the three axes. 
  vertices[0] = vec3( 0.5, phi,0); //++0 A
  vertices[1] = vec3( 0.5,-phi,0); //+-0 B
  vertices[2] = vec3(-0.5, phi,0); //-+0 C
  vertices[3] = vec3(-0.5,-phi,0); //--0 D

  vertices[4] = vec3(0, 0.5, phi); //0++ B
  vertices[5] = vec3(0, 0.5,-phi); //0+- D
  vertices[6] = vec3(0,-0.5, phi); //0-+ C
  vertices[7] = vec3(0,-0.5,-phi); //0-- A

  vertices[8]  = vec3( phi,0, 0.5); //+0+ D
  vertices[9]  = vec3( phi,0,-0.5); //+0- C
  vertices[10] = vec3(-phi,0, 0.5); //-0+ A
  vertices[11] = vec3(-phi,0,-0.5); //-0- B

  colors[0] = colors[7] = colors[10] = 0;
  colors[1] = colors[4] = colors[11] = 1;
  colors[2] = colors[6] = colors[9] = 2;
  colors[3] = colors[5] = colors[8] = 3;
}

vec3 getColor(int i) {
  if (i == 0) return vec3(1,0,0);
  if (i == 1) return vec3(1,1,0);
  if (i == 2) return vec3(0,1,0);
  if (i == 3) return vec3(0,0,1);
  return vec3(1,1,1);
}

bool intersectScene(Ray r, out Hit hit) {
  float t = iGlobalTime*0.5;
  vec3 p = -(0.3 + normalize(vec3(1,1,1))); // Origin of icosahedron
  mat3 m = mat3(cos(t),-sin(t),0,sin(t),cos(t),0,0,0,1);
  setVertices();
  bool found = false;
  for (int i = 0; i < 12; i++) {
    Sphere s = Sphere(0.5, p+m*vertices[i], getColor(colors[i]));
    Hit hits;
    if (intersectSphere(s,r,hits) && (!found || hits.t < hit.t)) {
      hit = hits;
      found = true;
    }
  }
  return found;
}

const vec3 light = normalize(vec3(0.0,1.0,1.0));
const float ambient = 0.6;
const float diffuse = 1.0-ambient;

vec4 solve(Ray r) {
  Hit hit;
  if (!intersectScene(r,hit)) {
    return vec4(0,0,0,1);
  } else {
    vec3 n = hit.n;
    if (dot(r.d,n) > 0.0) n *= -1.0;
    vec3 baseColor = hit.color;
    vec3 color = baseColor.xyz*(ambient+(1.0-ambient)*dot(light,n));
    float specular = pow(max(0.0,dot(reflect(light,n),r.d)),4.0);
    color += 0.7*specular*vec3(1.0,1.0,1.0);
    return vec4(sqrt(color),1.0);
  }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 uv = 2.0*fragCoord.xy/iResolution.xy - 1.0;
  vec3 p = vec3(-0.5, -1.0, 6.0);
  vec3 d = normalize(vec3(iResolution.x/iResolution.y * uv.x, uv.y, -3));
  fragColor = solve(Ray(p,d));
}

////////////////////////////////////////////////////////////////////////////////

void main() {
  mainImage(gl_FragColor, gl_FragCoord.xy);
}
