precision highp float;

#define M_PI 3.14159
#define M_1_PI (1.0/M_PI)
uniform ivec2 uWindow;
uniform float uAlpha;
uniform int uSeed;
uniform vec4 params1;
uniform vec4 params2;

int seed = 12345678;
float frand() {
  const float two31 = 2147483648.0;
  seed *= 1664525;
  seed += 1013904223;
  return abs(float(seed))/two31;
}

int hash(int a) {
    a = (a + 61) + (a/(256*256));
    a = a + (a*8);
    a = a + (a/16);
    a = a * 0x27d4eb2d;
    a = a + (a / (128*256));
    return a;
}

struct Ray {
  vec3 o, d;
};

struct Sphere {
  float rad;       // radius
  vec3 p, e, c;    // position, emission, color
  int refl;        // reflection type (DIFFuse, SPECular, REFRactive)
  float k;         //r^2 - p^2
};

// material types, used in radiance()
const int DIFF = 0;
const int SPEC = 1;
const int REFR = 2;

// Avoid rounding error, both in the calculation of the QE coeffs
// and in solving the QE itself.
float intersect(Sphere s, Ray r) { // returns distance, 0 if nohit
  // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
  float B = dot(s.p-r.o,r.d);
  //float det=b*b-dot(op,op)+s.rad*s.rad;
  // r2 - dot(op,op) = (r2-pp) + 2po - oo
  // Use precalculated k = r^2 - p^2
  float C = s.k + dot(r.o, 2.0*s.p-r.o);
  float disc = B*B+C;
  if (disc < 0.0) return 0.0;
  float t1,t2;
  if (B < 0.0) {
    t1 = B - sqrt(disc); t2 = -C/t1;
  } else {
    t2 = B + sqrt(disc); t1 = -C/t2;
  }
  const float eps = 1e-4;
  if (t1 > eps) return t1;
  if (t2 > eps) return t2;
  return 0.0;
}

#if defined TREES
vec3 tc = vec3(0.0588, 0.361, 0.0941);
vec3 sc = vec3(1,1,1)*.7;
const int nSpheres = 13;
const int nLights = 1;

Sphere spheres[nSpheres];
void setupSpheres() {
  // center 50 40.8 62
  // floor 0
  // back  0
  //  Sphere(1e5, vec3(50, 1e5+100, 0),  vec3(1,1,1)*1,vec3(),DIFF), //lite
  //  Sphere(1e5, vec3(50, -1e5, 0),  vec3(),vec3(.3,.3,.1),DIFF), //grnd
  //  Sphere(1e5, vec3(50, 1e5+100, 0),  vec3(0.761, 0.875, 1.00)*1.3,vec3(),DIFF),
  //  //lite
  //Scene: radius, position, emission, color, material
  spheres[0] = Sphere(1e5, vec3(50.0, 1e5+130.0, 0.0),  vec3(1.0,1.0,1.0)*1.3,vec3(0),DIFF,0.0); //lite
  spheres[1] = Sphere(1e2, vec3(50.0, -1e2+2.0, 47.0),  vec3(0),vec3(1.0,1.0,1.0)*.7,DIFF,0.0); //grnd

  spheres[2] = Sphere(1e4, vec3(50.0, -30.0, 300.0)+vec3(-sin(50.0*M_PI/180.0),0.0,cos(50.0*M_PI/180.0))*1e4, vec3(0), vec3(1.0,1.0,1.0)*.99,SPEC,0.0);// mirr L
  spheres[3] = Sphere(1e4, vec3(50.0, -30.0, 300.0)+vec3(sin(50.0*M_PI/180.0),0.0,cos(50.0*M_PI/180.0))*1e4,  vec3(0), vec3(1.0,1.0,1.0)*.99,SPEC,0.0);// mirr R
  spheres[4] = Sphere(1e4, vec3(50.0, -30., -50.0)+vec3(-sin(30.*M_PI/180.),0,-cos(30.*M_PI/180.))*1e4,vec3(0), vec3(1.,1.,1.)*.99,SPEC,0.0);// mirr FL
  spheres[5] = Sphere(1e4, vec3(50.0, -30., -50.0)+vec3(sin(30.*M_PI/180.),0,-cos(30.*M_PI/180.))*1e4, vec3(0), vec3(1.,1.,1.)*.99,SPEC,0.0);// mirr


  spheres[6] = Sphere(4., vec3(50.0,6.*.6,47.),   vec3(0),vec3(.13,.066,.033), DIFF,0.0);//"tree"
  spheres[7] = Sphere(16.,vec3(50.0,6.*2.+16.*.6,47.),   vec3(0), tc,  DIFF,0.0);//"tree"
  spheres[8] = Sphere(11.,vec3(50.0,6.*2.+16.*.6*2.+11.*.6,47.),   vec3(0), tc,  DIFF,0.0);//"tree"
  spheres[9] = Sphere(7., vec3(50.0,6.*2.+16.*.6*2.+11.*.6*2.+7.*.6,47.),   vec3(0), tc,  DIFF, 0.0);//"tree"

  spheres[10] = Sphere(15.5,vec3(50.0,1.8+6.*2.+16.*.6,47.),   vec3(0), sc,  DIFF,0.0);//"tree"
  spheres[11] = Sphere(10.5,vec3(50.0,1.8+6.*2.+16.*.6*2.+11.*.6,47.),   vec3(0), sc,  DIFF,0.0);//"tree"
  spheres[12] = Sphere(6.5, vec3(50.0,1.8+6.*2.+16.*.6*2.+11.*.6*2.+7.*.6,47.),   vec3(0), sc,  DIFF,0.0);//"tree"
  for (int i = 0; i < nSpheres; i++) {
    float c = sqrt(dot(spheres[i].p,spheres[i].p));
    spheres[i].k = (spheres[i].rad-c)*(spheres[i].rad+c); // rad^2 - pp
  }
}
#else
const int nSpheres = 9;
const int nLights = 1;
Sphere spheres[nSpheres];
void setupSpheres() {
  //Scene: radius, position, emission, color, material
  spheres[0] = Sphere(1.5, vec3(50.,81.6-16.5,81.6),vec3(4.,4.,4.)*100.,  vec3(0), DIFF, 0.0);//Lite
  spheres[1] = Sphere(1e5, vec3(-1e5+99.,40.8,81.6),vec3(0),vec3(.25,.25,.75),DIFF, 0.0);//Rght
  spheres[2] = Sphere(1e5, vec3(50.,40.8, 1e5),     vec3(0),vec3(.75,.75,.75),DIFF, 0.0);//Back
  spheres[3] = Sphere(1e5, vec3(50.,40.8,-1e5+170.),vec3(0),vec3(0),           DIFF, 0.0);//Frnt
  spheres[4] = Sphere(1e5, vec3(50., 1e5, 81.6),    vec3(0),vec3(.75,.75,.75),DIFF, 0.0);//Botm
  spheres[5] = Sphere(1e5, vec3(50.,-1e5+81.6,81.6),vec3(0),vec3(.75,.75,.75),DIFF, 0.0);//Top
  spheres[6] = Sphere(1e5, vec3( 1e5+1.,40.8,81.6), vec3(0),vec3(.75,.25,.25),DIFF, 0.0);//Left
  spheres[7] = Sphere(16.5,vec3(27.,16.5,47.),      vec3(0),vec3(1.,1.,1.)*.999, SPEC, 0.0);//Mirr
  spheres[8] = Sphere(16.5,vec3(73.,16.5,78.),      vec3(0),vec3(1.,1.,1.)*.999, REFR, 0.0);//Glas
  for (int i = 0; i < nSpheres; i++) {
    float c = sqrt(dot(spheres[i].p,spheres[i].p));
    spheres[i].k = (spheres[i].rad-c)*(spheres[i].rad+c); // rad^2 - pp
  }
}
#endif

// gamma correction
//float gamma(float x){ return pow(clamp(x,0.0,1.0),1.0/2.2); }
float gamma(float x){ return pow(clamp(x,0.0,1.0),1.0/2.2); }

const float inf = 1e20;
bool intersect(Ray r, out float t, out int id, out Sphere s){
  t = inf;
  for (int i = 0; i < nSpheres; i++) {
    float d = intersect(spheres[i],r);
    if (d != 0.0 && d<t){
      t = d; id = i; s = spheres[i];
    }
  }
  return t<inf;
}

vec3 radiance(Ray r) {
  vec3 acc = vec3(0), fact = vec3(1);
  float E = 1.0;
  for (int depth = 0; depth < 10; depth++) {
    float t;                               // distance to intersection
    int id = 0;                            // id of intersected object
    Sphere obj;
    if (!intersect(r, t, id, obj)) return acc; // if miss, return black
    vec3 x = r.o+r.d*t; // The intersection point
    vec3 n = normalize(x-obj.p);
    vec3 nl = dot(n,r.d)<0.0 ? n : -n;
    bool into = dot(n,nl)>0.0;                // Ray from outside going in?
    vec3 f = obj.c;
    float p = max(f.x, max(f.y,f.z)); // max refl
    if (depth>5 || p != 0.0) {
      if (frand() < p) {
        f *= 1.0/p;
      } else {
        return acc+fact*(obj.e*E);
      }
    }
    if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
      float r1 = 2.*M_PI*frand();
      float r2 = frand();
      float r2s = sqrt(r2);
      vec3 w = nl;
      vec3 u = normalize(cross(abs(w.x)>.1?vec3(0.,1.,0.):vec3(1.,0.0,0.),w));
      vec3 v = cross(w,u);
      vec3 d = normalize(u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1.-r2));

      // Loop over any lights
      vec3 e = vec3(0);
      for (int i = 0; i < nLights; i++){
        Sphere s = spheres[i];
        vec3 sw = s.p-x;
        vec3 su = normalize(cross(abs(sw.x)>.1?vec3(0.,1.,0.):vec3(1.,0.,0.),sw));
        vec3 sv = cross(sw,su);
        float cos_a_max = sqrt(1.-s.rad*s.rad/dot(x-s.p,x-s.p));
        float eps1 = frand(), eps2 = frand();
        float cos_a = 1.-eps1+eps1*cos_a_max;
        float sin_a = sqrt(1.-cos_a*cos_a);
        float phi = 2.*M_PI*eps2;
        vec3 l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
        l = normalize(l);
        Sphere tmp;
        if (intersect(Ray(x,l), t, id, tmp) && id==i) {  // shadow ray
          float omega = 2.*M_PI*(1.-cos_a_max);
          e += f*(s.e*dot(l,nl)*omega)*M_1_PI;  // 1/pi for brdf
        }
      }
      //return radiance(Ray(x,d),depth,Xi,acc+fact.mult(obj.e*E + e), fact.mult(f),0);
      r = Ray(x,d);
      acc += fact*(obj.e*E + e);
      fact *= f;
      E = 0.0;
    } else if (obj.refl == SPEC) {             // Ideal SPECULAR reflection
      r = Ray(x,r.d-n*2.*dot(n,r.d));
      acc += fact*(obj.e*E);
      fact *= f;
      E = 1.0;
    } else {
      Ray reflRay = Ray(x, r.d-n*2.*dot(n,r.d));     // Ideal dielectric REFRACTION
      float nc = 1.0;
      float nt = 1.5;
      float nnt = into ? nc/nt : nt/nc;
      float ddn = dot(r.d,nl);
      float cos2t = 1.0 - nnt*nnt*(1.0 - ddn*ddn);
      if (cos2t < 0.0) {   // Total internal reflection
        r = reflRay;
        acc += fact*(obj.e);
        fact *= f;
        E = 1.0;
      } else {
        vec3 tdir = normalize(r.d*nnt - n*((into?1.:-1.)*(ddn*nnt+sqrt(cos2t))));
        float a = nt-nc, b = nt+nc;
        float R0 = a*a/(b*b);
        float c = 1.0 - (into ? -ddn : dot(tdir,n));
        float Re = R0+(1.0-R0)*c*c*c*c*c;
        float Tr = 1.0-Re, P = 0.25 + 0.5*Re, RP = Re/P, TP = Tr/(1.-P);
        if (frand() < P) {
          r = reflRay;
          acc += fact*(obj.e);
          fact *= f*RP;
          E = 1.0;
        } else {
          r = Ray(x,tdir);
          acc += fact*(obj.e);
          fact *= f*TP;
          E = 1.0;
        }
      }
    }
  }
  return acc;
}

const int nsamples = 2;

float width, height;
void render(float x, float y) {
  float time = params2[3];
  seed = hash(seed+uSeed);
  seed = hash(seed+int(time*1000.0));
  seed = hash(seed+int(x));
  seed = hash(seed+int(y));
  seed = hash(seed+int(x*y));
  Ray cam = Ray(vec3(50.,52.,295.6), normalize(vec3(0.,-0.042612,-1.))); // cam pos, dir
  vec3 cx = vec3(width*.5135/height,0.0,0.0);
  vec3 cy = normalize(cross(cx,cam.d))*.5135;
  vec3 col = vec3(0);
  for (int sy = 0; sy<2; sy++){     // 2x2 subpixel rows
    for (int sx = 0; sx<2; sx++){   // 2x2 subpixel cols
      vec3 r = vec3(0);
      for (int s = 0; s < nsamples; s++){
        float r1 = 2.0*frand(), dx=r1<1.0 ? sqrt(r1)-1.0 : 1.0-sqrt(2.0-r1);
        float r2 = 2.0*frand(), dy=r2<1.0 ? sqrt(r2)-1.0 : 1.0-sqrt(2.0-r2);
        vec3 d = cx*( ( (float(sx)+0.5 + dx)/2.0 + x)/width - .5) +
          cy*( ( (float(sy)+0.5 + dy)/2.0 + y)/height - 0.5) + cam.d;
        r += radiance(Ray(cam.o+d*140.0,normalize(d)))*(1.0/float(nsamples));
      } // Camera rays are pushed ^^^^^ forward to start in interior
      col += 0.25*vec3(gamma(r.x),gamma(r.y),gamma(r.z));
    }
  }
  gl_FragColor = vec4(clamp(col,0.0,1.0),uAlpha);
}

void main(){
  float time = params2[3];
  width = float(uWindow.x);
  height = float(uWindow.y);
  setupSpheres();
  float x = gl_FragCoord.x;
  float y = gl_FragCoord.y;
  render(x,y);
}
