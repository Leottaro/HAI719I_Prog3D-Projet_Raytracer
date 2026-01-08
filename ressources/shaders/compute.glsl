#version 430 core
const float PI = 3.14159265359;
const uint UINT_MAX = 0xFFFFFFFF;
const float FLT_MAX = 3.402823466e+38;

// =====================================================================================================
// ========================================= INPUTS STRUCTURES =========================================
// =====================================================================================================

struct PhongSettings {
  bool ENABLED;
  uint SHADOW_RAYS;
};

struct MaterialSettings {
  bool ENABLE_MIRROR;
  bool ENABLE_GLASS;
  float AIR_INDEX_MEDIUM;
};

struct MeshSettings {
  bool ENABLE_INTERPOLATION;
};

struct KdTreeSettings {
  float EPSILON;
  uint MAX_LEAF_SIZE;
};

struct BonusSettings {
  bool ENABLE_TEXTURES;
};

struct Settings {
  float EPSILON;
  uint NSAMPLES;
  uint MAX_BOUNCES;
  uint SCREEN_WIDTH;
  uint SCREEN_HEIGHT;

  PhongSettings Phong;
  MaterialSettings Material;
  MeshSettings Mesh;
  KdTreeSettings KdTree;
  BonusSettings Bonus;
};

const uint Material_DiffUSE_PHONG = 1;
const uint Material_Glass = 2;
const uint Material_Mirror = 3;
struct Material {
  int image_id;
  vec3 ambient_material;
  vec3 diffuse_material;
  vec3 specular_material;
  float shininess;

  float index_medium;
  float transparency;

  uint type;
};

struct Sphere {
  vec3 m_center;
  float m_radius;
  Material material;
};

struct Square {
  vec3 m_normal;
  vec3 m_bottom_left;
  vec3 m_right_vector;
  vec3 m_up_vector;
  Material material;
};

const uint LightType_Spherical = 1;
const uint LightType_Quad = 2;
struct Light {
  vec3 material;
  bool isInCamSpace;
  uint type;
  Sphere sphere;
  Square quad;
  float powerCorrection;
};

struct MeshVertex {
  vec3 position;
  vec3 normal;
  float u;
  float v;
};

struct MeshTriangle {
  uint v0;
  uint v1;
  uint v2;
};

struct Mesh {
  Material material;
  uint nb_vertices;
  uint nb_triangles;
  MeshVertex vertices[1000];
  MeshTriangle triangles[1000];
};

// ======================================================================================================
// =============================================== INPUTS ===============================================
// ======================================================================================================

layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(binding = 0, rgba32f) uniform image2D imgOutput;

uniform dvec2 nearAndFarPlanes;
uniform dmat4 modelviewInverse;
uniform dmat4 projectionInverse;
uniform float min_t;
uniform float max_t;

uniform Settings settings;
uniform uint nb_spheres;
uniform uint nb_squares;
uniform uint nb_lights;
uniform uint nb_meshes;

uniform Sphere spheres[10];
uniform Square squares[10];
uniform Light lights[10];
uniform Mesh meshes[10];

ivec3 screen_coords = ivec3(gl_GlobalInvocationID.xyz);
uint seed = uint(screen_coords.y * settings.SCREEN_WIDTH * settings.NSAMPLES + screen_coords.x * settings.NSAMPLES + screen_coords.z);

// ======================================================================================================
// ========================================= STRUCTURES HELPERS =========================================
// ======================================================================================================

const uint MAX_BOUNCES = 100;
struct Ray {
  vec3 pos;
  vec3 dir;
  uint nb_bounces;
  float index_mediums[MAX_BOUNCES];
  uint object_types[MAX_BOUNCES];
  uint object_indices[MAX_BOUNCES];
};
Ray newRay(vec3 position, vec3 direction) {
  Ray ray;
  ray.pos = position;
  ray.dir = normalize(direction);
  ray.nb_bounces = 1;
  ray.index_mediums[0] = settings.Material.AIR_INDEX_MEDIUM;
  ray.object_types[0] = UINT_MAX;
  ray.object_indices[0] = UINT_MAX;
  return ray;
}
Ray newRay(vec3 position, vec3 direction, Ray input_ray) {
  Ray ray;
  ray.pos = position;
  ray.dir = normalize(direction);
  ray.nb_bounces = input_ray.nb_bounces;
  ray.index_mediums = input_ray.index_mediums;
  ray.object_types = input_ray.object_types;
  ray.object_indices = input_ray.object_indices;
  return ray;
}

vec2 getRayInOutIndexMediums(inout Ray ray, in float index_medium, in uint object_type, in uint object_index) {
  if(ray.object_types[ray.nb_bounces - 1] == object_type && ray.object_indices[ray.nb_bounces - 1] == object_index) {
    ray.nb_bounces--;
    return vec2(ray.index_mediums[ray.nb_bounces], ray.index_mediums[ray.nb_bounces - 1]);
  } else {
    ray.index_mediums[ray.nb_bounces] = index_medium;
    ray.object_types[ray.nb_bounces] = object_type;
    ray.object_indices[ray.nb_bounces] = object_index;
    ray.nb_bounces++;
    return vec2(ray.index_mediums[ray.nb_bounces - 2], ray.index_mediums[ray.nb_bounces - 1]);
  }
  // return vec2(nL, nT);
}

struct Triangle {
  vec3 m_normal;
  vec3 m_c[3];
  float area;
};
void updateAreaAndNormal(inout Triangle triangle) {
  vec3 nNotNormalized = cross(triangle.m_c[1] - triangle.m_c[0], triangle.m_c[2] - triangle.m_c[0]);
  float norm = nNotNormalized.length();
  triangle.m_normal = nNotNormalized / norm;
  triangle.area = norm / 2.;
}
Triangle newTriangle(in vec3 c0, in vec3 c1, in vec3 c2) {
  Triangle triangle;
  triangle.m_c[0] = c0;
  triangle.m_c[1] = c1;
  triangle.m_c[2] = c2;
  updateAreaAndNormal(triangle);
  return triangle;
}
bool isParallelTo(in Ray ray, in Triangle triangle) {
  return abs(dot(ray.dir, triangle.m_normal)) <= settings.EPSILON;
}
vec3 getIntersectionPointWithSupportPlane(in Ray ray, in Triangle triangle, out float t) {
  if(isParallelTo(ray, triangle)) {
    return vec3(0.);
  }

  t = -(dot(triangle.m_normal, ray.pos - triangle.m_c[0])) / dot(triangle.m_normal, ray.dir);
  return ray.pos + t * ray.dir;
}
vec3 computeBarycentricCoordinates(in vec3 p, in Triangle triangle) {
  return vec3(newTriangle(p, triangle.m_c[1], triangle.m_c[2]).area, //
  newTriangle(triangle.m_c[0], p, triangle.m_c[2]).area,  // 
  newTriangle(triangle.m_c[0], triangle.m_c[1], p).area //
  ) / triangle.area - vec3(settings.EPSILON);
}

struct RaySphereIntersection {
  bool intersectionExists;
  float t;
  float theta, phi;
  vec3 intersection;
  vec3 secondintersection;
  vec3 normal;
};
RaySphereIntersection newRaySphereIntersection() {
  RaySphereIntersection intersection;
  intersection.intersectionExists = false;
  intersection.t = FLT_MAX;
  return intersection;
}

struct RaySquareIntersection {
  bool intersectionExists;
  float t;
  float u;
  float v;
  vec3 intersection;
  vec3 normal;
};
RaySquareIntersection newRaySquareIntersection() {
  RaySquareIntersection intersection;
  intersection.intersectionExists = false;
  intersection.t = FLT_MAX;
  return intersection;
}

struct RayTriangleIntersection {
  bool intersectionExists;
  float t;
  float w0;
  float w1;
  float w2;
  float u;
  float v;
  uint tIndex;
  vec3 intersection;
  vec3 normal;
};
RayTriangleIntersection newRayTriangleIntersection() {
  RayTriangleIntersection intersection;
  intersection.intersectionExists = false;
  intersection.t = FLT_MAX;
  return intersection;
}

struct RaySceneIntersection {
  bool intersectionExists;
  float t;
  float u;
  float v;
  vec3 intersection;
  vec3 normal;
  Material material;
  uint typeOfIntersectedObject;
  uint objectIndex;
  RaySphereIntersection raySphereIntersection;
  RaySquareIntersection raySquareIntersection;
  RayTriangleIntersection rayMeshIntersection;
};
RaySceneIntersection newRaySceneIntersection(in float _max_t) {
  RaySceneIntersection intersection;
  intersection.intersectionExists = false;
  intersection.t = _max_t;
  return intersection;
}

const uint TriangleIntersection = 0;
const uint SphereIntersection = 1;
const uint SquareIntersection = 2;
const uint LightIntersection = 3;

// =====================================================================================================
// ============================================= FUNCTIONS =============================================
// =====================================================================================================

// HELPERS

// https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
// A single iteration of Bob Jenkins' One-At-A-Time hashing algorithm.
uint hash(uint x) {
  x += (x << 10u);
  x ^= (x >> 6u);
  x += (x << 3u);
  x ^= (x >> 11u);
  x += (x << 15u);
  return x;
}

  // Compound versions of the hashing algorithm I whipped together.
uint hash(uvec2 v) {
  return hash(v.x ^ hash(v.y));
}
uint hash(uvec3 v) {
  return hash(v.x ^ hash(v.y) ^ hash(v.z));
}
uint hash(uvec4 v) {
  return hash(v.x ^ hash(v.y) ^ hash(v.z) ^ hash(v.w));
}

  // Construct a float with half-open range [0:1] using low 23 bits.
  // All zeroes yields 0.0, all ones yields the next smallest representable value below 1.0.
float floatConstruct(uint m) {
  const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
  const uint ieeeOne = 0x3F800000u; // 1.0 in IEEE binary32

  m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
  m |= ieeeOne;                          // Add fractional part to 1.0

  float f = uintBitsToFloat(m);       // Range [1:2]
  return f - 1.0;                        // Range [0:1]
}

  // Pseudo-random value in half-open range [0:1].
float random(float x) {
  return floatConstruct(hash(floatBitsToUint(x)));
}
float random(vec2 v) {
  return floatConstruct(hash(floatBitsToUint(v)));
}
float random(vec3 v) {
  return floatConstruct(hash(floatBitsToUint(v)));
}
float random(vec4 v) {
  return floatConstruct(hash(floatBitsToUint(v)));
}

vec3 sphericalCoordinatesToEuclidean(in float theta, in float phi) {
  float sinPhi = sin(phi);
  float x = sinPhi * sin(theta);
  float y = cos(phi);
  float z = sinPhi * cos(theta);

  return vec3(x, y, z);
}

vec3 euclideanCoordinatesToSpherical(in vec3 pos) {
  float R = length(pos);
  float theta = atan(pos.x, pos.z); // azimuth around y-axis, 0..2π
  if(theta < 0.) {
    theta += 2. * PI;
  }

  float phi = acos(pos.y / R); // polar angle from +y axis, 0..π

  return vec3(theta, phi, R);
}

vec3 getLightCentralPos(in Light light) {
  if(light.type == LightType_Quad)
    return light.quad.m_bottom_left + 0.5 * (light.quad.m_right_vector + light.quad.m_up_vector);
  if(light.type == LightType_Spherical)
    return light.sphere.m_center;
}

// INTERSECTIONS

void setRaySphereIntersection(inout RaySceneIntersection result, in RaySphereIntersection intersection, in Material material, in uint object_i) {
  result.intersectionExists = true;
  result.t = intersection.t;
  result.u = intersection.theta / (2. * PI);
  result.v = intersection.phi / PI;
  result.intersection = intersection.intersection;
  result.material = material;
  result.normal = intersection.normal;
  result.typeOfIntersectedObject = SphereIntersection;
  result.objectIndex = object_i;
  result.raySphereIntersection = intersection;
}
void setRaySquareIntersection(inout RaySceneIntersection result, in RaySquareIntersection intersection, in Material material, in uint object_i) {
  result.intersectionExists = true;
  result.t = intersection.t;
  result.u = intersection.u;
  result.v = intersection.v;
  result.intersection = intersection.intersection;
  result.normal = intersection.normal;
  result.material = material;
  result.typeOfIntersectedObject = SquareIntersection;
  result.objectIndex = object_i;
  result.raySquareIntersection = intersection;
}
void setRayTriangleIntersection(inout RaySceneIntersection result, in RayTriangleIntersection intersection, in Material material, in uint object_i) {
  result.intersectionExists = true;
  result.t = intersection.t;
  result.u = intersection.u;
  result.v = intersection.v;
  result.intersection = intersection.intersection;
  result.normal = intersection.normal;
  result.material = material;
  result.typeOfIntersectedObject = TriangleIntersection;
  result.objectIndex = object_i;
  result.rayMeshIntersection = intersection;
}

RaySquareIntersection intersectSquare(in Ray ray, in Square square, in bool can_intersect_behind) {
  RaySquareIntersection intersection = RaySquareIntersection(false, 0., 0., 0., vec3(0.), vec3(0.));

  const vec3 O = ray.pos;
  const vec3 D = ray.dir;
  const vec3 C = square.m_bottom_left;
  const vec3 N = square.m_normal;

  // nous sommes derrière le carré donc pas d'intersection
  if(!can_intersect_behind && square.material.type != Material_Glass && dot(D, N) >= 0.) {
    return intersection;
  }

  float t = -(dot(N, O - C)) / dot(N, D);

  vec3 P = O + t * D;
  vec3 PC = P - C;
  float u = dot(PC, square.m_right_vector) / dot(square.m_right_vector, square.m_right_vector);
  float v = dot(PC, square.m_up_vector) / dot(square.m_up_vector, square.m_up_vector);

  if(settings.EPSILON <= u && u <= 1. && settings.EPSILON <= v && v <= 1.) {
    intersection.intersectionExists = true;
    intersection.t = t;
    intersection.u = u;
    intersection.v = v;
    intersection.intersection = P;
    intersection.normal = N;
  }

  return intersection;
}

RaySphereIntersection intersectSphere(in Ray ray, in Sphere sphere) {
  RaySphereIntersection intersection = newRaySphereIntersection();

  vec3 O = ray.pos;
  vec3 D = ray.dir;
  vec3 C = sphere.m_center;
  float r = sphere.m_radius;
  vec3 CO = O - C;

  float a = dot(D, D);
  float b = 2. * dot(D, CO);
  float c = dot(CO, CO) - r * r;
  float delta = b * b - 4. * a * c;

  if(delta < 0.) {
    return intersection;
  }

  float sqrt_delta = sqrt(delta);
  float t = (-b - sqrt_delta) / (2. * a);
  float t2 = (-b + sqrt_delta) / (2. * a);
  bool is_outside = t >= settings.EPSILON;

  intersection.intersectionExists = true;
  intersection.t = is_outside ? t : t2;
  intersection.intersection = is_outside ? O + t * D : O + t2 * D;
  intersection.normal = normalize(is_outside ? intersection.intersection - C : C - intersection.intersection);
  intersection.secondintersection = is_outside ? O + t2 * D : O + t * D;

  vec3 spherical_pos = euclideanCoordinatesToSpherical(intersection.normal);
  intersection.theta = spherical_pos.x;
  intersection.phi = spherical_pos.y;

  return intersection;
}

RayTriangleIntersection intersectTriangle(in Ray ray, in Triangle triangle) {
  RayTriangleIntersection intersection = newRayTriangleIntersection();

  // 1) check that the ray is not parallel to the triangle
  if(isParallelTo(ray, triangle)) {
    return intersection;
  }

  // calculate the intersectionRayTriangleIntersection
  float t = 0;
  vec3 p = getIntersectionPointWithSupportPlane(ray, triangle, t);

  // 2) check that the triangle is "in front of" the ray
  if(dot(p - ray.pos, ray.dir) < 0) {
    return intersection;
  }

  // 3) check that the intersection point is inside the triangle:
  vec3 ws = computeBarycentricCoordinates(p, triangle);
  if(any(lessThan(ws, vec3(0.))) || any(greaterThan(ws, vec3(1.))) || abs(dot(ws, vec3(1.)) - 1.) > 0.) { // TODO: 0 < w0+w1+w2 < 1
    return intersection;
  }

  // 4) Finally, if all conditions were met, then there is an intersection!
  intersection.intersectionExists = true;
  intersection.t = t;
  intersection.w0 = ws.x;
  intersection.w1 = ws.y;
  intersection.w2 = ws.z;
  // intersection.u;
  // intersection.v;
  // intersection.tIndex;
  intersection.intersection = p;
  // intersection.normal;
  return intersection;
}

RayTriangleIntersection intersectMesh(in Ray ray, in Mesh mesh) {
  RayTriangleIntersection closestIntersection = newRayTriangleIntersection();

  for(uint i = 0; i < mesh.nb_triangles; i++) {
    MeshTriangle mesh_triangle = mesh.triangles[i];
    MeshVertex v0 = mesh.vertices[mesh_triangle.v0];
    MeshVertex v1 = mesh.vertices[mesh_triangle.v1];
    MeshVertex v2 = mesh.vertices[mesh_triangle.v2];
    Triangle triangle = newTriangle(v0.position, v1.position, v2.position);
    RayTriangleIntersection intersection = intersectTriangle(ray, triangle);
    if(intersection.intersectionExists && intersection.t < closestIntersection.t) {
      if(settings.Mesh.ENABLE_INTERPOLATION) {
        intersection.normal = normalize(v0.normal * intersection.w0 + v1.normal * intersection.w1 + v2.normal * intersection.w2);
        intersection.u = v0.u * intersection.w0 + v1.u * intersection.w1 + v2.u * intersection.w2;
        intersection.v = v0.v * intersection.w0 + v1.v * intersection.w1 + v2.v * intersection.w2;
      } else {
        intersection.normal = normalize(v0.normal);
        intersection.u = v0.u;
        intersection.v = v0.v;
      }
      intersection.tIndex = i;
      closestIntersection = intersection;
    }
  }
  return closestIntersection;
}

// SCENE.H

RaySceneIntersection computeIntersection(in Ray ray, in float min_t, in float max_t, in bool intersect_lights) {
  RaySceneIntersection result = newRaySceneIntersection(max_t);

  for(uint i = 0; i < nb_spheres; i++) {
    RaySphereIntersection intersection = intersectSphere(ray, spheres[i]);
    if(intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
      setRaySphereIntersection(result, intersection, spheres[i].material, i);
    }
  }

  for(uint i = 0; i < nb_squares; i++) {
    RaySquareIntersection intersection = intersectSquare(ray, squares[i], false);
    if(intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
      setRaySquareIntersection(result, intersection, squares[i].material, i);
    }
  }

  for(uint i = 0; i < nb_meshes; i++) {
    RayTriangleIntersection intersection = intersectMesh(ray, meshes[i]);
    if(intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
      setRayTriangleIntersection(result, intersection, meshes[i].material, i);
    }
  }

  if(intersect_lights) {
    for(uint i = 0; i < nb_lights; i++) {
      if(lights[i].type == LightType_Quad) {
        RaySquareIntersection intersection = intersectSquare(ray, lights[i].quad, true);
        if(intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
          setRaySquareIntersection(result, intersection, lights[i].quad.material, i);
          result.typeOfIntersectedObject = LightIntersection;
        }
      } else if(lights[i].type == LightType_Spherical) {
        RaySphereIntersection intersection = intersectSphere(ray, lights[i].sphere);
        if(intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
          setRaySphereIntersection(result, intersection, lights[i].sphere.material, i);
          result.typeOfIntersectedObject = LightIntersection;
        }
      } else {
      // throw std::runtime_error("Light type not implemented in Scene::computeIntersection(...)");
      }
    }
  }

  return result;
}

Ray computeReflectionRay(in Ray ray, in RaySceneIntersection intersection) {
  // https://en.wikipedia.org/wiki/Specular_reflection#Vector_formulation
  vec3 di = ray.dir;
  vec3 dn = intersection.normal;

  vec3 v_reflect = di - 2. * dn * dot(di, dn);
  normalize(v_reflect);
  Ray reflection_ray = newRay(intersection.intersection, v_reflect, ray);
  return reflection_ray;
}

Ray computeRefractionRay(in Ray ray, in RaySceneIntersection intersection) {
  vec2 in_out_mediums = getRayInOutIndexMediums(ray, intersection.material.index_medium, intersection.typeOfIntersectedObject, intersection.objectIndex);

  // https://amrhmorsy.github.io/blog/2024/RefractionVectorCalculation/
  vec3 L = ray.dir;
  vec3 N = intersection.normal;

  vec3 LparallelN = dot(N, L) * N;
  vec3 LperpendicularN = L - LparallelN;
  float sin_thetaL = length(LperpendicularN);

  // TODO: https://en.wikipedia.org/wiki/Fresnel_equations
  if(sin_thetaL <= in_out_mediums.y / in_out_mediums.x) {
    float r = in_out_mediums.x / in_out_mediums.y;
    float c = dot(N, L);
    vec3 T = N * (-r * c - sqrt(1. - r * r * (1. - c * c))) + L * r; // v_refract
    Ray refraction_ray = newRay(intersection.intersection, T, ray);
    return refraction_ray;
  } else {
    return computeReflectionRay(ray, intersection);
  }
}

float computeShadowIndex(in Ray ray, in RaySceneIntersection intersection, in Light light) {
  uint shadow_count = 0;
  for(uint i = 0; i < settings.Phong.SHADOW_RAYS; i++) {
    vec3 sampled_pos;
    if(settings.Phong.SHADOW_RAYS == 1) {
      sampled_pos = getLightCentralPos(light);
    } else {
      if(light.type == LightType_Spherical) {
        float theta = 2. * PI * random(vec4(screen_coords, i * 3));
        float phi = PI * random(vec4(screen_coords, i * 3 + 1));
        float r = light.sphere.m_radius * sqrt(random(vec4(screen_coords, i * 3 + 2)));
        sampled_pos = light.sphere.m_center + r * sphericalCoordinatesToEuclidean(theta, phi);
      } else {
        float u = random(vec4(screen_coords, i * 3));
        float r = random(vec4(screen_coords, i * 3 + 1));
        sampled_pos = light.quad.m_bottom_left + u * light.quad.m_up_vector + r * light.quad.m_right_vector;
      }
    }
    vec3 direction = intersection.intersection - sampled_pos;
    Ray shadow_ray = newRay(sampled_pos, direction, ray);
    RaySceneIntersection shadow_intersection = computeIntersection(shadow_ray, settings.EPSILON, length(direction) - settings.EPSILON, true);
    if(shadow_intersection.intersectionExists && shadow_intersection.typeOfIntersectedObject != LightIntersection) {
      shadow_count++;
    }
  }
  return float(shadow_count) / float(settings.Phong.SHADOW_RAYS);
}

vec3 phong(in Ray ray, in RaySceneIntersection intersection) {
  const Material material = intersection.material;
  const vec3 P = intersection.intersection;
  // const vec3 kd = settings.Bonus.ENABLE_TEXTURES && material.image_id >= 0 && !images[material.image_id].data.empty() ? images[material.image_id].getPixel(intersection.u, intersection.v) : material.diffuse_material;
  const vec3 kd = material.diffuse_material;

  if(!settings.Phong.ENABLED) {
    return kd;
  }

  // https://en.wikipedia.org/wiki/Phong_reflection_model#Concepts AND
  // https://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_reflection_model
  const vec3 V = -1. * ray.dir;
  const vec3 N = intersection.normal;

  vec3 Ia = vec3(0.);
  vec3 Id = vec3(0.);
  vec3 Is = vec3(0.);

  const vec3 ka = material.ambient_material;
  const vec3 ks = material.specular_material;
  const float alpha = material.shininess;

  for(uint light_i = 0; light_i < nb_lights; light_i++) {
    Light light = lights[light_i];

    // bidouillage 2
    float ia = 1.;
    float id = 1.;
    float is = 1.;

    vec3 light_pos = getLightCentralPos(light);
    vec3 L = normalize(light_pos - P);

    float shadow_index = settings.Phong.SHADOW_RAYS > 0 ? computeShadowIndex(ray, intersection, light) : 0.;

    float LdotN = dot(L, N);
    vec3 R = normalize(2. * N * LdotN - L);

    float RdotV = dot(R, V);
    Ia += ka * ia;
    Id += LdotN > settings.EPSILON ? (1. - shadow_index) * kd * LdotN * id : vec3(0.);
    Is += RdotV > settings.EPSILON ? (1. - shadow_index) * ks * pow(RdotV, alpha) * is : vec3(0.);
  }
  return Ia + Id + Is;
}

vec3 rayTraceIterative(in Ray _ray, in float _min_t, in float _max_t) {
  vec3 color = vec3(0.);
  uint bounces = 0;

  for(bounces = 0; bounces <= MAX_BOUNCES; bounces++) {
    RaySceneIntersection intersection = computeIntersection(_ray, settings.EPSILON, _max_t, false);
    if(!intersection.intersectionExists) {
      color = vec3(0.);
      break;
    }

    if(settings.Material.ENABLE_GLASS && intersection.material.type == Material_Glass) {
      _ray = computeRefractionRay(_ray, intersection);
    } else if(settings.Material.ENABLE_MIRROR && intersection.material.type == Material_Mirror) {
      _ray = computeReflectionRay(_ray, intersection);
    } else {
      color = max(vec3(0.), min(vec3(1.), phong(_ray, intersection)));
      break;
    }

    _min_t = settings.EPSILON;
  }

  return color;
  // return vec3(float(bounces) / 255.);
}

void main() {
  // ray origin
  dvec4 origin_temp = modelviewInverse * dvec4(0., 0., 0., 1.);
  vec3 origin = vec3(origin_temp.xyz / origin_temp.w);

  // ray direction
  vec2 randoms = vec2(random(vec4(screen_coords, 0)), random(vec4(screen_coords, 1)));
  vec2 uv = (screen_coords.xy + randoms) / vec2(settings.SCREEN_WIDTH, settings.SCREEN_HEIGHT);
  dvec4 direction_temp = modelviewInverse * projectionInverse * dvec4(2. * uv.x - 1., -(2. * uv.y - 1.), nearAndFarPlanes.x, 1.);
  vec3 direction = normalize(vec3(direction_temp.xyz / direction_temp.w) - origin);

  Ray ray = newRay(origin, direction);
  vec3 color = imageLoad(imgOutput, screen_coords.xy).rgb;
  color += rayTraceIterative(ray, min_t, max_t) / settings.NSAMPLES;
  imageStore(imgOutput, screen_coords.xy, vec4(color, 1.));
}