#version 430 core
const float PI = 3.14159265359;

// ======================================================================================================
// ============================================= STRUCTURES =============================================
// ======================================================================================================

// Pour uniforms

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

  PhongSettings Phong;
  MaterialSettings Material;
  MeshSettings Mesh;
  KdTreeSettings KdTree;
  BonusSettings Bonus;
};

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
const uint Material_DiffUSE_PHONG = 0;
const uint Material_Glass = 1;
const uint Material_Mirror = 2;

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

struct Light {
  vec3 material;
  bool isInCamSpace;
  uint type;
  Sphere sphere;
  Square quad;
  float powerCorrection;
};
const uint LightType_Spherical = 1;
const uint LightType_Quad = 2;

// Pour calculs

struct Ray {
  vec3 pos;
  vec3 dir;
  uint nb_bounces;
  float index_mediums[10];
  uint object_types[10];
  uint object_indices[10];
};

struct RaySphereIntersection {
  bool intersectionExists;
  float t;
  float theta, phi;
  vec3 intersection;
  vec3 secondintersection;
  vec3 normal;
};
struct RaySquareIntersection {
  bool intersectionExists;
  float t;
  float u;
  float v;
  vec3 intersection;
  vec3 normal;
};
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
};
// const uint TriangleIntersection = 0;
const uint SphereIntersection = 1;
const uint SquareIntersection = 2;
const uint LightIntersection = 3;

// ======================================================================================================
// =============================================== INPUTS ===============================================
// ======================================================================================================

layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

layout(rgba32f, binding = 0) uniform image2D imgOutput;
uniform uint screen_width;
uniform uint screen_height;
uniform sampler2D randomTexture; // the texture has a width of (screen_width * NSAMPLES) and a height of (screen_height * (SHADOW_RAYS * 3 + 1))

uniform Settings settings;
uniform uint nb_spheres;
uniform Sphere spheres[10];
uniform uint nb_squares;
uniform Square squares[10];
uniform uint nb_lights;
uniform Light lights[10];
// uniform uint nb_meshes;
// uniform Mesh meshes;

uniform mat4 model_view;
uniform mat4 model_view_inverse;
uniform mat4 projection;
uniform mat4 projection_inverse;
uniform vec2 near_and_far_planes;

// =====================================================================================================
// ============================================= FUNCTIONS =============================================
// =====================================================================================================

// HELPERS

float sampleRandomValue(uint shadow_ray_i, uint sub_shadow_ray_i) {
  float y = gl_GlobalInvocationID.y * (settings.Phong.SHADOW_RAYS * 3 + 1) + 3 * shadow_ray_i + sub_shadow_ray_i;
  float x = gl_GlobalInvocationID.x * settings.NSAMPLES + gl_GlobalInvocationID.z;
  vec2 uv = vec2(x, y) / vec2(screen_width * settings.NSAMPLES, screen_height * (settings.Phong.SHADOW_RAYS + 1));
  return texture(randomTexture, uv).r;
}

float random2d(vec2 coord) {
  return fract(sin(dot(coord.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 cameraSpaceToWorldSpace(vec3 pos) {
  vec4 res = model_view_inverse * vec4(pos, 1.);
  return res.xyz / res.w;
}

vec3 screenSpaceToWorldSpace(vec2 uv) {
  // u et v sont entre 0 et 1 (0,0 est en haut a gauche de l'ecran)
  vec4 res = model_view_inverse * projection_inverse * vec4(2. * uv.x - 1, 2. * uv.y - 1, near_and_far_planes.x, 1.);
  return res.xyz / res.w;
}

vec3 sphericalCoordinatesToEuclidean(float theta, float phi) {
  float sinPhi = sin(phi);
  float x = sinPhi * sin(theta);
  float y = cos(phi);
  float z = sinPhi * cos(theta);

  return vec3(x, y, z);
}

vec3 euclideanCoordinatesToSpherical(vec3 pos) {
  float R = length(pos);
  float theta = atan(pos.x, pos.z); // azimuth around y-axis, 0..2π // TODO: atan2 ??
  if(theta < 0.0) {
    theta += 2.0 * PI;
  }

  float phi = acos(pos.y / R); // polar angle from +y axis, 0..π

  return vec3(theta, phi, R);
}

vec3 getLightCentralPos(Light light) {
  if(light.type == LightType_Quad)
    return light.quad.m_bottom_left + 0.5 * (light.quad.m_right_vector + light.quad.m_up_vector);
  if(light.type == LightType_Spherical)
    return light.sphere.m_center;
}

// INTERSECTIONS

void setRaySphereIntersection(inout RaySceneIntersection result, RaySphereIntersection intersection, Material material, uint object_i) {
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
void setRaySquareIntersection(inout RaySceneIntersection result, RaySquareIntersection intersection, Material material, uint object_i) {
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

RaySquareIntersection intersectSquare(Ray ray, Square square, bool can_intersect_behind) {
  RaySquareIntersection intersection = RaySquareIntersection(false, 0., 0., 0., vec3(0.), vec3(0.));

  const vec3 O = ray.pos;
  const vec3 D = ray.dir;
  const vec3 C = square.m_bottom_left;
  const vec3 N = square.m_normal;

  // nous sommes derrière le carré donc pas d'intersection
  if(!can_intersect_behind && square.material.type != Material_Glass && dot(D, N) >= 0) {
    return intersection;
  }

  float t = -(dot(N, O - C)) / dot(N, D);

  vec3 P = O + t * D;
  vec3 PC = P - C;
  float u = dot(PC, square.m_right_vector) / dot(square.m_right_vector, square.m_right_vector);
  float v = dot(PC, square.m_up_vector) / dot(square.m_up_vector, square.m_up_vector);

  if(settings.EPSILON <= u && u <= 1 && settings.EPSILON <= v && v <= 1) {
    intersection.intersectionExists = true;
    intersection.t = t;
    intersection.u = u;
    intersection.v = v;
    intersection.intersection = P;
    intersection.normal = N;
  }

  return intersection;
}

RaySphereIntersection intersectSphere(Ray ray, Sphere sphere) {
  RaySphereIntersection intersection = RaySphereIntersection(false, 0., 0., 0., vec3(0.), vec3(0.), vec3(0.));

  const vec3 O = ray.pos;
  const vec3 D = ray.dir;
  const vec3 C = sphere.m_center;
  float r = sphere.m_radius;
  vec3 CO = O - C;

  float a = dot(D, D);
  float b = 2 * dot(D, CO);
  float c = dot(CO, CO) - r * r;
  float delta = b * b - 4 * a * c;

  if(delta < 0) {
    return intersection;
  }

  float sqrt_delta = sqrt(delta);
  float t = (-b - sqrt_delta) / (2 * a);
  float t2 = (-b + sqrt_delta) / (2 * a);
  bool is_outside = t >= settings.EPSILON;
  // if (!is_outside) {
  //     return intersection;
  // }

  intersection.intersectionExists = true;
  intersection.t = is_outside ? t : t2;
  intersection.intersection = is_outside ? O + t * D : O + t2 * D;
  intersection.normal = normalize(is_outside ? intersection.intersection - C : C - intersection.intersection);
  intersection.secondintersection = is_outside ? O + t2 * D : O + t * D;

  vec3 spherical_pos = euclideanCoordinatesToSpherical(intersection.normal);
  intersection.theta = spherical_pos[0];
  intersection.phi = spherical_pos[1];

  return intersection;
}

// SCENE.H

RaySceneIntersection computeIntersection(Ray ray, float min_t, float max_t, bool intersect_lights) {
  RaySceneIntersection result;
  result.intersectionExists = false;
  result.t = max_t;

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

  // for(uint i = 0; i < nb_meshes; i++) {
  //   RayTriangleIntersection intersection = meshes[i].intersect(ray);
  //   if(intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
  //     setRayTriangleIntersection(result, intersection, meshes[i].material, i);
  //   }
  // }

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

Ray computeReflectionRay(Ray ray, RaySceneIntersection intersection) {
  // https://en.wikipedia.org/wiki/Specular_reflection#Vector_formulation
  const vec3 di = ray.dir;
  const vec3 dn = intersection.normal;

  vec3 v_reflect = normalize(di - 2. * dn * dot(di, dn));
  Ray reflection_ray = Ray(intersection.intersection, v_reflect, ray.nb_bounces, ray.index_mediums, ray.object_types, ray.object_indices);
  return reflection_ray;
}

Ray computeRefractionRay(Ray ray, RaySceneIntersection intersection) {
  float nL, nT;
  if(ray.object_types[ray.nb_bounces - 1] == intersection.typeOfIntersectedObject && ray.object_indices[ray.nb_bounces - 1] == intersection.objectIndex) {
    nL = ray.index_mediums[ray.nb_bounces - 1];
    ray.nb_bounces -= 1;
    nT = ray.index_mediums[ray.nb_bounces - 1];
  } else {
    nL = ray.index_mediums[ray.nb_bounces - 1];
    ray.index_mediums[ray.nb_bounces] = intersection.material.index_medium;
    ray.object_types[ray.nb_bounces] = intersection.typeOfIntersectedObject;
    ray.object_indices[ray.nb_bounces] = intersection.objectIndex;
    ray.nb_bounces += 1;
    nT = ray.index_mediums[ray.nb_bounces - 1];
  }

  // https://amrhmorsy.github.io/blog/2024/RefractionVectorCalculation/
  const vec3 L = ray.dir;
  const vec3 N = intersection.normal;

  vec3 LparallelN = dot(N, L) * N;
  vec3 LperpendicularN = L - LparallelN;
  float sin_thetaL = length(LperpendicularN);

  // TODO: https://en.wikipedia.org/wiki/Fresnel_equations
  if(sin_thetaL <= nT / nL) {
    float r = nL / nT;
    float c = dot(N, L);
    vec3 T = N * (-r * c - sqrt(1 - r * r * (1 - c * c))) + L * r; // v_refract
    Ray refraction_ray = Ray(intersection.intersection, T, ray.nb_bounces, ray.index_mediums, ray.object_types, ray.object_indices);
    return refraction_ray;
  } else {
    return computeReflectionRay(ray, intersection);
  }
}

float computeShadowIndex(Ray ray, RaySceneIntersection intersection, Light light) {
  int shadow_count = 0;
  for(uint i = 0; i < settings.Phong.SHADOW_RAYS; i++) {
    vec3 sampled_pos;
    if(settings.Phong.SHADOW_RAYS == 1) {
      sampled_pos = getLightCentralPos(light);
    } else {
      if(light.type == LightType_Spherical) {
        float theta = 2 * PI * sampleRandomValue(i, 0);
        float phi = PI * sampleRandomValue(i, 1);
        float r = light.sphere.m_radius * sampleRandomValue(i, 2);
        sampled_pos = light.sphere.m_center + r * sphericalCoordinatesToEuclidean(theta, phi);
      } else {
        float u = sampleRandomValue(i, 0);
        float r = sampleRandomValue(i, 1);
        sampled_pos = light.quad.m_bottom_left + u * light.quad.m_up_vector + r * light.quad.m_right_vector;
      }
    }
    vec3 direction = intersection.intersection - sampled_pos;
    Ray shadow_ray = Ray(sampled_pos, direction, ray.nb_bounces, ray.index_mediums, ray.object_types, ray.object_indices);
    RaySceneIntersection shadow_intersection = computeIntersection(shadow_ray, settings.EPSILON, length(direction) - settings.EPSILON, true);
    if(shadow_intersection.intersectionExists && shadow_intersection.typeOfIntersectedObject != LightIntersection) {
      shadow_count++;
    }
  }
  return float(shadow_count) / settings.Phong.SHADOW_RAYS;
}

vec3 phong(Ray ray, RaySceneIntersection intersection, uint NRemainingBounces) {
  const Material material = intersection.material;
  const vec3 P = intersection.intersection;
  // const vec3 kd = settings.Bonus.ENABLE_TEXTURES && material.image_id >= 0 && !images[material.image_id].data.empty() ? images[material.image_id].getPixel(intersection.u, intersection.v) : material.diffuse_material;
  const vec3 kd = material.diffuse_material;

  if(!settings.Phong.ENABLED) {
    return kd;
  }

  // https://en.wikipedia.org/wiki/Phong_reflection_model#Concepts AND
  // https://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_reflection_model
  const vec3 V = -1 * ray.dir;
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

vec3 rayTraceIterative(Ray ray, float min_t, float max_t) {
  vec3 color = vec3(0.); // Accumulated color
  uint bounces = 0;

  while(bounces < settings.MAX_BOUNCES) {
    RaySceneIntersection intersection = computeIntersection(ray, min_t, max_t, false);
    if(!intersection.intersectionExists) {
      break; // No intersection, terminate the loop
    }

    Material material = intersection.material;
    color += phong(ray, intersection, settings.MAX_BOUNCES - bounces);

    if(settings.Material.ENABLE_GLASS && material.type == Material_Glass) {
      ray = computeRefractionRay(ray, intersection);
    } else if(settings.Material.ENABLE_MIRROR && material.type == Material_Mirror) {
      ray = computeReflectionRay(ray, intersection);
    } else {
      break; // No further bounces
    }

    bounces++;
  }

  return color;
}

vec3 rayTrace(Ray rayStart, float min_t, float max_t) {
  return rayTraceIterative(rayStart, min_t, max_t);
}

void main() {
  ivec2 screenCoords = ivec2(gl_GlobalInvocationID.xy);
  vec2 uv = vec2(gl_GlobalInvocationID.xy) / gl_NumWorkGroups.xy;
  // float randomColor = texture(randomTexture, uv).r;
  // imageStore(imgOutput, screenCoords, vec4(uv, randomColor, 1.));

  vec3 camera_pos = cameraSpaceToWorldSpace(vec3(0.));
  vec3 direction = normalize(screenSpaceToWorldSpace(uv) - camera_pos); // TODO: nsamples

  Ray ray;
  ray.pos = camera_pos;
  ray.dir = direction;
  ray.nb_bounces = 0;
  vec3 color = rayTrace(ray, 0., 1000000.);
  imageStore(imgOutput, screenCoords, vec4(color, 1.));
}