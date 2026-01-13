#ifndef SCENE_H
#define SCENE_H

#include "Mesh.h"
#include "Settings.h"
#include "Sphere.h"
#include "Square.h"
#include "imageLoader.h"
#include <GL/glut.h>
#include <string>
#include <vector>

using namespace std;

enum LightType {
    LightType_Spherical = 1,
    LightType_Quad = 2
};

struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Sphere sphere;
    Square quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

    Vec3 getCentralPos() const {
        switch (type) {
        case LightType_Quad:
            return quad.m_bottom_left + 0.5 * (quad.m_right_vector + quad.m_up_vector);
        case LightType_Spherical:
            return sphere.m_center;
        default:
            throw std::runtime_error("Light type not implemented in getPos()");
        }
    }
};

enum IntersectedObjectType {
    TriangleIntersection = 0,
    SphereIntersection = 1,
    SquareIntersection = 2,
    LightIntersection = 3,
};

struct RaySceneIntersection {
    // private:
    bool intersectionExists;
    float t, u, v;
    Vec3 intersection;
    Vec3 normal;
    Material material;
    IntersectedObjectType typeOfIntersectedObject;
    unsigned int objectIndex;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;

public:
    RaySceneIntersection(float max_t = FLT_MAX) : intersectionExists(false), t(max_t) {}

    void setRaySphereIntersection(const RaySphereIntersection &intersection, const Material &material, int object_i) {
        this->intersectionExists = true;
        this->t = intersection.t;
        this->u = intersection.theta / (2. * M_PI);
        this->v = intersection.phi / M_PI;
        this->intersection = intersection.intersection;
        this->material = material;
        this->normal = intersection.normal;
        this->typeOfIntersectedObject = SphereIntersection;
        this->objectIndex = object_i;
        this->raySphereIntersection = intersection;
    }
    void setRaySquareIntersection(const RaySquareIntersection &intersection, const Material &material, int object_i) {
        this->intersectionExists = true;
        this->t = intersection.t;
        this->u = intersection.u;
        this->v = intersection.v;
        this->intersection = intersection.intersection;
        this->normal = intersection.normal;
        this->material = material;
        this->typeOfIntersectedObject = SquareIntersection;
        this->objectIndex = object_i;
        this->raySquareIntersection = intersection;
    }
    void setRayTriangleIntersection(const RayTriangleIntersection &intersection, const Material &material, int object_i) {
        this->intersectionExists = true;
        this->t = intersection.t;
        this->u = intersection.u;
        this->v = intersection.v;
        this->intersection = intersection.intersection;
        this->normal = intersection.normal;
        this->material = material;
        this->typeOfIntersectedObject = TriangleIntersection;
        this->objectIndex = object_i;
        this->rayMeshIntersection = intersection;
    }
};

class Scene {
public:
    vector<Mesh> meshes;
    vector<Sphere> spheres;
    vector<Square> squares;
    vector<Light> lights;
    vector<ppmLoader::ImageRGB> images;

private:
    RaySceneIntersection computeIntersection(Ray const &ray, float min_t, float max_t, bool intersect_lights) {
        RaySceneIntersection result = RaySceneIntersection(max_t);

        for (unsigned int i = 0; i < spheres.size(); i++) {
            RaySphereIntersection intersection = spheres[i].intersect(ray);
            if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                result.setRaySphereIntersection(intersection, spheres[i].material, i);
            }
        }

        for (unsigned int i = 0; i < squares.size(); i++) {
            RaySquareIntersection intersection = squares[i].intersect(ray, false);
            if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                result.setRaySquareIntersection(intersection, squares[i].material, i);
            }
        }

        for (unsigned int i = 0; i < meshes.size(); i++) {
            RayTriangleIntersection intersection = meshes[i].intersect(ray);
            if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                result.setRayTriangleIntersection(intersection, meshes[i].material, i);
            }
        }

        if (intersect_lights) {
            for (unsigned int i = 0; i < lights.size(); i++) {
                switch (lights[i].type) {
                case LightType_Quad: {
                    RaySquareIntersection intersection = lights[i].quad.intersect(ray, true);
                    if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                        result.setRaySquareIntersection(intersection, lights[i].quad.material, i);
                        result.typeOfIntersectedObject = LightIntersection;
                    }
                } break;
                case LightType_Spherical: {
                    RaySphereIntersection intersection = lights[i].sphere.intersect(ray);
                    if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                        result.setRaySphereIntersection(intersection, lights[i].sphere.material, i);
                        result.typeOfIntersectedObject = LightIntersection;
                    }
                } break;
                default:
                    throw std::runtime_error("Light type not implemented in Scene::computeIntersection(...)");
                }
            }
        }

        return result;
    }

    Ray computeReflectionRay(Ray const &ray, RaySceneIntersection &intersection) {
        // https://en.wikipedia.org/wiki/Specular_reflection#Vector_formulation
        Vec3 di = ray.direction();
        Vec3 dn = intersection.normal;

        Vec3 v_reflect = di - 2. * dn * Vec3::dot(di, dn);
        v_reflect.normalize();
        Ray reflection_ray = Ray(intersection.intersection, v_reflect, ray);
        return reflection_ray;
    }

    Ray computeRefractionRay(Ray ray, RaySceneIntersection &intersection) {
        float nL, nT;
        ray.getRayInOutIndexMediums(intersection.material.index_medium, intersection.typeOfIntersectedObject, intersection.objectIndex, nL, nT);

        // https://amrhmorsy.github.io/blog/2024/RefractionVectorCalculation/
        Vec3 L = ray.direction();
        Vec3 N = intersection.normal;

        Vec3 LparallelN = Vec3::dot(N, L) * N;
        Vec3 LperpendicularN = L - LparallelN;
        float sin_thetaL = LperpendicularN.norm();

        // TODO: https://en.wikipedia.org/wiki/Fresnel_equations
        if (sin_thetaL <= nT / nL) {
            float r = nL / nT;
            float c = Vec3::dot(N, L);
            Vec3 T = N * (-r * c - sqrt(1 - r * r * (1 - c * c))) + L * r; // v_refract
            Ray refraction_ray = Ray(intersection.intersection, T, ray);
            return refraction_ray;
        } else {
            return computeReflectionRay(ray, intersection);
        }
    }

    float computeShadowIndex(const Ray &ray, const RaySceneIntersection &intersection, const Light &light) {
        int shadow_count = 0;
        for (unsigned int i = 0; i < Settings::Phong::SHADOW_RAYS; i++) {
            Vec3 sampled_pos;
            if (Settings::Phong::SHADOW_RAYS == 1) {
                sampled_pos = light.getCentralPos();
            } else {
                if (light.type == LightType_Spherical) {
                    float theta = 2 * M_PI * ((float)rand()) / RAND_MAX;
                    float phi = M_PI * ((float)rand()) / RAND_MAX;
                    float r = light.sphere.m_radius * sqrt(((float)rand()) / RAND_MAX);
                    sampled_pos = light.sphere.m_center + r * SphericalCoordinatesToEuclidean(theta, phi);
                } else {
                    float u = ((float)rand()) / RAND_MAX;
                    float r = ((float)rand()) / RAND_MAX;
                    sampled_pos = light.quad.m_bottom_left + u * light.quad.m_up_vector + r * light.quad.m_right_vector;
                }
            }
            Vec3 direction = intersection.intersection - sampled_pos;
            Ray shadow_ray = Ray(sampled_pos, direction, ray);
            RaySceneIntersection shadow_intersection = computeIntersection(shadow_ray, Settings::EPSILON, direction.length() - Settings::EPSILON, true);
            if (shadow_intersection.intersectionExists && shadow_intersection.typeOfIntersectedObject != LightIntersection) {
                shadow_count++;
            }
        }
        return (float)shadow_count / Settings::Phong::SHADOW_RAYS;
    }

    Vec3 phong(Ray const &ray, RaySceneIntersection const &intersection) {
        const Material &material = intersection.material;
        const Vec3 &P = intersection.intersection;
        const Vec3 kd = Settings::Bonus::ENABLE_TEXTURES && material.image_id >= 0 && !images[material.image_id].data.empty() ? images[material.image_id].getPixel(intersection.u, intersection.v) : material.diffuse_material;

        if (!Settings::Phong::ENABLED) {
            return kd;
        }

        // https://en.wikipedia.org/wiki/Phong_reflection_model#Concepts AND
        // https://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_reflection_model
        const Vec3 V = -1 * ray.direction();
        const Vec3 &N = intersection.normal;

        Vec3 Ia = Vec3();
        Vec3 Id = Vec3();
        Vec3 Is = Vec3();

        const Vec3 &ka = material.ambient_material;
        const Vec3 &ks = material.specular_material;
        const float &alpha = material.shininess;

        for (unsigned int light_i = 0; light_i < lights.size(); light_i++) {
            Light const &light = lights[light_i];

            // bidouillage 2
            float ia = 1.;
            float id = 1.;
            float is = 1.;

            Vec3 light_pos = light.getCentralPos();
            Vec3 L = light_pos - P;
            L.normalize();

            float shadow_index = Settings::Phong::SHADOW_RAYS > 0 ? computeShadowIndex(ray, intersection, light) : 0.;

            float LdotN = Vec3::dot(L, N);
            Vec3 R = 2. * N * LdotN - L;
            R.normalize();

            float RdotV = Vec3::dot(R, V);
            Ia += ka * ia;
            Id += LdotN > Settings::EPSILON ? (1. - shadow_index) * kd * LdotN * id : Vec3();
            Is += RdotV > Settings::EPSILON ? (1. - shadow_index) * ks * pow(RdotV, alpha) * is : Vec3();
        }
        return Ia + Id + Is;
    }

public:
    Scene() {}

    Vec3 rayTraceRecursive(Ray const &ray, float min_t, float max_t, int NRemainingBounces) {
        RaySceneIntersection raySceneIntersection = computeIntersection(ray, min_t, max_t, false);
        if (!raySceneIntersection.intersectionExists) {
            return Vec3();
        }

        Material material = raySceneIntersection.material;
        Vec3 color;

        if (Settings::Material::ENABLE_GLASS && material.type == Material_Glass && NRemainingBounces > 0) {
            Ray refraction = computeRefractionRay(ray, raySceneIntersection);
            color = rayTraceRecursive(refraction, Settings::EPSILON, max_t, NRemainingBounces - 1);
        } else if (Settings::Material::ENABLE_MIRROR && material.type == Material_Mirror && NRemainingBounces > 0) {
            Ray reflection = computeReflectionRay(ray, raySceneIntersection);
            color = rayTraceRecursive(reflection, Settings::EPSILON, max_t, NRemainingBounces - 1);
        } else {
            color = phong(ray, raySceneIntersection);
        }

        color[0] = min(max(0.f, color[0]), 1.f);
        color[1] = min(max(0.f, color[1]), 1.f);
        color[2] = min(max(0.f, color[2]), 1.f);

        return color;
    }

    Vec3 rayTraceIterative(Ray ray, float min_t, float max_t, int max_bounces) {
        Vec3 color;
        int bounces = 0;

        for (bounces = 0; bounces <= max_bounces; bounces++) {
            RaySceneIntersection intersection = computeIntersection(ray, min_t, max_t, false);
            if (!intersection.intersectionExists) {
                color = Vec3();
                break;
            }

            if (Settings::Material::ENABLE_GLASS && intersection.material.type == Material_Glass) {
                ray = computeRefractionRay(ray, intersection);
            } else if (Settings::Material::ENABLE_MIRROR && intersection.material.type == Material_Mirror) {
                ray = computeReflectionRay(ray, intersection);
            } else {
                color = phong(ray, intersection);

                color[0] = std::min(std::max(0.f, color[0]), 1.f);
                color[1] = std::min(std::max(0.f, color[1]), 1.f);
                color[2] = std::min(std::max(0.f, color[2]), 1.f);

                break;
            }

            min_t = Settings::EPSILON;
        }

        return color;
        // return Vec3(float(bounces) / 255.);
    }

    Vec3 rayTrace(Ray const &rayStart, float min_t = Settings::EPSILON, float max_t = FLT_MAX) {
        // Vec3 color = rayTraceRecursive(rayStart, min_t, max_t, Settings::MAX_BOUNCES);
        Vec3 color = rayTraceIterative(rayStart, min_t, max_t, Settings::MAX_BOUNCES);
        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        images.clear();

        {
            images.resize(images.size() + 1);
            ppmLoader::ImageRGB &image = images[images.size() - 1];
            ppmLoader::load_ppm(image, "ressources/img/sphereTextures/s1.ppm");
        }

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.type = LightType_Spherical;
            light.sphere = Sphere(Vec3(-5, 5, 5), 2.5f);
            light.powerCorrection = 2.f;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        {
            spheres.resize(spheres.size() + 1);
            Sphere &s1 = spheres[spheres.size() - 1];
            s1.m_center = Vec3();
            s1.m_radius = 1.f;
            s1.build_arrays();
            s1.material.type = Material_DiffUSE_PHONG;
            s1.material.diffuse_material = Vec3(1., 1., 1.);
            s1.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s1.material.shininess = 20;
            s1.material.image_id = 0;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s1 = spheres[spheres.size() - 1];
            s1.m_center = Vec3(3., 0., 0.);
            s1.m_radius = 1.f;
            s1.build_arrays();
            s1.material.type = Material_DiffUSE_PHONG;
            s1.material.diffuse_material = Vec3(1., 0., 0.);
            s1.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s1.material.shininess = 20;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s2 = spheres[spheres.size() - 1];
            s2.m_center = Vec3(0., 3., 0.);
            s2.m_radius = 1.f;
            s2.build_arrays();
            s2.material.type = Material_DiffUSE_PHONG;
            s2.material.diffuse_material = Vec3(0., 1., 0.);
            s2.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s2.material.shininess = 20;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s2 = spheres[spheres.size() - 1];
            s2.m_center = Vec3(0., 0., 3.);
            s2.m_radius = 1.f;
            s2.build_arrays();
            s2.material.type = Material_DiffUSE_PHONG;
            s2.material.diffuse_material = Vec3(0., 0., 1.);
            s2.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s2.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        images.clear();

        {
            images.resize(images.size() + 1);
            ppmLoader::ImageRGB &image = images[images.size() - 1];
            ppmLoader::load_ppm(image, "ressources/img/test/128.ppm");
        }

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.sphere = Sphere(Vec3(-5, 5, 5), 2.5f);
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.type = Material_DiffUSE_PHONG;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 20;
            s.material.image_id = 0;
        }
    }

    void setup_cornell_box() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        images.clear();

        {
            images.resize(images.size() + 1);
            ppmLoader::ImageRGB &image = images[images.size() - 1];
            ppmLoader::load_ppm(image, "ressources/img/test/128.ppm");
        }

        // { // Light quad
        //     squares.resize(squares.size() + 1);
        //     Square &s = squares[squares.size() - 1];
        //     s.setQuad(Vec3(-2, 1.5, -2), Vec3(1., 0, 0.), Vec3(0., 0., 1.), 4, 4);
        //     s.build_arrays();
        //     s.material.diffuse_material = Vec3(1., 0., 1.);
        //     s.material.specular_material = Vec3(1., 0., 1.);
        //     s.material.shininess = 16;
        // }

        { // Light
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.sphere = Sphere(Vec3(0.0, 1.5, 0.0), 2.5f);
            light.quad = Square(Vec3(-2, 1.5, -2), Vec3(1., 0, 0.), Vec3(0., 0., 1.), 4, 4);
            light.type = LightType_Quad;
            light.powerCorrection = 2.f;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.image_id = 0;
            s.material.diffuse_material = Vec3(1., 0., 1.);
            s.material.specular_material = Vec3(1., 0., 1.);
            s.material.shininess = 16;
        }

        { // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
        }

        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.5, 0.5, 0.5);
            s.material.specular_material = Vec3(0.5, 0.5, 0.5);
            s.material.shininess = 16;
        }

        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 0.);
            s.material.shininess = 16;
        }

        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0.5, 0.5);
            s.material.specular_material = Vec3(1., 0.5, 0.5);
            s.material.shininess = 16;
        }

        { // Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }

        { // GLASS Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.;
            s.material.index_medium = 1.51;
        }

        { // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0., 0., 1.);
            s.material.specular_material = Vec3(0., 0., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }

    void setup_single_mesh() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        images.clear();

        {
            images.resize(images.size() + 1);
            ppmLoader::ImageRGB &image = images[images.size() - 1];
            ppmLoader::load_ppm(image, "ressources/img/test/128.ppm");
        }

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.sphere = Sphere(Vec3(-5, 5, 5), 2.5f);
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF(Settings::availableMeshToPath(Settings::Mesh::MESH));
            mesh.rotate_y(45);
            mesh.translate(Vec3(-0.5, 0., 0.));
            mesh.scale(Vec3(2.));
            mesh.build_arrays();
            mesh.material.type = Material_DiffUSE_PHONG;
            mesh.material.diffuse_material = Vec3(1., 0., 0.);
            mesh.material.specular_material = Vec3(0.2, 0.2, 0.2);
            mesh.material.shininess = 20;
        }

        { // Back wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(10., 10., 1.));
            s.translate(Vec3(0., 0., -5));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 0.);
            s.material.shininess = 16;
            s.material.image_id = 0;
        }
    }

    void setup_refraction_test() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        images.clear();

        {
            images.resize(images.size() + 1);
            ppmLoader::ImageRGB &image = images[images.size() - 1];
            ppmLoader::load_ppm(image, "ressources/img/test/1024.ppm");
        }

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.sphere = Sphere(Vec3(-5, 5, 5), 2.5f);
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { // Floor
            double size = 50.;
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(size, size, 1.));
            s.translate(Vec3(0., 0., -1));
            s.rotate_x(-60);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 0.);
            s.material.shininess = 16;
            s.material.image_id = 0;
        }

        { // Glass sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3();
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.ambient_material = Vec3(1., 1., 0.);
            s.material.diffuse_material = Vec3(1., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.;
            s.material.index_medium = 1.;
        }
    }

    void setup_showcase() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        images.clear();

        {
            images.resize(images.size() + 1);
            ppmLoader::ImageRGB &image = images[images.size() - 1];
            ppmLoader::load_ppm(image, "ressources/img/test/128.ppm");
        }

        { // Light
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.sphere = Sphere(Vec3(0.0, 1.5, 0.0), 2.5f);
            light.quad = Square(Vec3(-2, 1.5, -2), Vec3(1., 0, 0.), Vec3(0., 0., 1.), 4, 4);
            light.type = LightType_Quad;
            light.powerCorrection = 2.f;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.image_id = 0;
            s.material.diffuse_material = Vec3(1., 0., 1.);
            s.material.specular_material = Vec3(1., 0., 1.);
            s.material.shininess = 16;
        }

        { // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
        }

        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.5, 0.5, 0.5);
            s.material.specular_material = Vec3(0.5, 0.5, 0.5);
            s.material.shininess = 16;
        }

        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 0.);
            s.material.shininess = 16;
        }

        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0.5, 0.5);
            s.material.specular_material = Vec3(1., 0.5, 0.5);
            s.material.shininess = 16;
        }

        { // Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }

        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF(Settings::availableMeshToPath(Settings::Mesh::MESH));
            mesh.scale(Vec3(3.));
            mesh.rotate_y(90);
            mesh.rotate_z(-45);
            mesh.translate(Vec3(-0.5, 0.5, -1.5));
            mesh.build_arrays();
            mesh.material.type = Material_DiffUSE_PHONG;
            mesh.material.diffuse_material = Vec3(1., 0., 0.);
            mesh.material.specular_material = Vec3(0.2, 0.2, 0.2);
            mesh.material.shininess = 20;
        }

        { // GLASS Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-0.5, 0., 0.);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.;
            s.material.index_medium = 11. / 10.;
        }

        { // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.75);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0., 0., 1.);
            s.material.specular_material = Vec3(0., 0., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }
};

#endif
