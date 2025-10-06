#ifndef SCENE_H
#define SCENE_H

#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include <GL/glut.h>
#include <string>
#include <vector>

using namespace std;

enum LightType {
    LightType_Spherical,
    LightType_Quad
};

struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}
};

struct RaySceneIntersection {
    bool intersectionExists;
    float t;
    Vec3 intersection;
    Vec3 normal;
    Material material;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection(float max_t = FLT_MAX) : intersectionExists(false), t(max_t) {}
};

class Scene {
    vector<Mesh> meshes;
    vector<Sphere> spheres;
    vector<Square> squares;
    vector<Light> lights;

public:
    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for (unsigned int It = 0; It < meshes.size(); ++It) {
            Mesh const &mesh = meshes[It];
            mesh.draw();
        }
        for (unsigned int It = 0; It < spheres.size(); ++It) {
            Sphere const &sphere = spheres[It];
            sphere.draw();
        }
        for (unsigned int It = 0; It < squares.size(); ++It) {
            Square const &square = squares[It];
            square.draw();
        }
    }

    RaySceneIntersection computeIntersection(Ray const &ray, float min_t, float max_t) {
        RaySceneIntersection result = RaySceneIntersection(max_t);

        for (unsigned int i = 0; i < spheres.size(); i++) {
            RaySphereIntersection intersection = spheres[i].intersect(ray);
            if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                result.intersectionExists = true;
                result.t = intersection.t;
                result.intersection = intersection.intersection;
                result.normal = intersection.normal;
                result.material = spheres[i].material;
                result.typeOfIntersectedObject = 0;
                result.objectIndex = i;
                result.raySphereIntersection = intersection;
            }
        }

        for (unsigned int i = 0; i < squares.size(); i++) {
            RaySquareIntersection intersection = squares[i].intersect(ray);
            if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                result.intersectionExists = true;
                result.t = intersection.t;
                result.intersection = intersection.intersection;
                result.normal = intersection.normal;
                result.material = squares[i].material;
                result.typeOfIntersectedObject = 1;
                result.objectIndex = i;
                result.raySquareIntersection = intersection;
            }
        }

        for (unsigned int i = 0; i < meshes.size(); i++) {
            RayTriangleIntersection intersection = meshes[i].intersect(ray);
            if (intersection.intersectionExists && min_t < intersection.t && intersection.t < result.t) {
                result.intersectionExists = true;
                result.t = intersection.t;
                result.intersection = intersection.intersection;
                result.normal = intersection.normal;
                result.material = meshes[i].material;
                result.typeOfIntersectedObject = 2;
                result.objectIndex = i;
                result.rayMeshIntersection = intersection;
            }
        }

        return result;
    }

    Vec3 rayTraceRecursive(Ray ray, float min_t, float max_t, int shadow_samples, int NRemainingBounces) {
        RaySceneIntersection raySceneIntersection = computeIntersection(ray, min_t, max_t);
        if (!raySceneIntersection.intersectionExists) {
            return Vec3(0., 0., 0.);
        }
        Material material = raySceneIntersection.material;
        Vec3 P = raySceneIntersection.intersection;

        if (material.type == Material_Glass && NRemainingBounces > 0) {
            float nL, nT;
            if (ray.object_types[ray.object_types.size() - 1] == raySceneIntersection.typeOfIntersectedObject && ray.object_indices[ray.object_indices.size() - 1] == raySceneIntersection.objectIndex) {
                nL = ray.index_mediums[ray.index_mediums.size() - 1];
                ray.index_mediums.pop_back();
                ray.object_types.pop_back();
                ray.object_indices.pop_back();
                nT = ray.index_mediums[ray.index_mediums.size() - 1];
            } else {
                nL = ray.index_mediums[ray.index_mediums.size() - 1];
                ray.index_mediums.push_back(material.index_medium);
                ray.object_types.push_back(raySceneIntersection.typeOfIntersectedObject);
                ray.object_indices.push_back(raySceneIntersection.objectIndex);
                nT = ray.index_mediums[ray.index_mediums.size() - 1];
            }

            // https://amrhmorsy.github.io/blog/2024/RefractionVectorCalculation/
            Vec3 L = ray.direction();
            Vec3 N = raySceneIntersection.normal;

            Vec3 LparallelN = Vec3::dot(N, L) * N;
            Vec3 LperpendicularN = L - LparallelN;
            float sin_thetaL = LperpendicularN.norm();

            if (sin_thetaL <= nT / nL) { // TODO: https://en.wikipedia.org/wiki/Fresnel_equations
                float r = nL / nT;
                float c = Vec3::dot(N, L);
                Vec3 T = N * (-r * c - sqrt(1 - r * r * (1 - c * c))) + L * r; // v_refract
                Ray refraction_ray = Ray(P, T, ray.index_mediums, ray.object_types, ray.object_indices);
                return rayTraceRecursive(refraction_ray, 0.00001, max_t, shadow_samples, NRemainingBounces - 1);
            } else {
                Vec3 v_reflect = L - 2. * N * Vec3::dot(L, N);
                v_reflect.normalize();
                Ray reflection_ray = Ray(P, v_reflect, ray.index_mediums, ray.object_types, ray.object_indices);
                return rayTraceRecursive(reflection_ray, 0.00001, max_t, shadow_samples, NRemainingBounces - 1);
            }
        }

        if (material.type == Material_Mirror && NRemainingBounces > 0) {
            // https://en.wikipedia.org/wiki/Specular_reflection#Vector_formulation
            Vec3 di = ray.direction();
            Vec3 dn = raySceneIntersection.normal;

            Vec3 v_reflect = di - 2. * dn * Vec3::dot(di, dn);
            v_reflect.normalize();
            Ray reflection_ray = Ray(P, v_reflect, ray.index_mediums, ray.object_types, ray.object_indices);
            return rayTraceRecursive(reflection_ray, 0.00001, max_t, shadow_samples, NRemainingBounces - 1);
        }

        // https://en.wikipedia.org/wiki/Phong_reflection_model#Concepts AND
        // https://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_reflection_model

        const Vec3 V = -1 * ray.direction();
        const Vec3 &N = raySceneIntersection.normal;

        Vec3 I = Vec3(0., 0., 0.);
        Vec3 ka = material.ambient_material;
        Vec3 kd = material.diffuse_material;
        Vec3 ks = material.specular_material;
        float alpha = material.shininess;

        float one_over_shadow_sample = 1. / (float)shadow_samples;

        for (Light light : lights) {
            // bidouillage 2
            float ia = 1.;
            float id = 1.;
            float is = 1.;

            Vec3 L = light.pos - P;
            float L_norm = L.norm();
            L /= L_norm;

            float shadow_factor = 1.;
            Vec3 arbitrary = Vec3();
            arbitrary[L.getMaxAbsoluteComponent()] = 1.;
            Vec3 u = Vec3::cross(arbitrary, L);
            Vec3 v = Vec3::cross(u, L);
            for (int i = 0; i < shadow_samples; i++) {
                float radius = light.radius * sqrt(((float)rand()) / RAND_MAX);
                float angle = 2 * M_PI * ((float)rand()) / RAND_MAX;
                Vec3 sampled_pos = light.pos + radius * cos(angle) * u + radius * sin(angle) * v;
                Ray shadow_ray = Ray(P, sampled_pos - P, ray.index_mediums, ray.object_types, ray.object_indices);
                RaySceneIntersection shadow_intersection = computeIntersection(shadow_ray, 0.00001, L_norm - 0.00001);
                if (shadow_intersection.intersectionExists && (shadow_intersection.typeOfIntersectedObject != raySceneIntersection.typeOfIntersectedObject || shadow_intersection.objectIndex != raySceneIntersection.objectIndex)) {
                    shadow_factor -= one_over_shadow_sample;
                }
            }

            float LdotN = Vec3::dot(L, N);
            Vec3 R = 2. * N * LdotN - L;
            R.normalize();

            float RdotV = Vec3::dot(R, V);
            Vec3 Ia = ka * ia;
            Vec3 Id = LdotN > 0.00001 ? kd * LdotN * id : Vec3();
            Vec3 Is = RdotV > 0.00001 ? ks * pow(RdotV, alpha) * is : Vec3();
            I += shadow_factor * (Ia + Id + Is);
        }

        return I;
    }

    Vec3 rayTrace(Ray const &rayStart, float min_t = 0.00001, float max_t = FLT_MAX, int shadow_samples = 10) {
        Vec3 color = rayTraceRecursive(rayStart, min_t, max_t, shadow_samples, 100);
        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        {
            spheres.resize(spheres.size() + 1);
            Sphere &s1 = spheres[spheres.size() - 1];
            s1.m_center = Vec3(1., 0., 0.);
            s1.m_radius = 1.f;
            s1.build_arrays();
            s1.material.type = Material_Diffuse_Blinn_Phong;
            s1.material.diffuse_material = Vec3(1., 0., 0.);
            s1.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s1.material.shininess = 20;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s2 = spheres[spheres.size() - 1];
            s2.m_center = Vec3(-1., 0., 0.);
            s2.m_radius = 1.f;
            s2.build_arrays();
            s2.material.type = Material_Diffuse_Blinn_Phong;
            s2.material.diffuse_material = Vec3(0., 1., 0);
            s2.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s2.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
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
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3(0.8, 0.8, 0.8);
            s.material.specular_material = Vec3(0.8, 0.8, 0.8);
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
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
            s.material.diffuse_material = Vec3(0., 0., 1.);
            s.material.specular_material = Vec3(0., 0., 1.);
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
            s.material.diffuse_material = Vec3(0., 1., 0.);
            s.material.specular_material = Vec3(0., 1., 0.);
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
            s.material.diffuse_material = Vec3(0., 1., 1.);
            s.material.specular_material = Vec3(0., 1., 1.);
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
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
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
            s.material.diffuse_material = Vec3(1., 0., 1.);
            s.material.specular_material = Vec3(1., 0., 1.);
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
            s.material.diffuse_material = Vec3(1.0, 1.0, 0.);
            s.material.specular_material = Vec3(1.0, 1.0, 0.);
            s.material.shininess = 16;
        }

        { // GLASS Sphere

            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
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
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }

    void setup_single_mesh(const string &filename) {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes[meshes.size() - 1];
            mesh.loadOFF(filename);
            mesh.scale(Vec3(2., 2., 2.));
            mesh.build_arrays();
            mesh.material.type = Material_Diffuse_Blinn_Phong;
            mesh.material.diffuse_material = Vec3(1., 0., 0.);
            mesh.material.specular_material = Vec3(0.2, 0.2, 0.2);
            mesh.material.shininess = 20;
        }
    }

    void setup_refraction_test() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        size_t width = 50;
        size_t depth = 50;
        squares.resize(width * depth);
        for (size_t z = 0; z < depth; z++) {
            for (size_t x = 0; x < width; x++) {
                size_t i = z * width + x;
                Vec3 color = (x + z) % 2 == 0 ? Vec3(0.46, 0.59, 0.34) : Vec3(0.93, 0.93, 0.82);

                Square &s = squares[i];
                s.setQuad(Vec3((float)x - width / 2., -1., (float)z - depth / 2.), Vec3(0., 0, 1.), Vec3(1., 0., 0.), 1., 1.);
                s.rotate_x(35);
                s.build_arrays();
                s.material.diffuse_material = Vec3(color);
                s.material.specular_material = Vec3(color);
                s.material.shininess = 16;
            }
        }

        { // Glass sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 0., 0.);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.ambient_material = Vec3(1., 1., 0.);
            s.material.diffuse_material = Vec3(1., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.;
            s.material.index_medium = 1.51;
        }
    }
};

#endif
