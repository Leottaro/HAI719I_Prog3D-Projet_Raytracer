#ifndef Sphere_H
#define Sphere_H
#include "Mesh.h"
#include "Vec3.h"
#include <cmath>
#include <vector>

struct RaySphereIntersection {
    bool intersectionExists;
    float t;
    float theta, phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
    RaySphereIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

static Vec3 SphericalCoordinatesToEuclidean(float theta, float phi) {
    float sinPhi = sinf(phi);
    float x = sinPhi * sinf(theta);
    float y = cosf(phi);
    float z = sinPhi * cosf(theta);

    return Vec3(x, y, z);
}

static Vec3 EuclideanCoordinatesToSpherical(Vec3 xyz) {
    float R = xyz.length();
    float theta = atan2(xyz[0], xyz[2]); // azimuth around y-axis, 0..2π
    if (theta < 0.0f)
        theta += 2.0f * M_PI;

    float phi = acos(xyz[1] / R); // polar angle from +y axis, 0..π

    return Vec3(theta, phi, R);
}

class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c, float r) : Mesh(), m_center(c), m_radius(r) {}

    void build_arrays() {
        unsigned int nTheta = 20, nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi);
        normalsArray.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);
        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt) {
            float u = (float)(thetaIt) / (float)(nTheta - 1);
            float theta = u * 2 * M_PI;
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi - 1);
                float phi = v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean(theta, phi);
                positions_array[3 * vertexIndex + 0] = m_center[0] + m_radius * xyz[0];
                positions_array[3 * vertexIndex + 1] = m_center[1] + m_radius * xyz[1];
                positions_array[3 * vertexIndex + 2] = m_center[2] + m_radius * xyz[2];
                normalsArray[3 * vertexIndex + 0] = xyz[0];
                normalsArray[3 * vertexIndex + 1] = xyz[1];
                normalsArray[3 * vertexIndex + 2] = xyz[2];
                uvs_array[2 * vertexIndex + 0] = u;
                uvs_array[2 * vertexIndex + 1] = v;
            }
        }
        triangles_array.clear();
        for (unsigned int thetaIt = 0; thetaIt < nTheta - 1; ++thetaIt) {
            for (unsigned int phiIt = 0; phiIt < nPhi - 1; ++phiIt) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt + 1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt + 1) * nTheta;
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuV);
            }
        }
    }

    void translate(Vec3 const &translation) {
        m_center += translation;
        Mesh::translate(translation);
    }

    void apply_transformation_matrix(Mat3 transform) {
        m_center = transform * m_center;
        Mesh::apply_transformation_matrix(transform);
    }

    void scale(Vec3 const &scale) {
        Mat3 scale_matrix(scale[0], 0., 0.,
                          0., scale[1], 0.,
                          0., 0., scale[2]); // Matrice de transformation de mise à l'échelle
        apply_transformation_matrix(scale_matrix);
    }

    void rotate_x(float angle) {
        float x_angle = angle * M_PI / 180.;
        Mat3 x_rotation(1., 0., 0.,
                        0., cos(x_angle), -sin(x_angle),
                        0., sin(x_angle), cos(x_angle));
        apply_transformation_matrix(x_rotation);
    }

    void rotate_y(float angle) {
        float y_angle = angle * M_PI / 180.;
        Mat3 y_rotation(cos(y_angle), 0., sin(y_angle),
                        0., 1., 0.,
                        -sin(y_angle), 0., cos(y_angle));
        apply_transformation_matrix(y_rotation);
    }

    void rotate_z(float angle) {
        float z_angle = angle * M_PI / 180.;
        Mat3 z_rotation(cos(z_angle), -sin(z_angle), 0.,
                        sin(z_angle), cos(z_angle), 0.,
                        0., 0., 1.);
        apply_transformation_matrix(z_rotation);
    }

    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection = RaySphereIntersection();

        /*
        equation paramétrique du rayon : R(t) = O + tD    (avec O l'origine et D la direction)
        un point X est sur la sphère si ||X - C||^2 = r^2 (avec C le centre et r le radius)
        on sait que pour tout V, ||V||^2 = V·V donc l'équation de la sphère deviens (X-C)·(X-C) = r^2
        On remplace avec un point de la droite : (O+tD-C)·(O+tD-C) = r^2
        On développe : t^2 * ||D||^2 + t * 2(D·(O-C)) + ||O-C||^2 - r^2 = 0
        On a une équation du second degré avec a = ||D||^2 = 1 (car d normalisé), b = 2*(D·(O-C)) et c = ||O-C||^2 - r^2
        */

        const Vec3 &O = ray.origin();
        const Vec3 &D = ray.direction();
        const Vec3 &C = this->m_center;
        float r = this->m_radius;
        Vec3 CO = O - C;

        float a = D.squareNorm();
        float b = 2 * Vec3::dot(D, CO);
        float c = CO.squareNorm() - r * r;
        float delta = b * b - 4 * a * c;

        if (delta < 0) {
            return intersection;
        }

        float sqrt_delta = sqrt(delta);
        float t = (-b - sqrt_delta) / (2 * a);
        float t2 = (-b + sqrt_delta) / (2 * a);
        bool is_outside = t >= constants::general::EPSILON;
        // if (!is_outside) {
        //     return intersection;
        // }

        intersection.intersectionExists = true;
        intersection.t = is_outside ? t : t2;
        intersection.intersection = is_outside ? O + t * D : O + t2 * D;
        intersection.normal = is_outside ? intersection.intersection - C : C - intersection.intersection;
        intersection.secondintersection = is_outside ? O + t2 * D : O + t * D;

        intersection.normal.normalize();
        Vec3 spherical_pos = EuclideanCoordinatesToSpherical(intersection.normal);
        intersection.theta = spherical_pos[0];
        intersection.phi = spherical_pos[1];

        return intersection;
    }
};
#endif
