#ifndef SQUARE_H
#define SQUARE_H
#include "Mesh.h"
#include "Vec3.h"
#include <cmath>
#include <vector>

struct RaySquareIntersection {
    bool intersectionExists;
    float t;
    float u, v;
    Vec3 intersection;
    Vec3 normal;
    RaySquareIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
           float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
                 float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector, upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector * width;
        m_up_vector = m_up_vector * height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;
        vertices[0].u = uMin;
        vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;
        vertices[1].u = uMax;
        vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;
        vertices[2].u = uMax;
        vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;
        vertices[3].u = uMin;
        vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;
    }

    void translate(Vec3 const &translation) {
        m_bottom_left += translation;
        Mesh::translate(translation);
    }

    void apply_transformation_matrix(Mat3 transform) {
        m_normal = transform * m_normal;
        m_bottom_left = transform * m_bottom_left;
        m_right_vector = transform * m_right_vector;
        m_up_vector = transform * m_up_vector;
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

    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersection = RaySquareIntersection();

        /*
        equation paramétrique du rayon : R(t) = O + tD      (avec O l'origine et D la direction)
        Un point X est sur le plan si : N·(X-C) = 0         (avec C un point du plan et N sa normale)
        On remplace par le point du rayon : N·(O + tD - C) = 0
        On réarrange bien : t = - N·(O-C) / N·D
        On a P = O+tD l'intersection du rayon au plan.
        Notons qu'un vecteur A peut se décomposer en deux vecteurs u et v de cette manière: A = u * (A·U/U·U) + v * (A·V/V·V)
        On peut donc retrouver les u et v de notre point P : u = PC·U/U·U et v = PC·V/V·V
        Le point P est dans le rectangle il se décompose en un u et v dans [0; 1]
        */

        const Vec3 &O = ray.origin();
        const Vec3 &D = ray.direction();
        const Vec3 &C = this->m_bottom_left;
        const Vec3 &N = this->m_normal;

        // nous sommes derrière le carré donc pas d'intersection
        if (material.type != Material_Glass && Vec3::dot(D, N) > 0) {
            return intersection;
        }

        float t = -(Vec3::dot(N, O - C)) / Vec3::dot(N, D);

        Vec3 P = O + t * D;
        Vec3 PC = P - C;
        float u = Vec3::dot(PC, this->m_right_vector) / this->m_right_vector.squareNorm();
        float v = Vec3::dot(PC, this->m_up_vector) / this->m_up_vector.squareNorm();

        if (0 <= u && u <= 1 && 0 <= v && v <= 1) {
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.u = u;
            intersection.v = v;
            intersection.intersection = P;
            intersection.normal = N;
        }

        return intersection;
    }
};
#endif