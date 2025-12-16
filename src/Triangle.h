#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Plane.h"
#include "Ray.h"
#include "Vec3.h"
#include <cfloat>

struct RayTriangleIntersection {
    bool intersectionExists;
    float t;
    float w0, w1, w2;
    float u, v;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
    RayTriangleIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

class Triangle {
private:
    Vec3 m_normal;
    float area;

public:
    Vec3 m_c[3];

    Triangle() {}
    Triangle(Vec3 const &c0, Vec3 const &c1, Vec3 const &c2) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]);
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    Vec3 const &normal() const { return m_normal; }

    Vec3 getCentroid() const {
        return (m_c[0] + m_c[1] + m_c[2]) / 3.;
    }

    Vec3 projectOnSupportPlane(Vec3 const &p) const {
        return Plane(m_c[0], m_normal).project(p);
    }

    float squareDistanceToSupportPlane(Vec3 const &p) const {
        return Plane(m_c[0], m_normal).squareDistance(p);
    }
    float distanceToSupportPlane(Vec3 const &p) const {
        return Plane(m_c[0], m_normal).distance(p);
    }

    bool isParallelTo(Line const &L) const {
        return Plane(m_c[0], m_normal).isParallelTo(L);
    }

    Vec3 getIntersectionPointWithSupportPlane(Line const &L, float &t) const {
        return Plane(m_c[0], m_normal).getIntersectionPoint(L, t);
    }

    static float computeArea(Vec3 const &p0, Vec3 const &p1, Vec3 const &p2) {
        Vec3 crossProduct = Vec3::cross(p1 - p0, p2 - p0);
        return 0.5f * crossProduct.length();
    }

    void computeBarycentricCoordinates(Vec3 const &p, float &u0, float &u1, float &u2) const {
        u0 = Triangle(p, m_c[1], m_c[2]).area / area - Settings::EPSILON;
        u1 = Triangle(m_c[0], p, m_c[2]).area / area - Settings::EPSILON;
        u2 = Triangle(m_c[0], m_c[1], p).area / area - Settings::EPSILON;
    }

    RayTriangleIntersection getIntersection(Ray const &ray) const {
        RayTriangleIntersection intersection;
        intersection.intersectionExists = false;

        // 1) check that the ray is not parallel to the triangle
        if (this->isParallelTo(ray)) {
            return intersection;
        }

        // calculate the intersection
        float t = 0;
        Vec3 P = getIntersectionPointWithSupportPlane(ray, t);

        // 2) check that the triangle is "in front of" the ray
        if (Vec3::dot(P - ray.origin(), ray.direction()) < 0) {
            return intersection;
        }

        // 3) check that the intersection point is inside the triangle:
        float w0, w1, w2;
        computeBarycentricCoordinates(P, w0, w1, w2);
        if (w0 < 0 || 1 < w0 || w1 < 0 || 1 < w1 || w2 < 0 || 1 < w2 || w0 + w1 + w2 < 0 || 1 < w0 + w1 + w2) {
            return intersection;
        }

        // 4) Finally, if all conditions were met, then there is an intersection!
        intersection.intersectionExists = true;
        intersection.t = t;
        intersection.w0 = w0;
        intersection.w1 = w1;
        intersection.w2 = w2;
        // intersection.u;
        // intersection.v;
        // intersection.tIndex;
        intersection.intersection = P;
        // intersection.normal;
        return intersection;
    }
};
#endif
