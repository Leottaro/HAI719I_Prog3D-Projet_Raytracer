#ifndef PLANE_H
#define PLANE_H
#include "Line.h"
#include "Vec3.h"
#include <cfloat>
class Plane {
private:
    Vec3 m_center, m_normal;

public:
    Plane() {}
    Plane(Vec3 const &c, Vec3 const &n) {
        m_center = c;
        m_normal = n;
        m_normal.normalize();
    }
    void setCenter(Vec3 const &c) { m_center = c; }
    void setNormal(Vec3 const &n) {
        m_normal = n;
        m_normal.normalize();
    }
    Vec3 const &center() const { return m_center; }
    Vec3 const &normal() const { return m_normal; }
    Vec3 project(Vec3 const &p) const {
        return p - squareDistance(p) * m_normal;
    }
    float squareDistance(Vec3 const &p) const { return (project(p) - p).squareLength(); }
    float distance(Vec3 const &p) const { return sqrt(squareDistance(p)); }
    bool isParallelTo(Line const &L) const {
        float dot = Vec3::dot(L.direction(), m_normal);
        return fabs(dot) <= Settings::EPSILON;
    }
    Vec3 getIntersectionPoint(Line const &L, float &t) const {
        if (isParallelTo(L)) {
            return Vec3();
        }
        const Vec3 &O = L.origin();
        const Vec3 &D = L.direction();
        const Vec3 &N = m_normal;
        t = -(Vec3::dot(N, O - m_center)) / Vec3::dot(N, D);
        return O + t * D;
    }
};
#endif
