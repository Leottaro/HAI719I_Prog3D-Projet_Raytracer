#ifndef KDTREE_H
#define KDTREE_H

#include "Material.h"
#include "Ray.h"
#include "Triangle.h"
#include "Vec3.h"
#include <GL/glut.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

using namespace std;

class BoundingBox {
protected:
    Vec3 min;
    Vec3 max;
    // vector<Plane> planes;

public:
    BoundingBox() {}
    BoundingBox(Vec3 const &min, Vec3 const &max) : min(min), max(max) /*, planes({
                                                                             Plane(max, Vec3(1., 0., 0.)),
                                                                             Plane(max, Vec3(0., 1., 0.)),
                                                                             Plane(max, Vec3(0., 0., 1.)),
                                                                             Plane(min, Vec3(-1., 0., 0.)),
                                                                             Plane(min, Vec3(0., -1., 0.)),
                                                                             Plane(min, Vec3(0., 0., -1.)),
                                                                         })*/
    {}

    bool isInside(Vec3 const &v) const {
        return min[0] - 0.00001 < v[0] && min[1] - 0.00001 < v[1] && min[2] - 0.00001 < v[2] &&
               max[0] + 0.00001 > v[0] && max[1] + 0.00001 > v[1] && max[2] + 0.00001 > v[2];
    }

    bool intersect(Ray const &line) const {
        vector<Plane> planes = {
            Plane(max, Vec3(1., 0., 0.)),
            Plane(max, Vec3(0., 1., 0.)),
            Plane(max, Vec3(0., 0., 1.))};
        for (size_t i = 0; i < 3; i++) {
            if (planes[i].isParallelTo(line)) {
                continue;
            }
            float t;
            Vec3 intersection = planes[i].getIntersectionPoint(line, t);
            if (isInside(intersection)) {
                return true;
            }
        }
        return false;
    }
};

class KdTriangle {
public:
    Vec3 v0, v1, v2;
    size_t triangle_index;
    Vec3 centroid;

    KdTriangle() : v0(), v1(), v2(), triangle_index(), centroid() {}
    KdTriangle(Vec3 v0, Vec3 v1, Vec3 v2, size_t index) : v0(v0), v1(v1), v2(v2), triangle_index(index), centroid((v0 + v1 + v2) / 3.) {}

    bool isBefore(size_t axis, float pos) const {
        return v0[axis] - 0.00001 <= pos ||
               v1[axis] - 0.00001 <= pos ||
               v2[axis] - 0.00001 <= pos;
    }

    bool isAfter(size_t axis, float pos) const {
        return v0[axis] + 0.00001 >= pos ||
               v1[axis] + 0.00001 >= pos ||
               v2[axis] + 0.00001 >= pos;
    }
};

// -------------------------------------------
// KdTree
// -------------------------------------------

class KdTree {
protected:
    BoundingBox bounding_box;

    size_t split_axis;
    float split_pos;

    bool is_leaf;
    vector<KdTriangle> kd_triangles;

    KdTree *left = nullptr;
    KdTree *right = nullptr;

public:
    KdTree() {}
    KdTree(vector<KdTriangle> kd_triangles, Vec3 min, Vec3 max, unsigned short split_axis, size_t max_triangles_per_leaf) : bounding_box(min, max), split_axis(split_axis), is_leaf(kd_triangles.size() <= max_triangles_per_leaf) {
        assert(split_axis >= 0 && split_axis <= 2);
        assert(!kd_triangles.empty());
        for (const auto &triangle : kd_triangles) {
            assert(triangle.centroid[0] == triangle.centroid[0]); // Check for NaN
            assert(triangle.centroid[1] == triangle.centroid[1]);
            assert(triangle.centroid[2] == triangle.centroid[2]);
        }

        if (this->is_leaf) {
            this->kd_triangles = kd_triangles;
            return;
        }

        // sort the triangles and split by the median
        size_t function_split_axis = split_axis;
        assert(function_split_axis >= 0 && function_split_axis <= 2);
        auto sort_function = [function_split_axis](KdTriangle a, KdTriangle b) { return a.centroid[function_split_axis] < b.centroid[function_split_axis]; };
        stable_sort(kd_triangles.begin(), kd_triangles.end(), sort_function);
        size_t median_i = kd_triangles.size() / 2;
        this->split_pos = kd_triangles[median_i].centroid[split_axis];

        // split the Triangles
        vector<KdTriangle> left_triangles;
        vector<KdTriangle> right_triangles;
        for (const KdTriangle &triangle : kd_triangles) {
            if (triangle.isBefore(this->split_axis, this->split_pos)) {
                left_triangles.push_back(triangle);
            }
            if (triangle.isAfter(this->split_axis, this->split_pos)) {
                right_triangles.push_back(triangle);
            }
        }

        // If triangles are unsplittable don't split the node
        if (left_triangles.empty() || right_triangles.empty() || kd_triangles.size() == left_triangles.size() || kd_triangles.size() == right_triangles.size()) {
            this->is_leaf = true;
            this->kd_triangles = kd_triangles;
            return;
        }

        Vec3 new_min = Vec3(min);
        new_min[split_axis] = fmax(new_min[split_axis], this->split_pos);
        Vec3 new_max = Vec3(max);
        new_max[split_axis] = fmin(new_max[split_axis], this->split_pos);
        size_t new_split_axis = (split_axis + 1) % 3;

        if (!left_triangles.empty()) {
            this->left = new KdTree(left_triangles, min, new_max, new_split_axis, max_triangles_per_leaf);
        }
        if (!right_triangles.empty()) {
            this->right = new KdTree(right_triangles, new_min, max, new_split_axis, max_triangles_per_leaf);
        }
    }

    ~KdTree() {
        if (this->left != nullptr) {
            delete this->left;
            this->left = nullptr;
        }
        if (this->right != nullptr) {
            delete this->right;
            this->right = nullptr;
        }
    }

    KdTree(const KdTree &) = delete;
    KdTree &operator=(const KdTree &) = delete;

    bool intersect(Ray const &ray, KdTriangle &res) const {
        if (!this->bounding_box.intersect(ray)) {
            return false;
        }

        if (this->is_leaf) {
            float t = FLT_MAX;
            KdTriangle closest_kd_triangle;
            for (KdTriangle kd_triangle : kd_triangles) {
                Triangle triangle = Triangle(kd_triangle.v0, kd_triangle.v1, kd_triangle.v2);
                RayTriangleIntersection intersection = triangle.getIntersection(ray);
                if (intersection.intersectionExists && intersection.t < t) {
                    closest_kd_triangle = kd_triangle;
                    t = intersection.t;
                }
            }
            if (t == FLT_MAX) {
                return false;
            }
            res = closest_kd_triangle;
            return true;
        }

        KdTree *first = nullptr, *second = nullptr;
        if (ray.direction()[split_axis] > 0) {
            first = this->left;
            second = this->right;
        } else {
            first = this->right;
            second = this->left;
        }

        if (first != nullptr && first->intersect(ray, res)) {
            return true;
        }
        if (second != nullptr && second->intersect(ray, res)) {
            return true;
        }

        return false;
    }
};

#endif