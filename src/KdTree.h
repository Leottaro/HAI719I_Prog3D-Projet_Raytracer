#ifndef KDTREE_H
#define KDTREE_H

#include "Constants.h"
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

class BoundingBox {
protected:
    // vector<Plane> planes;

public:
    Vec3 min;
    Vec3 max;

    BoundingBox() {}
    BoundingBox(Vec3 const &min, Vec3 const &max) : min(min), max(max) {}

    bool isInside(Vec3 const &v) const {
        return min[0] - 0.00001 < v[0] && min[1] - 0.00001 < v[1] && min[2] - 0.00001 < v[2] &&
               max[0] + 0.00001 > v[0] && max[1] + 0.00001 > v[1] && max[2] + 0.00001 > v[2];
    }

    bool intersect(Line const &line, float &min_t) const {
        vector<Plane> planes = {
            Plane(max, Vec3(1., 0., 0.)),
            Plane(max, Vec3(0., 1., 0.)),
            Plane(max, Vec3(0., 0., 1.)),
        };
        min_t = FLT_MAX;
        for (size_t i = 0; i < 3; i++) {
            if (planes[i].isParallelTo(line)) {
                continue;
            }
            float t;
            planes[i].getIntersectionPoint(line, t);
            if (t < min_t) {
                min_t = t;
            }
        }
        return min_t != FLT_MAX;
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
    KdTree(vector<KdTriangle> kd_triangles, BoundingBox bounding_box, unsigned short split_axis) {
        this->bounding_box = bounding_box;
        this->split_axis = split_axis;
        this->is_leaf = kd_triangles.size() <= constants::kdtree::MAX_LEAF_SIZE;

        assert(constants::kdtree::MAX_LEAF_SIZE > 0);
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

        // Generate the bounding boxes
        Vec3 left_max = Vec3(bounding_box.max);
        left_max[split_axis] = this->split_pos;
        Vec3 right_min = Vec3(bounding_box.min);
        right_min[split_axis] = this->split_pos;
        BoundingBox left_bounding_box = BoundingBox(bounding_box.min, left_max);
        BoundingBox right_bounding_box = BoundingBox(right_min, bounding_box.max);

        // split the Triangles
        vector<KdTriangle> left_triangles;
        vector<KdTriangle> right_triangles;
        for (const KdTriangle &triangle : kd_triangles) {
            bool is_before = triangle.isBefore(this->split_axis, this->split_pos);
            bool is_after = triangle.isAfter(this->split_axis, this->split_pos);
            if (is_before) {
                left_triangles.push_back(triangle);
            }
            if (is_after) {
                right_triangles.push_back(triangle);
            }
        }

        // If triangles are unsplittable don't split the node
        if (kd_triangles.size() == left_triangles.size() || kd_triangles.size() == right_triangles.size()) {
            this->is_leaf = true;
            this->kd_triangles = kd_triangles;
            return;
        }

        size_t new_split_axis = (split_axis + 1) % 3;
        if (!left_triangles.empty()) {
            this->left = new KdTree(left_triangles, left_bounding_box, new_split_axis);
        }
        if (!right_triangles.empty()) {
            this->right = new KdTree(right_triangles, right_bounding_box, new_split_axis);
        }
    }

    // ~KdTree() {
    //     if (this->is_leaf) {
    //         return;
    //     }

    //     if (this->left != nullptr) {
    //         delete this->left;
    //         this->left = nullptr;
    //     }
    //     if (this->right != nullptr) {
    //         delete this->right;
    //         this->right = nullptr;
    //     }
    // }

    bool intersect(Ray const &ray, KdTriangle &res) const {
        float tttttttt;
        if (!this->bounding_box.intersect(ray, tttttttt)) {
            return false;
        }

        if (this->is_leaf) {
            float closest_t = FLT_MAX;
            KdTriangle closest_kd_triangle;
            for (const KdTriangle &kd_triangle : this->kd_triangles) {
                Triangle triangle = Triangle(kd_triangle.v0, kd_triangle.v1, kd_triangle.v2);
                RayTriangleIntersection intersection = triangle.getIntersection(ray);
                if (intersection.intersectionExists && intersection.t < closest_t) {
                    closest_kd_triangle = KdTriangle(kd_triangle);
                    closest_t = intersection.t;
                }
            }
            if (closest_t == FLT_MAX) {
                return false;
            }
            res = closest_kd_triangle;
            return true;
        }

        KdTree *first = nullptr, *second = nullptr;
        float t_first = FLT_MAX, t_second = FLT_MAX;

        // Calculate distances to child bounding boxes
        if (this->left && this->left->bounding_box.intersect(ray, t_first)) {
            first = this->left;
        }
        if (this->right && this->right->bounding_box.intersect(ray, t_second)) {
            second = this->right;
        }

        // Order children based on distance
        if (t_first > t_second) {
            std::swap(first, second);
        }

        // Traverse children in order
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