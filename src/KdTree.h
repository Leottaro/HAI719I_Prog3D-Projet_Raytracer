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
#include <memory>
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
        return v0[axis] - constants::kdtree::EPSILON <= pos ||
               v1[axis] - constants::kdtree::EPSILON <= pos ||
               v2[axis] - constants::kdtree::EPSILON <= pos;
    }

    bool isAfter(size_t axis, float pos) const {
        return v0[axis] + constants::kdtree::EPSILON >= pos ||
               v1[axis] + constants::kdtree::EPSILON >= pos ||
               v2[axis] + constants::kdtree::EPSILON >= pos;
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
        return min[0] - constants::kdtree::EPSILON < v[0] && min[1] - constants::kdtree::EPSILON < v[1] && min[2] - constants::kdtree::EPSILON < v[2] &&
               max[0] + constants::kdtree::EPSILON > v[0] && max[1] + constants::kdtree::EPSILON > v[1] && max[2] + constants::kdtree::EPSILON > v[2];
    }

    bool intersect(const Ray &ray, float &tmin, float &tmax) const {
        // https://www.rose-hulman.edu/class/cs/csse451/AABB/#:~:text=Axis%2DAligned%20Bounding%20Boxes%20(AABBs,bound%20and%20a%20maximum%20bound.
        Vec3 const &origin = ray.origin();
        Vec3 const &direction = ray.direction();

        tmin = (min[0] - origin[0]) / direction[0];
        tmax = (max[0] - origin[0]) / direction[0];
        if (tmin > tmax)
            std::swap(tmin, tmax);

        float tymin = (min[1] - origin[1]) / direction[1];
        float tymax = (max[1] - origin[1]) / direction[1];
        if (tymin > tymax)
            std::swap(tymin, tymax);

        if ((tmin - constants::kdtree::EPSILON > tymax) || (tymin - constants::kdtree::EPSILON > tmax))
            return false;
        tmin = std::max(tmin, tymin);
        tmax = std::min(tmax, tymax);

        float tzmin = (min[2] - origin[2]) / direction[2];
        float tzmax = (max[2] - origin[2]) / direction[2];
        if (tzmin > tzmax)
            std::swap(tzmin, tzmax);

        if ((tmin - constants::kdtree::EPSILON > tzmax) || (tzmin - constants::kdtree::EPSILON > tmax))
            return false;
        tmin = std::max(tmin, tzmin);
        tmax = std::min(tmax, tzmax);

        return true;
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

    std::unique_ptr<KdTree> left = nullptr;
    std::unique_ptr<KdTree> right = nullptr;

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
            this->left = make_unique<KdTree>(left_triangles, left_bounding_box, new_split_axis);
        }
        if (!right_triangles.empty()) {
            this->right = make_unique<KdTree>(right_triangles, right_bounding_box, new_split_axis);
        }
    }

    bool intersect(Ray const &ray, KdTriangle &res) const {
        float tmin = 0., tmax = 0.;
        if (!this->bounding_box.intersect(ray, tmin, tmax)) {
            return false;
        }

        if (this->is_leaf) {
            return this->intersect_leaf(ray, res);
        } else {
            return this->intersect_nonleaf(ray, res);
        }
    }

    bool intersect_leaf(Ray const &ray, KdTriangle &res) const {
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

    bool intersect_nonleaf(Ray const &ray, KdTriangle &res) const {
        KdTree *first = nullptr, *second = nullptr;
        float t_first_min = FLT_MAX, t_first_max = FLT_MAX;
        float t_second_min = FLT_MAX, t_second_max = FLT_MAX;

        // Check intersection with child bounding boxes
        if (this->left && this->left->bounding_box.intersect(ray, t_first_min, t_first_max)) {
            first = this->left.get();
        }
        if (this->right && this->right->bounding_box.intersect(ray, t_second_min, t_second_max)) {
            second = this->right.get();
        }

        // Determine traversal order
        if (t_first_min > t_second_min) {
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
