#ifndef RAY_H
#define RAY_H
#include "Line.h"
#include "Settings.h"

class Ray : public Line {
public:
    vector<float> index_mediums;
    vector<unsigned int> object_types;
    vector<unsigned int> object_indices;
    Ray() {}
    Ray(Vec3 const &o, Vec3 const &d) : Line(o, d), index_mediums({Settings::Material::AIR_INDEX_MEDIUM}), object_types({UINT32_MAX}), object_indices({UINT32_MAX}) {}
    Ray(Vec3 const &o, Vec3 const &d, Ray const &ray) : Line(o, d) {
        vector<float> im_copy = ray.index_mediums;
        vector<unsigned int> ot_copy = ray.object_types;
        vector<unsigned int> oi_copy = ray.object_indices;
        index_mediums = im_copy;
        object_types = ot_copy;
        object_indices = oi_copy;
    }
};
#endif