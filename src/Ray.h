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
    Ray(Vec3 const &o, Vec3 const &d, vector<float> const &im, vector<unsigned int> const &ot, vector<unsigned int> const &oi) : Line(o, d) {
        vector<float> im_copy = im;
        vector<unsigned int> ot_copy = ot;
        vector<unsigned int> oi_copy = oi;
        index_mediums = im_copy;
        object_types = ot_copy;
        object_indices = oi_copy;
    }
};
#endif