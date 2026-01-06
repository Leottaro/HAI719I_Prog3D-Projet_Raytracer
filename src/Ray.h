#ifndef RAY_H
#define RAY_H
#include "Line.h"
#include "Settings.h"

class Ray : public Line {
private:
    unsigned int nb_mediums;
    vector<float> index_mediums;
    vector<unsigned int> object_types;
    vector<unsigned int> object_indices;

public:
    Ray() {}
    Ray(Vec3 const &o, Vec3 const &d) : Line(o, d), nb_mediums(1), index_mediums({Settings::Material::AIR_INDEX_MEDIUM}), object_types({UINT32_MAX}), object_indices({UINT32_MAX}) {}
    Ray(Vec3 const &o, Vec3 const &d, Ray const &ray) : Line(o, d) {
        vector<float> im_copy = ray.index_mediums;
        vector<unsigned int> ot_copy = ray.object_types;
        vector<unsigned int> oi_copy = ray.object_indices;
        nb_mediums = ray.nb_mediums;
        index_mediums = im_copy;
        object_types = ot_copy;
        object_indices = oi_copy;
    }

    void getRayInOutIndexMediums(float index_medium, unsigned int object_type, unsigned int object_index, float &nL, float &nT) {
        if (object_types[nb_mediums - 1] == object_type && object_indices[nb_mediums - 1] == object_index) {
            nb_mediums--;
            nL = index_mediums[nb_mediums];
            nT = index_mediums[nb_mediums - 1];
            index_mediums.pop_back();
            object_types.pop_back();
            object_indices.pop_back();
        } else {
            index_mediums.push_back(index_medium);
            object_types.push_back(object_type);
            object_indices.push_back(object_index);
            nb_mediums++;
            nL = index_mediums[nb_mediums - 2];
            nT = index_mediums[nb_mediums - 1];
        }
        // return vec2(nL, nT);
    }
};
#endif