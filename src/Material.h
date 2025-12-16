#ifndef MATERIAL_H
#define MATERIAL_H

#include "Vec3.h"
#include "imageLoader.h"
#include <GL/glut.h>
#include <cmath>

enum MaterialType {
    Material_DiffUSE_PHONG,
    Material_Glass,
    Material_Mirror
};

struct Material {
    int image_id;
    Vec3 ambient_material;
    Vec3 diffuse_material;
    Vec3 specular_material;
    double shininess;

    float index_medium;
    float transparency;

    MaterialType type;

    Material() {
        type = Material_DiffUSE_PHONG;
        image_id = -1;
        transparency = 0.0;
        index_medium = Settings::Material::AIR_INDEX_MEDIUM;
        ambient_material = Vec3();
    }
};

#endif // MATERIAL_H