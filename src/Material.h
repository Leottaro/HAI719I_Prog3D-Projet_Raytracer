#ifndef MATERIAL_H
#define MATERIAL_H

#include "Vec3.h"
#include "imageLoader.h"
#include <GL/glut.h>
#include <cmath>

enum MaterialType {
    Material_DiffUSE_PHONG = 1,
    Material_Glass = 2,
    Material_Mirror = 3
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
        diffuse_material = Vec3(.5);
        specular_material = Vec3(1.);
        ambient_material = Vec3();
        shininess = 16.;
    }
};

#endif // MATERIAL_H