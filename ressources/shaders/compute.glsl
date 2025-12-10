#version 430 core

// struct Material {
//     int image_id;
//     vec3 ambient_material;
//     vec3 diffuse_material;
//     vec3 specular_material;
//     float shininess;

//     float index_medium;
//     float transparency;

//     int type;
// };

// struct Sphere {
//     vec3 m_center;
//     float m_radius;
//     Material material;
// };

// struct Square {
//     vec3 m_normal;
//     vec3 m_bottom_left;
//     vec3 m_right_vector;
//     vec3 m_up_vector;
//     Material material;
// };

// struct Light {
//     Vec3 material;
//     bool isInCamSpace;
//     int type;
//     Sphere sphere;
//     Square quad;
//     float powerCorrection;
// };

// const int Material_DiffUSE_PHONG = 0;
// const int Material_Glass = 1;
// const int Material_Mirror = 2;

// const int LightType_Spherical = 1;
// const int LightType_Quad = 2;

layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

uniform image2D imgOutput;
// uniform int nb_spheres;
// uniform Sphere spheres[10];
// uniform int nb_squares;
// uniform Square squares[10];
// const int Material_DiffUSE_PHONG = 0;
// const int Material_Glass = 1;
// const int Material_Mirror = 2;
// uniform int nb_lights;
// uniform Light lights[10];

void main() {
    vec4 value = vec4(0.0, 0.0, 0.0, 1.0);
    vec2 texelCoord = vec2(gl_GlobalInvocationID.xy);

    value.x = texelCoord.x / gl_NumWorkGroups.x;
    value.y = texelCoord.y / gl_NumWorkGroups.y;

    imageStore(imgOutput, texelCoord, value);
}