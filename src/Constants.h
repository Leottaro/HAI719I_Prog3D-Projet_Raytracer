#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {

namespace general {
constexpr unsigned int NSAMPLES = 1;
// constexpr char MESH_PATH[] = "data/triangle.off";
constexpr char MESH_PATH[] = "data/nefertiti.off";
// constexpr char MESH_PATH[] = "data/unit_sphere_n.off";
constexpr int MAX_BOUNCES = 100;
} // namespace general

namespace phong {
constexpr bool ENABLED = true;
constexpr unsigned int SHADOW_RAYS = 1; // 16 default, 1 for hard shadows. Need PHONG enabled.
} // namespace phong

namespace materials {
constexpr bool ENABLE_MIRROR = true;
constexpr bool ENABLE_GLASS = true;
constexpr float AIR_INDEX_MEDIUM = 1.;
} // namespace materials

namespace kdtree {
constexpr unsigned int MAX_LEAF_SIZE = 32; // 0 to disable the kdtree, 32 default
}

} // namespace constants

#endif
