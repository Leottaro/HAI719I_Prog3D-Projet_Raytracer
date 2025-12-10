#ifndef CONSTANTS_H
#define CONSTANTS_H

// TODO: rayons proches toucheront les mêmes objets (les indexer selon la touche)
// TODO: profondeur de champ (considerer une lentille au lancement des pixels)

namespace constants {

namespace general {
constexpr float EPSILON = 1e-4;       // 1e-4 default
constexpr unsigned int NSAMPLES = 16; // 16 default
constexpr int MAX_BOUNCES = 1000;     // 100 default
} // namespace general

namespace phong {
constexpr bool ENABLED = true;           // true default
constexpr unsigned int SHADOW_RAYS = 16; // 16 default, 1 for hard shadows. Need PHONG enabled.
} // namespace phong

namespace materials {
constexpr bool ENABLE_MIRROR = true;   // true default
constexpr bool ENABLE_GLASS = true;    // true default
constexpr float AIR_INDEX_MEDIUM = 1.; // 1. default
} // namespace materials

namespace kdtree {
constexpr float EPSILON = 1e-8;            // 1e-8 default
constexpr unsigned int MAX_LEAF_SIZE = 32; // 32 default, 0 to disable the kdtree
} // namespace kdtree

namespace scenes {
constexpr char MESH_PATH[] = "data/nefertiti.off"; // "data/nefertiti.off" default
} // namespace scenes

} // namespace constants

#endif
