#ifndef SETTINGS_H
#define SETTINGS_H

#include <ostream>
#include <string>

// TODO: rayons proches toucheront les mêmes objets (les indexer selon la touche)
// TODO: profondeur de champ (considerer une lentille au lancement des pixels)

class Settings {
public:
    static constexpr unsigned int NB_PRESETS = 11;
    enum class Presets {
        PHASE_1,
        PHASE_2_NO_SHADOWS,
        PHASE_2_HARD_SHADOWS,
        PHASE_2_SOFT_SHADOWS,
        PHASE_3_REFLECTION,
        PHASE_3_NO_INTERPOLATION,
        PHASE_3_INTERPOLATION,
        PHASE_4_REFRACTION,
        PHASE_4_KDTREE,
        PICTURE,
        BENCHMARKS,
    };

    enum class AvailableMeshes {
        NONE = 0,            // 0 triangles
        TRIANGLE = 1,        // 1 triangle
        NEFERTITI = 562,     // 562 triangles
        UNIT_SPHERE_N = 722, // 722 triangles
        KNOT = 4160,         // 4160 triangles
        FEMUR = 7798,        // 7798 triangles
        PEGASO = 30658,      // 30658 triangles
        FEMALE01 = 368369,   // 368369 triangles
    };
    static constexpr AvailableMeshes AVAILABLE_MESHES[8] = {AvailableMeshes::NONE, AvailableMeshes::TRIANGLE, AvailableMeshes::NEFERTITI, AvailableMeshes::UNIT_SPHERE_N, AvailableMeshes::KNOT, AvailableMeshes::FEMUR, AvailableMeshes::PEGASO, AvailableMeshes::FEMALE01};
    inline static std::string availableMeshToPath(Settings::AvailableMeshes mesh) {
        switch (mesh) {
        case Settings::AvailableMeshes::FEMALE01:
            return "ressources/models/female01.off";
        case Settings::AvailableMeshes::FEMUR:
            return "ressources/models/femur.off";
        case Settings::AvailableMeshes::KNOT:
            return "ressources/models/knot.off";
        case Settings::AvailableMeshes::NEFERTITI:
            return "ressources/models/nefertiti.off";
        case Settings::AvailableMeshes::PEGASO:
            return "ressources/models/pegaso.off";
        case Settings::AvailableMeshes::TRIANGLE:
            return "ressources/models/triangle.off";
        case Settings::AvailableMeshes::UNIT_SPHERE_N:
            return "ressources/models/unit_sphere_n.off";
        default:
            return "";
        }
    }

    inline static unsigned int selected_renderer;
    inline static Presets selected_preset;

    static constexpr float EPSILON = 1e-4;
    inline static unsigned int NSAMPLES;
    inline static unsigned int MAX_BOUNCES;
    inline static unsigned int SCREEN_WIDTH = 480;
    inline static unsigned int SCREEN_HEIGHT = 480;

    class Phong {
    public:
        inline static bool ENABLED;
        inline static unsigned int SHADOW_RAYS;
    };

    class Material {
    public:
        inline static bool ENABLE_MIRROR;
        inline static bool ENABLE_GLASS;
        inline static float AIR_INDEX_MEDIUM;
    };

    class Mesh {
    public:
        inline static AvailableMeshes MESH;
        inline static bool ENABLE_INTERPOLATION;
    };

    class KdTree {
    public:
        static constexpr float EPSILON = 1e-8;
        inline static unsigned int MAX_LEAF_SIZE;
    };

    class Bonus {
    public:
        inline static bool ENABLE_TEXTURES;
    };

    inline static void applyPreset() {
        switch (selected_preset) {
        case Presets::PHASE_1:
            selected_renderer = 2;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NONE;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = false;
            Phong::SHADOW_RAYS = 16;
            Material::ENABLE_MIRROR = false;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;

        case Presets::PHASE_2_NO_SHADOWS:
            selected_renderer = 2;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NONE;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 0;
            Material::ENABLE_MIRROR = false;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_2_HARD_SHADOWS:
            selected_renderer = 2;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NONE;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 1;
            Material::ENABLE_MIRROR = false;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_2_SOFT_SHADOWS:
            selected_renderer = 2;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NONE;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 16;
            Material::ENABLE_MIRROR = false;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_3_REFLECTION:
            selected_renderer = 2;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NONE;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 16;
            Material::ENABLE_MIRROR = true;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_3_NO_INTERPOLATION:
            selected_renderer = 3;
            NSAMPLES = 4;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NEFERTITI;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 4;
            Material::ENABLE_MIRROR = true;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_3_INTERPOLATION:
            selected_renderer = 3;
            NSAMPLES = 4;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NEFERTITI;
            Mesh::ENABLE_INTERPOLATION = true;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 4;
            Material::ENABLE_MIRROR = true;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_4_REFRACTION:
            selected_renderer = 4;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::NEFERTITI;
            Mesh::ENABLE_INTERPOLATION = true;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 16;
            Material::ENABLE_MIRROR = true;
            Material::ENABLE_GLASS = true;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        case Presets::PHASE_4_KDTREE:
            selected_renderer = 3;
            NSAMPLES = 16;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::KNOT;
            Mesh::ENABLE_INTERPOLATION = true;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 16;
            Material::ENABLE_MIRROR = true;
            Material::ENABLE_GLASS = true;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 16;
            Bonus::ENABLE_TEXTURES = false;
            break;

        case Presets::PICTURE:
            selected_renderer = 2;
            NSAMPLES = 16;
            MAX_BOUNCES = 10000;
            Mesh::MESH = AvailableMeshes::PEGASO;
            Mesh::ENABLE_INTERPOLATION = true;
            Phong::ENABLED = true;
            Phong::SHADOW_RAYS = 32;
            Material::ENABLE_MIRROR = true;
            Material::ENABLE_GLASS = true;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 16;
            Bonus::ENABLE_TEXTURES = true;
            break;

        case Presets::BENCHMARKS:
            NSAMPLES = 1;
            MAX_BOUNCES = 100;
            Mesh::MESH = AvailableMeshes::TRIANGLE;
            Mesh::ENABLE_INTERPOLATION = false;
            Phong::ENABLED = false;
            Phong::SHADOW_RAYS = 0;
            Material::ENABLE_MIRROR = false;
            Material::ENABLE_GLASS = false;
            Material::AIR_INDEX_MEDIUM = 1.0;
            KdTree::MAX_LEAF_SIZE = 0;
            Bonus::ENABLE_TEXTURES = false;
            break;
        }
    }

    inline static void applyPreset(Presets preset) {
        selected_preset = preset;
        applyPreset();
    }
};

inline std::ostream &operator<<(std::ostream &os, const Settings::Presets &preset) {
    switch (preset) {
    case Settings::Presets::PHASE_1:
        os << "PHASE_1";
        break;
    case Settings::Presets::PHASE_2_NO_SHADOWS:
        os << "PHASE_2_NO_SHADOWS";
        break;
    case Settings::Presets::PHASE_2_HARD_SHADOWS:
        os << "PHASE_2_HARD_SHADOWS";
        break;
    case Settings::Presets::PHASE_2_SOFT_SHADOWS:
        os << "PHASE_2_SOFT_SHADOWS";
        break;
    case Settings::Presets::PHASE_3_REFLECTION:
        os << "PHASE_3_REFLECTION";
        break;
    case Settings::Presets::PHASE_3_NO_INTERPOLATION:
        os << "PHASE_3_NO_INTERPOLATION";
        break;
    case Settings::Presets::PHASE_3_INTERPOLATION:
        os << "PHASE_3_INTERPOLATION";
        break;
    case Settings::Presets::PHASE_4_REFRACTION:
        os << "PHASE_4_REFRACTION";
        break;
    case Settings::Presets::PHASE_4_KDTREE:
        os << "PHASE_4_KDTREE";
        break;
    case Settings::Presets::PICTURE:
        os << "PICTURE";
        break;
    case Settings::Presets::BENCHMARKS:
        os << "BENCHMARKS";
        break;
    }
    return os;
}

inline std::ostream &operator<<(std::ostream &os, const Settings::AvailableMeshes &mesh) {
    switch (mesh) {
    case Settings::AvailableMeshes::NONE:
        os << "NONE";
        break;
    case Settings::AvailableMeshes::TRIANGLE:
        os << "TRIANGLE";
        break;
    case Settings::AvailableMeshes::NEFERTITI:
        os << "NEFERTITI";
        break;
    case Settings::AvailableMeshes::UNIT_SPHERE_N:
        os << "UNIT_SPHERE_N";
        break;
    case Settings::AvailableMeshes::KNOT:
        os << "KNOT";
        break;
    case Settings::AvailableMeshes::FEMUR:
        os << "FEMUR";
        break;
    case Settings::AvailableMeshes::PEGASO:
        os << "PEGASO";
        break;
    case Settings::AvailableMeshes::FEMALE01:
        os << "FEMALE01";
        break;
    }
    return os;
}

#endif