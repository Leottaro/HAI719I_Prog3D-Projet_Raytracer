#ifndef SHADER_PROGRAM_H
#define SHADER_PROGRAM_H

#include <GL/glext.h>
#include <GLES2/gl2.h>
#include <memory>
#include <string>

class ShaderProgram {
public:
    /// Create the program. A valid OpenGL context must be active.
    ShaderProgram();

    virtual ~ShaderProgram();

    /// Generate a minimal shader program, made of one vertex shader and one fragment shader
    static std::shared_ptr<ShaderProgram> genComputeShaderProgram(const std::string &computeShaderFilename);

    /// OpenGL identifier of the program
    inline GLuint id() { return m_id; }

    /// Loads and compile a shader from a text file, before attaching it to a program
    void loadShader(GLenum type, const std::string &shaderFilename);

    /// The main GPU program is ready to be handle streams of polygons
    inline void link() { glLinkProgram(m_id); }

    /// Activate the program
    inline void use() { glUseProgram(m_id); }

    /// Desactivate the current program
    inline static void stop() { glUseProgram(0); }

    inline GLuint getLocation(const std::string &name) {
        return glGetUniformLocation(m_id, name.c_str());
    }

    inline void set(const std::string &name, int value) {
        glUniform1i(getLocation(name.c_str()), value);
    }

    inline void set(const std::string &name, GLuint value) {
        glUniform1i(getLocation(name.c_str()), value);
    }

    inline void set(const std::string &name, float value) {
        glUniform1f(getLocation(name.c_str()), value);
    }

    inline void set(const std::string &name, const Vec3 &value) {
        glUniform3fv(getLocation(name.c_str()), 1, value.valuePtr());
    }

    inline void set(const std::string &name, const Material &value) {
        set(name + ".image_id", value.image_id);
        set(name + ".ambient_material", value.ambient_material);
        set(name + ".diffuse_material", value.diffuse_material);
        set(name + ".specular_material", value.specular_material);
        set(name + ".shininess", (float)value.shininess);
        set(name + ".index_medium", value.index_medium);
        set(name + ".transparency", value.transparency);
        set(name + ".type", (int)value.type);
    }

    inline void set(const std::string &name, const Sphere &value) {
        set(name + ".m_center", value.m_center);
        set(name + ".m_radius", value.m_radius);
        set(name + ".material", value.material);
    }

    inline void set(const std::string &name, const Square &value) {
        set(name + ".m_normal", value.m_normal);
        set(name + ".m_bottom_left", value.m_bottom_left);
        set(name + ".m_right_vector", value.m_right_vector);
        set(name + ".m_up_vector", value.m_up_vector);
        set(name + ".material", value.material);
    }

    inline void set(const std::string &name, const Light &value) {
        set(name + ".material", value.material);
        set(name + ".isInCamSpace", value.isInCamSpace);
        set(name + ".type", value.type);
        set(name + ".sphere", value.sphere);
        set(name + ".quad", value.quad);
        set(name + ".powerCorrection", value.powerCorrection);
    }

    // inline void set(const std::string &name, const ppmLoader::ImageRGB &value) {
    // }

private:
    /// Loads the content of an ASCII file in a standard C++ string
    std::string file2String(const std::string &filename);

    GLuint m_id = 0;
};

#endif // SHADER_PROGRAM_H