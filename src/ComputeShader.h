#ifndef COMPUTE_SHADER_H
#define COMPUTE_SHADER_H

#include <GL/glext.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Vec3.h"
#include "Material.h"
#include "Sphere.h"
#include "Square.h"
#include "Scene.h"

struct ShaderMaterial {
    float ambient_material[3];
    float shininess;
    
    float diffuse_material[3];
    float index_medium;
    
    float specular_material[3];
    float transparency;
    
    int image_id;
    unsigned int type;
    float _pad0[2];
    
    ShaderMaterial() {}

    void assign(Material const &material) {
        shininess = material.shininess;
        index_medium = material.index_medium;
        transparency = material.transparency;
        image_id = material.image_id;
        type = static_cast<unsigned int>(material.type);
        for (unsigned char i = 0; i < 3; i++) {
            ambient_material[i] = material.ambient_material[i];
            diffuse_material[i] = material.diffuse_material[i];
            specular_material[i] = material.specular_material[i];
        }
    }
};

struct ShaderSphere {
    float m_center[3];
    float m_radius;
    ShaderMaterial material;

    ShaderSphere() {}

    void assign(Sphere const &sphere) {
        m_radius = sphere.m_radius;
        material.assign(sphere.material);
        for (unsigned char i = 0; i < 3; i++) {
            m_center[i] = sphere.m_center[i];
        }
    }
};

struct ShaderSquare {
    float m_normal[4]; // vec3 + padding
    float m_bottom_left[4]; // vec3 + padding
    float m_right_vector[4]; // vec3 + padding
    float m_up_vector[4]; // vec3 + padding
    ShaderMaterial material;

    ShaderSquare() {}

    void assign(Square const &square) {
        material.assign(square.material);
        for (unsigned char i = 0; i < 3; i++) {
            m_normal[i] = square.m_normal[i];
            m_bottom_left[i] = square.m_bottom_left[i];
            m_right_vector[i] = square.m_right_vector[i];
            m_up_vector[i] = square.m_up_vector[i];
        }
    }
};

struct ShaderLight {
    float material[3];
    float powerCorrection;

    unsigned int type;
    bool isInCamSpace;
    float _pad0[2];

    ShaderSphere sphere;
    ShaderSquare quad;

    ShaderLight() {}

    void assign(Light const &light) {
        powerCorrection = light.powerCorrection;
        type = light.type;
        isInCamSpace = light.isInCamSpace;
        sphere.assign(light.sphere);
        quad.assign(light.quad);
        for (unsigned char i = 0; i < 3; i++) {
            material[i] = light.material[i];
        }
    }
};

// https://learnopengl.com/Guest-Articles/2022/Compute-Shaders/Introduction
class ComputeShader {
private:
    unsigned int m_shader_id;
    int success;
    char infoLog[512];

public:
    // constructor reads and builds the shader
    ComputeShader(const std::string &_compute_path) {
        // load
        std::ifstream compute_file;
        compute_file.open(_compute_path);
        if (!compute_file.is_open()) {
            std::cerr << "ERROR::COMPUTE_SHADER::FILE_NOT_SUCCESFULLY_OPENED: " << _compute_path << std::endl;
            return;
        }

        std::stringstream compute_stream;
        compute_stream << compute_file.rdbuf();
        std::string compute_code_str = compute_stream.str();
        const char *compute_code = compute_code_str.c_str();
        compute_file.close();

        // compile
        unsigned int compute = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(compute, 1, &compute_code, NULL);
        glCompileShader(compute);
        glGetShaderiv(compute, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(compute, 512, NULL, infoLog);
            std::cerr << "ERROR::COMPUTE_SHADER::COMPILATION_FAILED\n"
                      << infoLog << std::endl;
        };

        m_shader_id = glCreateProgram();
        glAttachShader(m_shader_id, compute);
        glLinkProgram(m_shader_id);
        glGetProgramiv(m_shader_id, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(m_shader_id, 512, NULL, infoLog);
            std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n"
                      << infoLog << std::endl;
        }

        glDeleteShader(compute);
    }

    void use() {
        glUseProgram(m_shader_id);
    }

    void execute(GLuint num_groups_x, GLuint num_groups_y, GLuint num_groups_z) const {
        glDispatchCompute(num_groups_x, num_groups_y, num_groups_z);
    }

    void set(const std::string &name, GLboolean value) const {
        glUniform1i(glGetUniformLocation(m_shader_id, name.c_str()), (int)value);
    }

    void set(const std::string &name, GLint value) const {
        glUniform1i(glGetUniformLocation(m_shader_id, name.c_str()), value);
    }

    void set(const std::string &name, GLuint value) const {
        glUniform1ui(glGetUniformLocation(m_shader_id, name.c_str()), value);
    }

    void set(const std::string &name, GLfloat value) const {
        glUniform1f(glGetUniformLocation(m_shader_id, name.c_str()), value);
    }

    void set(const std::string &name, const Vec3 &value) const {
        glUniform3fv(glGetUniformLocation(m_shader_id, name.c_str()), 1, value.valuePtr());
    }

    void set(const std::string &name, const GLdouble v0, const GLdouble v1) const {
        glUniform2d(glGetUniformLocation(m_shader_id, name.c_str()), v0, v1);
    }

    void set(const std::string &name, const GLdouble value[]) const {
        glUniformMatrix4dv(glGetUniformLocation(m_shader_id, name.c_str()), 1, false, value);
    }

    void set(const std::string &name, const Material &value) const {
        set(name + ".image_id", value.image_id);
        set(name + ".ambient_material", value.ambient_material);
        set(name + ".diffuse_material", value.diffuse_material);
        set(name + ".specular_material", value.specular_material);
        set(name + ".shininess", (float)value.shininess);
        set(name + ".index_medium", value.index_medium);
        set(name + ".transparency", value.transparency);
        set(name + ".type", (GLuint)value.type);
    }

    void set(const std::string &name, const Sphere &value) const {
        set(name + ".m_center", value.m_center);
        set(name + ".m_radius", value.m_radius);
        set(name + ".material", value.material);
    }

    void set(const std::string &name, const Square &value) const {
        set(name + ".m_normal", value.m_normal);
        set(name + ".m_bottom_left", value.m_bottom_left);
        set(name + ".m_right_vector", value.m_right_vector);
        set(name + ".m_up_vector", value.m_up_vector);
        set(name + ".material", value.material);
    }

    void set(const std::string &name, const Light &value) const {
        set(name + ".material", value.material);
        set(name + ".isInCamSpace", value.isInCamSpace);
        set(name + ".type", (GLuint)value.type);
        set(name + ".sphere", value.sphere);
        set(name + ".quad", value.quad);
        set(name + ".powerCorrection", value.powerCorrection);
    }
};

#endif