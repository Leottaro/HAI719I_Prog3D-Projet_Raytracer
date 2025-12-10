#ifndef RENDERER_H
#define RENDERER_H

#include "Scene.h"
#include "ShaderProgram.h"
#include "Vec3.h"
#include "matrixUtilities.h"
#include <GL/glew.h>
#include <GL/glext.h>
#include <chrono>
#include <iomanip>
#include <string>
#include <vector>

struct Renderer : Scene {
public:
    Renderer() {}

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for (unsigned int i = 0; i < meshes.size(); ++i) {
            Mesh const &mesh = meshes[i];
            mesh.draw();
        }
        for (unsigned int i = 0; i < spheres.size(); ++i) {
            Sphere const &sphere = spheres[i];
            sphere.draw();
        }
        for (unsigned int i = 0; i < squares.size(); ++i) {
            Square const &square = squares[i];
            square.draw();
        }
    }

    vector<Vec3> rayTraceFromCameraCPU(int _width, int _height, int _max_t) {
        cout << "Ray tracing a " << _width << " x " << _height << "image :" << endl;
        Vec3 pos, dir;
        int n_pixels = _width * _height;
        vector<Vec3> image(n_pixels, Vec3(0, 0, 0));

        auto begin = chrono::high_resolution_clock::now();

        for (int y = 0; y < _height; y++) {
            for (int x = 0; x < _width; x++) {
                int pixel_i = y * _width + x + 1;

                float percent = (float)pixel_i / n_pixels;
                int64_t currently_elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count();
                float remaining_ms = (currently_elapsed / percent) * (1. - percent);

                cout << "\r\tCalculating pixel " << pixel_i << " of " << n_pixels << " (" << fixed << setprecision(2) << 100.f * percent << "% completed) ~" << remaining_ms / 1000. << "s remaining     " << flush;
                for (unsigned int s = 0; s < constants::general::NSAMPLES; ++s) {
                    float u = ((float)(x) + (float)(rand()) / (float)(RAND_MAX)) / _width;
                    float v = ((float)(y) + (float)(rand()) / (float)(RAND_MAX)) / _height;
                    // this is a random uv that belongs to the pixel xy.
                    screen_space_to_world_space_ray(u, v, pos, dir);
                    Vec3 color = rayTrace(Ray(pos, dir), 0., _max_t);
                    image[x + y * _width] += color;
                }
                image[x + y * _width] /= (float)constants::general::NSAMPLES;
            }
        }
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(end - begin).count();
        cout << endl
             << "\tDone in " << elapsed << " seconds." << endl;
        return image;
    }

    // =====================================================================================================
    // ================================================ GPU ================================================
    // =====================================================================================================
private:
    std::shared_ptr<ShaderProgram> m_shaderProgramPtr;
    GLuint m_texture_out;
    unsigned int m_width;
    unsigned int m_height;

    void free() {
        glDeleteTextures(1, &m_texture_out);
        m_texture_out = 0;
        m_width = 0;
        m_height = 0;
    }

public:
    bool init(unsigned int _width, unsigned int _height) {
        if (m_width != 0 || m_height != 0) {
            free();
        }
        m_width = _width;
        m_height = _height;

        glGenTextures(1, &m_texture_out);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_texture_out);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, m_width, m_height, 0, GL_RGBA, GL_FLOAT, NULL);
    }

    void updateScene() {
        int nb_spheres = spheres.size();
        m_shaderProgramPtr->set("nb_spheres", nb_spheres);
        for (unsigned int i = 0; i < nb_spheres; i++) {
            m_shaderProgramPtr->set("spheres[" + std::to_string(i) + "]", nb_spheres);
        }
        int nb_squares = squares.size();
        m_shaderProgramPtr->set("nb_squares", nb_squares);
        for (unsigned int i = 0; i < nb_squares; i++) {
            m_shaderProgramPtr->set("squares[" + std::to_string(i) + "]", nb_squares);
        }
        int nb_meshes = meshes.size();
        m_shaderProgramPtr->set("nb_meshes", nb_meshes);
        for (unsigned int i = 0; i < nb_meshes; i++) {
            m_shaderProgramPtr->set("meshes[" + std::to_string(i) + "]", nb_meshes);
        }
        int nb_lights = lights.size();
        m_shaderProgramPtr->set("nb_lights", nb_lights);
        for (unsigned int i = 0; i < nb_lights; i++) {
            m_shaderProgramPtr->set("lights[" + std::to_string(i) + "]", nb_lights);
        }
    }

    void execute(int _max_t) {
        m_shaderProgramPtr->use();
        updateScene();
        m_shaderProgramPtr->stop();
    }

    vector<Vec3> getImageOut() {
        glBindTexture(GL_TEXTURE_2D, m_texture_out);
        std::vector<Vec3> image(m_width * m_height);

        std::vector<float> textureData(m_width * m_height * 4); // RGBA format
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, textureData.data());

        for (int y = 0; y < m_height; ++y) {
            for (int x = 0; x < m_width; ++x) {
                int index = (y * m_width + x) * 4; // RGBA has 4 components
                image[y * m_width + x] = Vec3(textureData[index], textureData[index + 1], textureData[index + 2]);
            }
        }

        return image;
    }

    vector<Vec3> rayTraceFromCameraGPU(int _width, int _height, int _max_t) {
        init(_width, _height);
        updateScene();
        execute(_max_t);
        return getImageOut();
    }
};

#endif