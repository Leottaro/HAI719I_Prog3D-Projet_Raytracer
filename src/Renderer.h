#ifndef RENDERER_H
#define RENDERER_H

#include <GL/glew.h>
#include <GL/glext.h>
#include <GL/gl.h>

#include "Scene.h"
#include "Vec3.h"
#include "matrixUtilities.h"
#include "ComputeShader.h"

#include <chrono>
#include <iomanip>
#include <string>
#include <vector>
#include <cstring>

struct Renderer : Scene {
public:
    // =====================================================================================================
    // ================================================ CPU ================================================
    // =====================================================================================================

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
        cout << "Ray tracing a " << _width << " x " << _height << " image (CPU) :" << endl;
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
                for (unsigned int s = 0; s < Settings::NSAMPLES; ++s) {
                    float u = ((float)(x) + (float)(rand()) / (float)(RAND_MAX)) / _width;
                    float v = ((float)(y) + (float)(rand()) / (float)(RAND_MAX)) / _height;
                    // this is a random uv that belongs to the pixel xy.
                    screen_space_to_world_space_ray(u, v, pos, dir);
                    Vec3 color = rayTrace(Ray(pos, dir), 0., _max_t);
                    image[x + y * _width] += color;
                }
                image[x + y * _width] /= (float)Settings::NSAMPLES;
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
    ComputeShader m_shader;

    void updateSceneUniforms() {
        unsigned int nb_spheres = spheres.size();
        m_shader.set("nb_spheres", nb_spheres);
        for (unsigned int i = 0; i < nb_spheres; i++) {
            m_shader.set("spheres[" + std::to_string(i) + "]", spheres[i]);
        }
        unsigned int nb_squares = squares.size();
        m_shader.set("nb_squares", nb_squares);
        for (unsigned int i = 0; i < nb_squares; i++) {
            m_shader.set("squares[" + std::to_string(i) + "]", squares[i]);
        }
        // unsigned int nb_meshes = meshes.size();
        // m_shader.set("nb_meshes", nb_meshes);
        // for (unsigned int i = 0; i < nb_meshes; i++) {
        //     m_shader.set("meshes[" + std::to_string(i) + "]", meshes[i]);
        // }
        unsigned int nb_lights = lights.size();
        m_shader.set("nb_lights", nb_lights);
        for (unsigned int i = 0; i < nb_lights; i++) {
            m_shader.set("lights[" + std::to_string(i) + "]", lights[i]);
        }
    }

public:
    vector<Vec3> rayTraceFromCameraGPU(int _width, int _height, int _max_t) {
        cout << "Ray tracing a " << _width << " x " << _height << " image (GPU) :" << endl;
        auto begin = chrono::high_resolution_clock::now();

        GLuint texture;
        glGenTextures(1, &texture);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, _width, _height, 0, GL_RGBA, GL_FLOAT, NULL);
        glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

        m_shader.use();

        m_shader.execute(_width, _height, 1);

        vector<Vec3> image(_width * _height);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, image.data());
        glDeleteTextures(1, &texture);

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << " milliseconds." << endl;

        return image;
    }

    Renderer() : m_shader("ressources/shaders/compute.glsl") {}
};

#endif