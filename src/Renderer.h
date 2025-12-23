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
    unsigned int m_width, m_height;
    Vec3 m_camera_pos;

    GLuint m_ray_texture, m_random_texture, m_out_texture;

    void createTextures() {
        // rays data
        vector<float> rayTextureData(m_width * Settings::NSAMPLES * m_height * 4);
        Vec3 dir;
        for (unsigned int y = 0; y < m_height; y++) {
            for (unsigned int x = 0; x < m_width; x++) {
                for (unsigned int z = 0; z < Settings::NSAMPLES; ++z) {
                    unsigned int i = y * m_width * Settings::NSAMPLES + x * Settings::NSAMPLES + z;
                    float u = ((float)(x) + (float)(rand()) / (float)(RAND_MAX)) / m_width;
                    float v = ((float)(y) + (float)(rand()) / (float)(RAND_MAX)) / m_height;
                    screen_space_to_world_space_ray(u, v, m_camera_pos, dir); // TODO: fix: ça casse la scène
                    rayTextureData[i * 4] = dir[0];
                    rayTextureData[i * 4 + 1] = dir[1];
                    rayTextureData[i * 4 + 2] = dir[2];
                    rayTextureData[i * 4 + 3] = 0.;
                }
            }
        }

        // randoms data
        vector<float> randomTextureData(m_width * Settings::NSAMPLES * m_height * (Settings::Phong::SHADOW_RAYS + 1) * 4);
        for (unsigned int i = 0; i < m_width * Settings::NSAMPLES * m_height * (Settings::Phong::SHADOW_RAYS + 1); i++) {
            randomTextureData[i * 4] = float(rand()) / RAND_MAX;
            randomTextureData[i * 4 + 1] = float(rand()) / RAND_MAX;
            randomTextureData[i * 4 + 2] = float(rand()) / RAND_MAX;
            randomTextureData[i * 4 + 3] = 0.;
        }

        // Create the output texture
        glGenTextures(1, &m_out_texture);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_out_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, m_width, m_height, 0, GL_RGBA, GL_FLOAT, NULL);
        glBindImageTexture(0, m_out_texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);

        // Create the rays directions
        glGenTextures(1, &m_ray_texture);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, m_ray_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, m_width * Settings::NSAMPLES, m_height, 0, GL_RGBA, GL_FLOAT, rayTextureData.data());
        glBindImageTexture(1, m_ray_texture, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);

        // Create random texture
        glGenTextures(1, &m_random_texture);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, m_random_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, m_width * Settings::NSAMPLES, m_height * (Settings::Phong::SHADOW_RAYS + 1), 0, GL_RGBA, GL_FLOAT, randomTextureData.data());
        glBindImageTexture(2, m_random_texture, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);
    }

    void updateSettings() {
        m_shader.set("settings.EPSILON", Settings::EPSILON);
        m_shader.set("settings.NSAMPLES", Settings::NSAMPLES);
        m_shader.set("settings.MAX_BOUNCES", Settings::MAX_BOUNCES);

        m_shader.set("settings.Phong.ENABLED", Settings::Phong::ENABLED);
        m_shader.set("settings.Phong.SHADOW_RAYS", Settings::Phong::SHADOW_RAYS);

        m_shader.set("settings.Material.ENABLE_MIRROR", Settings::Material::ENABLE_MIRROR);
        m_shader.set("settings.Material.ENABLE_GLASS", Settings::Material::ENABLE_GLASS);
        m_shader.set("settings.Material.AIR_INDEX_MEDIUM", Settings::Material::AIR_INDEX_MEDIUM);

        m_shader.set("settings.Mesh.ENABLE_INTERPOLATION", Settings::Mesh::ENABLE_INTERPOLATION);

        m_shader.set("settings.KdTree.EPSILON", Settings::KdTree::EPSILON);
        m_shader.set("settings.KdTree.MAX_LEAF_SIZE", Settings::KdTree::MAX_LEAF_SIZE);

        m_shader.set("settings.Bonus.ENABLE_TEXTURES", Settings::Bonus::ENABLE_TEXTURES);
    }

    void updateSceneUniforms() {
        // objects uniforms
        GLuint nb_spheres = spheres.size();
        m_shader.set("nb_spheres", nb_spheres);
        for (unsigned int i = 0; i < nb_spheres; i++) {
            m_shader.set("spheres[" + std::to_string(i) + "]", spheres[i]);
        }
        GLuint nb_squares = squares.size();
        m_shader.set("nb_squares", nb_squares);
        for (unsigned int i = 0; i < nb_squares; i++) {
            m_shader.set("squares[" + std::to_string(i) + "]", squares[i]);
        }
        // GLuint nb_meshes = meshes.size();
        // m_shader.set("nb_meshes", nb_meshes);
        // for (unsigned int i = 0; i < nb_meshes; i++) {
        //     m_shader.set("meshes[" + std::to_string(i) + "]", meshes[i]);
        // }
        GLuint nb_lights = lights.size();
        m_shader.set("nb_lights", nb_lights);
        for (unsigned int i = 0; i < nb_lights; i++) {
            m_shader.set("lights[" + std::to_string(i) + "]", lights[i]);
        }
    }

public:
    vector<Vec3> rayTraceFromCameraGPU(unsigned int _width, unsigned int _height, int _max_t) {
        auto begin = chrono::high_resolution_clock::now();
        cout << "Ray tracing a " << m_width << " x " << m_height << " image (GPU) :" << endl;

        m_width = _width;
        m_height = _height;

        cout << "\tGenerating textures..." << endl;
        createTextures();

        // use and execute the shader
        cout << "\tUpdating uniforms..." << endl;
        m_shader.use();
        m_shader.set("imgOutput", 0);
        m_shader.set("rayTexture", 1);
        m_shader.set("randomTexture", 2);

        m_shader.set("screen_width", m_width);
        m_shader.set("screen_height", m_height);
        m_shader.set("camera_pos", m_camera_pos);

        updateSettings();
        updateSceneUniforms();

        cout << "\tExecuting shader..." << endl;
        m_shader.execute(m_width, m_height, Settings::NSAMPLES);

        // retrieve the image
        cout << "\tRetrieving the image..." << endl;
        vector<Vec3> image(m_width * m_height);
        glBindTexture(GL_TEXTURE_2D, m_out_texture);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, image.data());
        glDeleteTextures(1, &m_out_texture);
        glDeleteTextures(1, &m_ray_texture);
        glDeleteTextures(1, &m_random_texture);

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << " milliseconds." << endl;

        return image;
    }

    Renderer() : m_shader("ressources/shaders/compute.glsl") {}
};

#endif