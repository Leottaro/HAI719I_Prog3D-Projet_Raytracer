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

    vector<Vec3> rayTraceFromCameraCPU(float _min_t, float _max_t) {
        cout << "Ray tracing a " << Settings::SCREEN_WIDTH << " x " << Settings::SCREEN_HEIGHT << " image (CPU) :" << endl;
        Vec3 pos = cameraSpaceToWorldSpace(Vec3(0.));
        int n_pixels = Settings::SCREEN_WIDTH * Settings::SCREEN_HEIGHT;
        vector<Vec3> image(n_pixels, Vec3(0, 0, 0));

        auto begin = chrono::high_resolution_clock::now();

        for (unsigned int y = 0; y < Settings::SCREEN_HEIGHT; y++) {
            for (unsigned int x = 0; x < Settings::SCREEN_WIDTH; x++) {
                int pixel_i = y * Settings::SCREEN_WIDTH + x + 1;

                float percent = (float)pixel_i / n_pixels;
                int64_t currently_elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count();
                float remaining_ms = (currently_elapsed / percent) * (1. - percent);

                cout << "\r\tCalculating pixel " << pixel_i << " of " << n_pixels << " (" << fixed << setprecision(2) << 100.f * percent << "% completed) ~" << remaining_ms / 1000. << "s remaining     " << flush;
                for (unsigned int s = 0; s < Settings::NSAMPLES; ++s) {
                    float u = ((float)(x) + (float)(rand()) / (float)(RAND_MAX)) / Settings::SCREEN_WIDTH;
                    float v = ((float)(y) + (float)(rand()) / (float)(RAND_MAX)) / Settings::SCREEN_HEIGHT;
                    Vec3 dir = screen_space_to_worldSpace(u, v) - pos;
                    Vec3 color = rayTrace(Ray(pos, dir), _min_t, _max_t);
                    image[x + y * Settings::SCREEN_WIDTH] += color;
                }
                image[x + y * Settings::SCREEN_WIDTH] /= (float)Settings::NSAMPLES;
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
    GLuint m_out_texture;
    GLuint m_spheres_ssbo, m_squares_ssbo, m_lights_ssbo, m_mesh_vertices_ssbo, m_mesh_triangles_ssbo, m_meshes_ssbo, m_images_array;

    void createTextures() {
        // Create the output texture
        cout << "\tCreating texture..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        glGenTextures(1, &m_out_texture);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_out_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, Settings::SCREEN_WIDTH, Settings::SCREEN_HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
        glBindImageTexture(0, m_out_texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);

        glGenTextures(1, &m_images_array);
        glActiveTexture(GL_TEXTURE7);
        glBindTexture(GL_TEXTURE_2D_ARRAY, m_images_array);
        unsigned int nb_images = images.size();
        for (size_t i = 0; i < nb_images; i++) {
            vector<float> image_data(images[i].w * images[i].h * 4);
            for (size_t j = 0; j < images[i].w * images[i].h; j++) {
                image_data[j * 4] = float(images[i].data[j].r) / 255.;
                image_data[j * 4 + 1] = float(images[i].data[j].g) / 255.;
                image_data[j * 4 + 2] = float(images[i].data[j].b) / 255.;
                image_data[j * 4 + 3] = 1.;
            }
            glTexStorage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGBA32F, images[i].w, images[i].h, 1);
            glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, i, images[i].w, images[i].h, 1, GL_RGBA, GL_FLOAT, image_data.data());
        }

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void deleteTextures() {
        cout << "\tDeleting textures..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        glDeleteTextures(1, &m_out_texture);

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void createSceneBuffers() {
        cout << "\tCreating buffer..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        unsigned int nb_spheres = spheres.size();
        vector<ShaderSphere> shader_spheres(nb_spheres);
        for (size_t i = 0; i < nb_spheres; i++) {
            shader_spheres[i].assign(spheres[i]);
        }
        glGenBuffers(1, &m_spheres_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_spheres_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 16 + nb_spheres * sizeof(ShaderSphere), NULL, GL_DYNAMIC_DRAW);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(unsigned int), &nb_spheres);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 16, nb_spheres * sizeof(ShaderSphere), shader_spheres.data());
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_spheres_ssbo);

        unsigned int nb_squares = squares.size();
        vector<ShaderSquare> shader_squares(nb_squares);
        for (size_t i = 0; i < nb_squares; i++) {
            shader_squares[i].assign(squares[i]);
        }
        glGenBuffers(1, &m_squares_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_squares_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 16 + nb_squares * sizeof(ShaderSquare), NULL, GL_DYNAMIC_DRAW);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(unsigned int), &nb_squares);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 16, nb_squares * sizeof(ShaderSquare), shader_squares.data());
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_squares_ssbo);

        unsigned int nb_lights = lights.size();
        vector<ShaderLight> shader_lights(nb_lights);
        for (size_t i = 0; i < nb_lights; i++) {
            shader_lights[i].assign(lights[i]);
        }
        glGenBuffers(1, &m_lights_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_lights_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 16 + nb_lights * sizeof(ShaderLight), NULL, GL_DYNAMIC_DRAW);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(unsigned int), &nb_lights);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 16, nb_lights * sizeof(ShaderLight), shader_lights.data());
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, m_lights_ssbo);

        unsigned int nb_meshes = meshes.size();
        vector<ShaderMeshVertex> shader_mesh_vertices;
        vector<ShaderMeshTriangle> shader_mesh_triangles;
        vector<ShaderMesh> shader_meshes(nb_meshes);
        for (size_t i = 0; i < nb_meshes; i++) {
            vector<MeshVertex> const &vertices = meshes[i].vertices;
            vector<MeshTriangle> const &triangles = meshes[i].triangles;
            Material const &material = meshes[i].material;

            shader_meshes[i].nb_vertices = vertices.size();
            shader_meshes[i].vertices_offset = shader_mesh_vertices.size();
            shader_meshes[i].nb_triangles = triangles.size();
            shader_meshes[i].triangles_offset = shader_mesh_triangles.size();
            shader_meshes[i].material.assign(material);

            vector<ShaderMeshVertex> translated_vertices(vertices.size());
            for (size_t j = 0; j < vertices.size(); j++) {
                translated_vertices[j].assign(vertices[j]);
            }
            shader_mesh_vertices.reserve(shader_mesh_vertices.size() + distance(translated_vertices.begin(), translated_vertices.end()));
            shader_mesh_vertices.insert(shader_mesh_vertices.end(), translated_vertices.begin(), translated_vertices.end());

            vector<ShaderMeshTriangle> translated_triangles(triangles.size());
            for (size_t j = 0; j < triangles.size(); j++) {
                translated_triangles[j].assign(triangles[j]);
            }
            shader_mesh_triangles.reserve(shader_mesh_triangles.size() + distance(translated_triangles.begin(), translated_triangles.end()));
            shader_mesh_triangles.insert(shader_mesh_triangles.end(), translated_triangles.begin(), translated_triangles.end());
        }

        glGenBuffers(1, &m_mesh_vertices_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_mesh_vertices_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, shader_mesh_vertices.size() * sizeof(ShaderMeshVertex), shader_mesh_vertices.data(), GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, m_mesh_vertices_ssbo);

        glGenBuffers(1, &m_mesh_triangles_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_mesh_triangles_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, shader_mesh_triangles.size() * sizeof(ShaderMeshTriangle), shader_mesh_triangles.data(), GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, m_mesh_triangles_ssbo);

        glGenBuffers(1, &m_meshes_ssbo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_meshes_ssbo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, 16 + nb_meshes * sizeof(ShaderMesh), NULL, GL_DYNAMIC_DRAW);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(unsigned int), &nb_meshes);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 16, nb_meshes * sizeof(ShaderMesh), shader_meshes.data());
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, m_meshes_ssbo);

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void deleteSceneBuffers() {
        cout << "\tDeleting buffers..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        glDeleteBuffers(1, &m_spheres_ssbo);
        glDeleteBuffers(1, &m_squares_ssbo);
        glDeleteBuffers(1, &m_lights_ssbo);
        glDeleteBuffers(1, &m_mesh_vertices_ssbo);
        glDeleteBuffers(1, &m_mesh_triangles_ssbo);
        glDeleteBuffers(1, &m_meshes_ssbo);

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void updateGeneralUniforms(float _min_t, float _max_t) {
        GLdouble projection[16];
        GLdouble projectionInverse[16];
        glMatrixMode(GL_PROJECTION);
        glGetDoublev(GL_PROJECTION_MATRIX, projection);
        gluInvertMatrix(projection, projectionInverse);
        GLdouble nearAndFarPlanes[2];
        glGetDoublev(GL_DEPTH_RANGE, nearAndFarPlanes);
        GLdouble modelview[16];
        GLdouble modelviewInverse[16];
        glMatrixMode(GL_MODELVIEW);
        glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
        gluInvertMatrix(modelview, modelviewInverse);

        // m_shader.set("imgOutput", 0);
        m_shader.set("nearAndFarPlanes", nearAndFarPlanes[0], nearAndFarPlanes[1]);
        m_shader.set("projectionInverse", projectionInverse);
        m_shader.set("modelviewInverse", modelviewInverse);
        m_shader.set("min_t", _min_t);
        m_shader.set("max_t", _max_t);
        m_shader.set("images", 7);
        // m_shader.set("IMAGE", 8);

        // for (size_t i = 0; i < images.size(); i++) {
        //     m_shader.set("images[" + std::to_string(i) + "]", m_images[i]);
        // }
    }

    void updateSettingsUniforms() {
        m_shader.set("settings.EPSILON", Settings::EPSILON);
        m_shader.set("settings.NSAMPLES", Settings::NSAMPLES);
        m_shader.set("settings.MAX_BOUNCES", Settings::MAX_BOUNCES);
        m_shader.set("settings.SCREEN_WIDTH", Settings::SCREEN_WIDTH);
        m_shader.set("settings.SCREEN_HEIGHT", Settings::SCREEN_HEIGHT);

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

    void updateUniforms(float _min_t, float _max_t) {
        cout << "\tUpdating uniforms..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        m_shader.use();
        updateGeneralUniforms(_min_t, _max_t);
        updateSettingsUniforms();

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void executeShader() {
        cout << "\tExecuting shader..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        m_shader.execute(Settings::SCREEN_WIDTH, Settings::SCREEN_HEIGHT, Settings::NSAMPLES);

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    vector<Vec3> retrieveImage() {
        cout << "\tRetrieving the image..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        vector<Vec3> image(Settings::SCREEN_WIDTH * Settings::SCREEN_HEIGHT);
        glBindTexture(GL_TEXTURE_2D, m_out_texture);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, image.data());

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
        return image;
    }

public:
    vector<Vec3> rayTraceFromCameraGPU(float _min_t, float _max_t) {
        auto begin = chrono::high_resolution_clock::now();
        cout << "Ray tracing a " << Settings::SCREEN_WIDTH << " x " << Settings::SCREEN_HEIGHT << " image (GPU) :" << endl;

        createTextures();
        createSceneBuffers();
        updateUniforms(_min_t, _max_t);
        executeShader();
        vector<Vec3> image = retrieveImage();
        deleteSceneBuffers();
        deleteTextures();

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tTotal time: " << elapsed << "ms" << endl;

        return image;
    }

    Renderer() : m_shader("ressources/shaders/compute.glsl") {}
};

#endif