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
private:
    vector<Vec3> image = vector<Vec3>();
    unsigned int image_w = 0, image_h = 0;

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

    void rayTraceFromCameraCPU(float _min_t, float _max_t, int64_t &elapsed_ms) {
        cout << "Ray tracing a " << Settings::SCREEN_WIDTH << " x " << Settings::SCREEN_HEIGHT << " image (CPU) :" << endl;
        Vec3 pos = cameraSpaceToWorldSpace(Vec3(0.));
        int n_pixels = Settings::SCREEN_WIDTH * Settings::SCREEN_HEIGHT;
        image = vector<Vec3>(n_pixels, Vec3(0, 0, 0));
        image_w = Settings::SCREEN_WIDTH;
        image_h = Settings::SCREEN_HEIGHT;

        auto begin = chrono::high_resolution_clock::now();

        for (unsigned int y = 0; y < image_h; y++) {
            for (unsigned int x = 0; x < image_w; x++) {
                int pixel_i = y * image_w + x + 1;

                float percent = (float)pixel_i / n_pixels;
                int64_t currently_elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - begin).count();
                float remaining_ms = (currently_elapsed / percent) * (1. - percent);

                cout << "\r\tCalculating pixel " << pixel_i << " of " << n_pixels << " (" << fixed << setprecision(2) << 100.f * percent << "% completed) ~" << remaining_ms / 1000. << "s remaining     " << flush;
                for (unsigned int s = 0; s < Settings::NSAMPLES; ++s) {
                    float u = ((float)(x) + (float)(rand()) / (float)(RAND_MAX)) / image_w;
                    float v = ((float)(y) + (float)(rand()) / (float)(RAND_MAX)) / image_h;
                    Vec3 dir = screen_space_to_worldSpace(u, v) - pos;
                    Vec3 color = rayTrace(Ray(pos, dir), _min_t, _max_t);
                    image[x + y * image_w] += color;
                }
                image[x + y * image_w] /= (float)Settings::NSAMPLES;
            }
        }
        auto end = chrono::high_resolution_clock::now();
        elapsed_ms = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << endl
             << "\tDone in " << elapsed_ms << "ms." << endl;
    }

    void writeImage(string filename) {
        ofstream f(filename.c_str(), ios::binary);
        if (f.fail()) {
            cout << "Could not open file: " << filename << endl;
            return;
        }
        f << "P3" << endl
          << image_w << " " << image_h << endl
          << 255 << endl;
        for (unsigned int i = 0; i < image_w * image_h; i++)
            f << (int)(255.f * min<float>(1.f, image[i][0])) << " " << (int)(255.f * min<float>(1.f, image[i][1])) << " " << (int)(255.f * min<float>(1.f, image[i][2])) << " ";
        f << endl;
        f.close();
    }

    // =====================================================================================================
    // ================================================ GPU ================================================
    // =====================================================================================================

private:
    GLuint m_out_texture;
    GLuint m_spheres_ssbo, m_squares_ssbo, m_lights_ssbo, m_mesh_vertices_ssbo, m_mesh_triangles_ssbo, m_meshes_ssbo;

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

    void updateGeneralUniforms(ComputeShader &shader, float _min_t, float _max_t) {
        shader.set("imgOutput", 0);

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

        shader.set("nearAndFarPlanes", nearAndFarPlanes[0], nearAndFarPlanes[1]);
        shader.set("projectionInverse", projectionInverse);
        shader.set("modelviewInverse", modelviewInverse);

        shader.set("min_t", _min_t);
        shader.set("max_t", _max_t);
    }

    void updateSettingsUniforms(ComputeShader &shader) {
        shader.set("settings.EPSILON", Settings::EPSILON);
        shader.set("settings.NSAMPLES", Settings::NSAMPLES);
        shader.set("settings.MAX_BOUNCES", Settings::MAX_BOUNCES);
        shader.set("settings.SCREEN_WIDTH", Settings::SCREEN_WIDTH);
        shader.set("settings.SCREEN_HEIGHT", Settings::SCREEN_HEIGHT);

        shader.set("settings.Phong.ENABLED", Settings::Phong::ENABLED);
        shader.set("settings.Phong.SHADOW_RAYS", Settings::Phong::SHADOW_RAYS);

        shader.set("settings.Material.ENABLE_MIRROR", Settings::Material::ENABLE_MIRROR);
        shader.set("settings.Material.ENABLE_GLASS", Settings::Material::ENABLE_GLASS);
        shader.set("settings.Material.AIR_INDEX_MEDIUM", Settings::Material::AIR_INDEX_MEDIUM);

        shader.set("settings.Mesh.ENABLE_INTERPOLATION", Settings::Mesh::ENABLE_INTERPOLATION);

        shader.set("settings.KdTree.EPSILON", Settings::KdTree::EPSILON);
        shader.set("settings.KdTree.MAX_LEAF_SIZE", Settings::KdTree::MAX_LEAF_SIZE);

        shader.set("settings.Bonus.ENABLE_TEXTURES", Settings::Bonus::ENABLE_TEXTURES);
    }

    void updateUniforms(ComputeShader &shader, float _min_t, float _max_t) {
        cout << "\tUpdating uniforms..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        shader.use();
        updateGeneralUniforms(shader, _min_t, _max_t);
        updateSettingsUniforms(shader);

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void executeShader(ComputeShader &shader) {
        cout << "\tExecuting shader..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        shader.execute(Settings::SCREEN_WIDTH, Settings::SCREEN_HEIGHT, Settings::NSAMPLES);

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

    void retrieveImage() {
        cout << "\tRetrieving the image..." << flush;
        auto begin = chrono::high_resolution_clock::now();

        image.resize(Settings::SCREEN_WIDTH * Settings::SCREEN_HEIGHT);
        image_w = Settings::SCREEN_WIDTH;
        image_h = Settings::SCREEN_HEIGHT;
        glBindTexture(GL_TEXTURE_2D, m_out_texture);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, image.data());

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tDone in " << elapsed << "ms" << endl;
    }

public:
    void rayTraceFromCameraGPU(ComputeShader &shader, float _min_t, float _max_t, int64_t &elapsed_ms) {
        auto begin = chrono::high_resolution_clock::now();
        cout << "Ray tracing a " << Settings::SCREEN_WIDTH << " x " << Settings::SCREEN_HEIGHT << " image (GPU) :" << endl;

        createTextures();
        createSceneBuffers();
        updateUniforms(shader, _min_t, _max_t);
        executeShader(shader);
        retrieveImage();
        deleteSceneBuffers();
        deleteTextures();

        glFinish();
        auto end = chrono::high_resolution_clock::now();
        elapsed_ms = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
        cout << "\tTotal time: " << elapsed_ms << "ms" << endl;
    }

    Renderer() {}
};

#endif