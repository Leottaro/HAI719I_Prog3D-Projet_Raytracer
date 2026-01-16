#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "src/Camera.h"
#include "src/Renderer.h"
#include "src/Scene.h"
#include "src/Settings.h"
#include "src/Vec3.h"
#include "src/matrixUtilities.h"
#include <GL/glut.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

class Benchmark {
private:
    static bool openBenchmarkCSV(string filename, string headers, ofstream &writer) {
        if (system("mkdir -p benchmarks") != 0) {
            cerr << "Error: Unable to create directory benchmarks." << endl;
            return false;
        }

        string filepath = "benchmarks/" + filename + ".csv";

        ifstream reader(filepath);
        bool is_empty = !reader.is_open() || reader.peek() == EOF;
        reader.close();

        writer = ofstream(filepath, ios::app);
        if (!writer.is_open()) {
            cerr << "Error: Unable to open file \"" << filepath << "\" for writing." << endl;
            return false;
        }
        if (is_empty) {
            writer << headers << endl;
        }

        return true;
    }

public:
    static void meshBenchmark(vector<unsigned int> const &max_leaf_sizes, Renderer &renderer) {
        ofstream writer;
        if (!openBenchmarkCSV("raw_mesh_benchmark", "Mesh,MaxLeafSize,ElapsedMs", writer))
            return;

        int64_t elapsed_ms;
        for (Settings::AvailableMeshes mesh : Settings::AVAILABLE_MESHES) {
            for (unsigned int max_leaf_size : max_leaf_sizes) {

                Settings::applyPreset(Settings::Presets::BENCHMARKS);
                // Settings::NSAMPLES = 1;
                // Settings::MAX_BOUNCES = 100;
                Settings::Mesh::MESH = mesh;
                // Settings::Mesh::ENABLE_INTERPOLATION = false;
                // Settings::Phong::ENABLED = false;
                // Settings::Phong::SHADOW_RAYS = 0;
                // Settings::Material::ENABLE_MIRROR = false;
                // Settings::Material::ENABLE_GLASS = false;
                // Settings::Material::AIR_INDEX_MEDIUM = 1.0;
                Settings::KdTree::MAX_LEAF_SIZE = max_leaf_size;
                // Settings::Bonus::ENABLE_TEXTURES = false;<< renderer.meshes[0].getNbVertices() << ";" << renderer.meshes[0].getNbTriangles() <<

                cout << mesh << " (" << max_leaf_size << ")" << endl;
                renderer.setup_mesh_benchmark();
                renderer.rayTraceFromCameraCPU(Settings::EPSILON, FLT_MAX, elapsed_ms);
                renderer.writeImage("last_benchmark.ppm");
                writer << Settings::Mesh::MESH << "," << Settings::KdTree::MAX_LEAF_SIZE << "," << elapsed_ms << endl;

                if (static_cast<unsigned int>(mesh) <= max_leaf_size) {
                    break;
                }
            }
        }
    }

    static void gpuBenchmark(vector<unsigned int> const &nbs_objects, ComputeShader &shader, Renderer &renderer) {
        ofstream writer;
        if (!openBenchmarkCSV("raw_gpu_benchmark", "NbObjects,CPUElapsedMs,GPUElapsedMs", writer))
            return;

        int64_t cpu_elapsed_ms, gpu_elapsed_ms;
        for (unsigned int nb_objects : nbs_objects) {
            Settings::applyPreset(Settings::Presets::BENCHMARKS);
            Settings::NSAMPLES = 1;
            Settings::MAX_BOUNCES = 100;
            // Settings::Mesh::MESH = Settings::AvailableMeshes::TRIANGLE; // No meshes in this test
            Settings::Mesh::ENABLE_INTERPOLATION = true;
            Settings::Phong::ENABLED = true;
            Settings::Phong::SHADOW_RAYS = 16;
            Settings::Material::ENABLE_MIRROR = true;
            Settings::Material::ENABLE_GLASS = true;
            Settings::Material::AIR_INDEX_MEDIUM = 1.0;
            // Settings::KdTree::MAX_LEAF_SIZE = 0; // No kdTree GPU side
            // Settings::Bonus::ENABLE_TEXTURES = false; // No textures GPU side

            cout << "benchmarking " << nb_objects << " objects..." << endl;
            renderer.setup_gpu_benchmark(nb_objects);
            renderer.rayTraceFromCameraCPU(Settings::EPSILON, FLT_MAX, cpu_elapsed_ms);
            renderer.writeImage("last_benchmark.ppm");
            renderer.rayTraceFromCameraGPU(shader, Settings::EPSILON, FLT_MAX, gpu_elapsed_ms);
            renderer.writeImage("last_benchmark.ppm");
            writer << nb_objects << "," << cpu_elapsed_ms << "," << gpu_elapsed_ms << endl;
        }
    }
};

#endif