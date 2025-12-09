#ifndef MESH_H
#define MESH_H

#include "Constants.h"
#include "KdTree.h"
#include "Material.h"
#include "Ray.h"
#include "Triangle.h"
#include "Vec3.h"
#include <GL/glut.h>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using namespace std;

// -------------------------------------------
// Basic Mesh class
// -------------------------------------------

struct MeshVertex {
    inline MeshVertex() {}
    inline MeshVertex(const Vec3 &_p, const Vec3 &_n) : position(_p), normal(_n), u(0), v(0) {}
    inline MeshVertex(const MeshVertex &vertex) : position(vertex.position), normal(vertex.normal), u(vertex.u), v(vertex.v) {}
    inline virtual ~MeshVertex() {}
    inline MeshVertex &operator=(const MeshVertex &vertex) {
        position = vertex.position;
        normal = vertex.normal;
        u = vertex.u;
        v = vertex.v;
        return (*this);
    }
    // membres :
    Vec3 position; // une position
    Vec3 normal;   // une normale
    float u, v;    // coordonnees uv
};

struct MeshTriangle {
    inline MeshTriangle() {
        v[0] = v[1] = v[2] = 0;
    }
    inline MeshTriangle(const MeshTriangle &t) {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
    }
    inline MeshTriangle(unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
    }
    unsigned int &operator[](unsigned int iv) { return v[iv]; }
    unsigned int operator[](unsigned int iv) const { return v[iv]; }
    inline virtual ~MeshTriangle() {}
    inline MeshTriangle &operator=(const MeshTriangle &t) {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
        return (*this);
    }
    // membres :
    unsigned int v[3];
};

class Mesh {
protected:
    void build_positions_array() {
        positions_array.resize(3 * vertices.size());
        for (unsigned int v = 0; v < vertices.size(); ++v) {
            positions_array[3 * v + 0] = vertices[v].position[0];
            positions_array[3 * v + 1] = vertices[v].position[1];
            positions_array[3 * v + 2] = vertices[v].position[2];
        }
    }
    void build_normals_array() {
        normalsArray.resize(3 * vertices.size());
        for (unsigned int v = 0; v < vertices.size(); ++v) {
            normalsArray[3 * v + 0] = vertices[v].normal[0];
            normalsArray[3 * v + 1] = vertices[v].normal[1];
            normalsArray[3 * v + 2] = vertices[v].normal[2];
        }
    }
    void build_UVs_array() {
        uvs_array.resize(2 * vertices.size());
        for (unsigned int vert = 0; vert < vertices.size(); ++vert) {
            uvs_array[2 * vert + 0] = vertices[vert].u;
            uvs_array[2 * vert + 1] = vertices[vert].v;
        }
    }
    void build_triangles_array() {
        triangles_array.resize(3 * triangles.size());
        for (unsigned int t = 0; t < triangles.size(); ++t) {
            triangles_array[3 * t + 0] = triangles[t].v[0];
            triangles_array[3 * t + 1] = triangles[t].v[1];
            triangles_array[3 * t + 2] = triangles[t].v[2];
        }
    }
    void build_kd_tree() {
        vector<KdTriangle> passed_triangles(this->triangles.size());
        Vec3 min = Vec3(FLT_MAX);
        Vec3 max = Vec3(-FLT_MAX);

        for (size_t i = 0; i < this->triangles.size(); i++) {
            KdTriangle kd_triangle = KdTriangle(this->vertices[this->triangles[i][0]].position,
                                                this->vertices[this->triangles[i][1]].position,
                                                this->vertices[this->triangles[i][2]].position,
                                                i);
            min = Vec3(
                std::min({min[0], kd_triangle.v0[0], kd_triangle.v1[0], kd_triangle.v2[0]}),
                std::min({min[1], kd_triangle.v0[1], kd_triangle.v1[1], kd_triangle.v2[1]}),
                std::min({min[2], kd_triangle.v0[2], kd_triangle.v1[2], kd_triangle.v2[2]}));
            max = Vec3(
                std::max({max[0], kd_triangle.v0[0], kd_triangle.v1[0], kd_triangle.v2[0]}),
                std::max({max[1], kd_triangle.v0[1], kd_triangle.v1[1], kd_triangle.v2[1]}),
                std::max({max[2], kd_triangle.v0[2], kd_triangle.v1[2], kd_triangle.v2[2]}));
            passed_triangles[i] = kd_triangle;
        }

        this->kdtree = make_unique<KdTree>(passed_triangles, BoundingBox(min, max), (max - min).getMaxAbsoluteComponent());
    }

public:
    vector<MeshVertex> vertices;
    vector<MeshTriangle> triangles;

    vector<float> positions_array;
    vector<float> normalsArray;
    vector<float> uvs_array;
    vector<unsigned int> triangles_array;

    Material material;
    std::unique_ptr<KdTree> kdtree = nullptr;

    void loadOFF(const string &filename);
    void recomputeNormals();
    void centerAndScaleToUnit();
    void scaleUnit();

    virtual void build_arrays() {
        recomputeNormals();
        build_positions_array();
        build_normals_array();
        build_UVs_array();
        build_triangles_array();
        if (constants::kdtree::MAX_LEAF_SIZE > 0) {
            build_kd_tree();
        }
    }

    void translate(Vec3 const &translation) {
        for (unsigned int v = 0; v < vertices.size(); ++v) {
            vertices[v].position += translation;
        }
    }

    void apply_transformation_matrix(Mat3 transform) {
        for (unsigned int v = 0; v < vertices.size(); ++v) {
            vertices[v].position = transform * vertices[v].position;
        }

        //        recomputeNormals();
        //        build_positions_array();
        //        build_normals_array();
    }

    void scale(Vec3 const &scale) {
        Mat3 scale_matrix(scale[0], 0., 0.,
                          0., scale[1], 0.,
                          0., 0., scale[2]); // Matrice de transformation de mise à l'échelle
        apply_transformation_matrix(scale_matrix);
    }

    void rotate_x(float angle) {
        float x_angle = angle * M_PI / 180.;
        Mat3 x_rotation(1., 0., 0.,
                        0., cos(x_angle), -sin(x_angle),
                        0., sin(x_angle), cos(x_angle));
        apply_transformation_matrix(x_rotation);
    }

    void rotate_y(float angle) {
        float y_angle = angle * M_PI / 180.;
        Mat3 y_rotation(cos(y_angle), 0., sin(y_angle),
                        0., 1., 0.,
                        -sin(y_angle), 0., cos(y_angle));
        apply_transformation_matrix(y_rotation);
    }

    void rotate_z(float angle) {
        float z_angle = angle * M_PI / 180.;
        Mat3 z_rotation(cos(z_angle), -sin(z_angle), 0.,
                        sin(z_angle), cos(z_angle), 0.,
                        0., 0., 1.);
        apply_transformation_matrix(z_rotation);
    }

    void draw() const {
        if (triangles_array.size() == 0)
            return;
        GLfloat material_color[4] = {material.diffuse_material[0],
                                     material.diffuse_material[1],
                                     material.diffuse_material[2],
                                     1.0};

        GLfloat material_specular[4] = {material.specular_material[0],
                                        material.specular_material[1],
                                        material.specular_material[2],
                                        1.0};

        GLfloat material_ambient[4] = {material.ambient_material[0],
                                       material.ambient_material[1],
                                       material.ambient_material[2],
                                       1.0};

        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 3 * sizeof(float), (GLvoid *)(normalsArray.data()));
        glVertexPointer(3, GL_FLOAT, 3 * sizeof(float), (GLvoid *)(positions_array.data()));
        glDrawElements(GL_TRIANGLES, triangles_array.size(), GL_UNSIGNED_INT, (GLvoid *)(triangles_array.data()));
    }

    RayTriangleIntersection intersect(Ray const &ray) const {
        return constants::kdtree::MAX_LEAF_SIZE > 0 ? intersect_kdtree(ray) : intersect_no_kdtree(ray);
    }

    RayTriangleIntersection intersect_no_kdtree(Ray const &ray) const {
        RayTriangleIntersection closestIntersection = RayTriangleIntersection();

        // Creer un objet Triangle pour chaque face
        // Vous constaterez des problemes de précision
        // solution : ajouter un facteur d'échelle lors de la création du Triangle :
        // float triangleScaling = 1.000001;
        size_t n = triangles.size();

        for (size_t i = 0; i < n; i++) {
            MeshTriangle mesh_triangle = triangles[i];
            MeshVertex v0 = vertices[mesh_triangle[0]];
            MeshVertex v1 = vertices[mesh_triangle[1]];
            MeshVertex v2 = vertices[mesh_triangle[2]];
            Triangle triangle = Triangle(v0.position, v1.position, v2.position);
            RayTriangleIntersection intersection = triangle.getIntersection(ray);
            if (intersection.intersectionExists && intersection.t < closestIntersection.t) {
                intersection.tIndex = i;
                intersection.normal = v0.normal * intersection.w0 + v1.normal * intersection.w1 + v2.normal * intersection.w2;
                intersection.normal.normalize();
                intersection.u = v0.u * intersection.w0 + v1.u * intersection.w1 + v2.u * intersection.w2;
                intersection.v = v0.v * intersection.w0 + v1.v * intersection.w1 + v2.v * intersection.w2;
                closestIntersection = intersection;
            }
        }
        return closestIntersection;
    }

    RayTriangleIntersection intersect_kdtree(Ray const &ray) const {
        KdTriangle kd_triangle;
        bool intersect = kdtree->intersect(ray, kd_triangle);
        if (!intersect) {
            return RayTriangleIntersection();
        }

        Triangle triangle = Triangle(kd_triangle.v0, kd_triangle.v1, kd_triangle.v2);
        MeshTriangle mesh_triangle = this->triangles[kd_triangle.triangle_index];
        MeshVertex mesh_v0 = this->vertices[mesh_triangle[0]];
        MeshVertex mesh_v1 = this->vertices[mesh_triangle[1]];
        MeshVertex mesh_v2 = this->vertices[mesh_triangle[2]];

        RayTriangleIntersection intersection = triangle.getIntersection(ray);
        intersection.tIndex = kd_triangle.triangle_index;
        intersection.normal = mesh_v0.normal * intersection.w0 + mesh_v1.normal * intersection.w1 + mesh_v2.normal * intersection.w2;
        intersection.normal.normalize();
        intersection.u = mesh_v0.u * intersection.w0 + mesh_v1.u * intersection.w1 + mesh_v2.u * intersection.w2;
        intersection.v = mesh_v0.v * intersection.w0 + mesh_v1.v * intersection.w1 + mesh_v2.v * intersection.w2;
        return intersection;
    }
};

#endif
