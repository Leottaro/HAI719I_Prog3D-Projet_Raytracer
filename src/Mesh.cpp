#include "Mesh.h"
#include <fstream>
#include <iostream>

using namespace std;

void Mesh::loadOFF(const string &filename) {
    ifstream in(filename.c_str());
    if (!in)
        exit(EXIT_FAILURE);
    string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    // cout << "loading mesh at \"" << filename << "\"" << endl
    //      << "siseV: " << sizeV << endl
    //      << "sizeT: " << sizeT << endl;
    vertices.resize(sizeV);
    triangles.resize(sizeT);
    for (unsigned int i = 0; i < sizeV; i++) {
        in >> vertices[i].position;
        // cout << "position: " << vertices[i].position;
        if (!(in.peek() == '\n' || in.peek() == '\r' || in.eof())) {
            in >> vertices[i].normal;
            // cout << ", normal: " << vertices[i].normal;
        }
        // cout << endl;
    }
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        // cout << "Triangle: ";
        for (unsigned int j = 0; j < 3; j++) {
            in >> triangles[i].v[j];
            // cout << triangles[i].v[j] << ", ";
        }
        if (!(in.peek() == '\n' || in.peek() == '\r' || in.eof())) {
            string restOfLine;
            getline(in, restOfLine);
            // cout << "and some things";
        }
        // cout << endl;
    }
    in.close();
}

void Mesh::recomputeNormals() {
    for (unsigned int i = 0; i < vertices.size(); i++)
        vertices[i].normal = Vec3();
    for (unsigned int i = 0; i < triangles.size(); i++) {
        Vec3 e01 = vertices[triangles[i].v[1]].position - vertices[triangles[i].v[0]].position;
        Vec3 e02 = vertices[triangles[i].v[2]].position - vertices[triangles[i].v[0]].position;
        Vec3 n = Vec3::cross(e01, e02);
        n.normalize();
        for (unsigned int j = 0; j < 3; j++)
            vertices[triangles[i].v[j]].normal += n;
    }
    for (unsigned int i = 0; i < vertices.size(); i++)
        vertices[i].normal.normalize();
}

void Mesh::centerAndScaleToUnit() {
    Vec3 c(0, 0, 0);
    for (unsigned int i = 0; i < vertices.size(); i++)
        c += vertices[i].position;
    c /= vertices.size();
    float maxD = (vertices[0].position - c).length();
    for (unsigned int i = 0; i < vertices.size(); i++) {
        float m = (vertices[i].position - c).length();
        if (m > maxD)
            maxD = m;
    }
    for (unsigned int i = 0; i < vertices.size(); i++)
        vertices[i].position = (vertices[i].position - c) / maxD;
}
