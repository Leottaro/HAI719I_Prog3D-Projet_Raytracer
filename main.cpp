// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

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
#include <string>
#include <vector>

using namespace std;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX = 0, lastY = 0, lastZoom = 0;
static unsigned int FPS = 0;
static bool fullScreen = false;

vector<Renderer> renderers;

vector<pair<Vec3, Vec3>> rays;

void printUsage() {
    cerr << endl
         << "gMini: a minimal OpenGL/GLUT application" << endl
         << "for 3D graphics." << endl
         << "Author : Tamy Boubekeur (http://www.labri.fr/~boubek)" << endl
         << endl
         << "Usage : ./gmini [<file.off>]" << endl
         << "Keyboard commands" << endl
         << "------------------" << endl
         << " ?: Print help" << endl
         << " p: increase current preset" << endl
         << " P: decrease current preset" << endl
         << " w: Toggle Wireframe Mode" << endl
         << " g: Toggle Gouraud Shading Mode" << endl
         << " f: Toggle full screen mode" << endl
         << " <drag>+<left button>: rotate model" << endl
         << " <drag>+<right button>: move model" << endl
         << " <drag>+<middle button>: zoom" << endl
         << " q, <esc>: Quit" << endl
         << endl;
}

void usage() {
    printUsage();
    exit(EXIT_FAILURE);
}

// ------------------------------------
void initLight() {
    GLfloat light_position[4] = {0.0, 1.5, 0.0, 1.0};
    GLfloat color[4] = {1.0, 1.0, 1.0, 1.0};
    GLfloat ambient[4] = {1.0, 1.0, 1.0, 1.0};

    glLightfv(GL_LIGHT1, GL_POSITION, light_position);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT1, GL_SPECULAR, color);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
}

void initRenderers() {
    renderers.resize(5);
    renderers[0].setup_single_sphere();
    renderers[1].setup_single_square();
    renderers[2].setup_cornell_box();
    renderers[3].setup_single_mesh();
    renderers[4].setup_refraction_test();
}

void init() {
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cerr << "GLEW Initialization Error: " << glewGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }

    camera.resize(Settings::SCREEN_WIDTH, Settings::SCREEN_HEIGHT);
    initLight();
    initRenderers();
    // glCullFace (GL_BACK);
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
}

// ------------------------------------
// Replace the code of this
// functions for cleaning memory,
// closing sockets, etc.
// ------------------------------------

void clear() {
}

// ------------------------------------
// Replace the code of this
// functions for alternative rendering.
// ------------------------------------

void draw() {
    glEnable(GL_LIGHTING);
    renderers[Settings::selected_renderer].draw();

    // draw rays : (for debug)
    //  cout << rays.size() << endl;
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glLineWidth(6);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);
    for (unsigned int r = 0; r < rays.size(); ++r) {
        glVertex3f(rays[r].first[0], rays[r].first[1], rays[r].first[2]);
        glVertex3f(rays[r].second[0], rays[r].second[1], rays[r].second[2]);
    }
    glEnd();
}

void display() {
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply();
    draw();
    glFlush();
    glutSwapBuffers();
}

void idle() {
    static float lastTime = glutGet((GLenum)GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet((GLenum)GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char winTitle[64];
        sprintf(winTitle, "Raytracer - FPS: %d", FPS);
        glutSetWindowTitle(winTitle);
        lastTime = currentTime;
    }
    glutPostRedisplay();
}

void writePPM(vector<Vec3> image, string filename, int width, int height) {
    ofstream f(filename.c_str(), ios::binary);
    if (f.fail()) {
        cout << "Could not open file: " << filename << endl;
        return;
    }
    f << "P3" << endl
      << width << " " << height << endl
      << 255 << endl;
    for (int i = 0; i < width * height; i++)
        f << (int)(255.f * min<float>(1.f, image[i][0])) << " " << (int)(255.f * min<float>(1.f, image[i][1])) << " " << (int)(255.f * min<float>(1.f, image[i][2])) << " ";
    f << endl;
    f.close();
}

void key(unsigned char keyPressed, int x, int y) {
    Vec3 pos, dir;
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow(Settings::SCREEN_WIDTH, Settings::SCREEN_HEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen();
            fullScreen = true;
        }
        break;
    case 'q':
    case 27:
        clear();
        exit(0);
        break;
    case 'w':
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        if (polygonMode[0] != GL_FILL)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        break;

    case 'r':
        camera.apply();
        rays.clear();
        writePPM(
            renderers[Settings::selected_renderer].rayTraceFromCameraCPU(width, height, camera.getFarPlane()),
            "rendu.ppm",
            width,
            height);
        break;
    case 'R':
        camera.apply();
        rays.clear();
        writePPM(
            renderers[Settings::selected_renderer].rayTraceFromCameraGPU(width, height, camera.getFarPlane()),
            "renduGPU.ppm",
            width,
            height);
        break;
    case '+':
        Settings::selected_renderer = (Settings::selected_renderer + 1) % renderers.size();
        break;
    case 'p':
        Settings::selected_preset = static_cast<Settings::Presets>((static_cast<int>(Settings::selected_preset) + 1) % Settings::NB_PRESETS);
        Settings::applySelectedPreset();
        init();
        cout << "Settings preset set to: " << Settings::selected_preset << endl;
        break;
    case 'P':
        Settings::selected_preset = static_cast<Settings::Presets>((static_cast<int>(Settings::selected_preset) - 1 + Settings::NB_PRESETS) % Settings::NB_PRESETS);
        Settings::applySelectedPreset();
        init();
        cout << "Settings preset set to: " << Settings::selected_preset << endl;
        break;
    default:
        printUsage();
        break;
    }
    idle();
}

void mouse(int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate(x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle();
}

void motion(int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate(x, y);
    } else if (mouseMovePressed == true) {
        camera.move((x - lastX) / static_cast<float>(Settings::SCREEN_WIDTH), (lastY - y) / static_cast<float>(Settings::SCREEN_HEIGHT), 0.0);
        lastX = x;
        lastY = y;
    } else if (mouseZoomPressed == true) {
        camera.zoom(float(y - lastZoom) / Settings::SCREEN_HEIGHT);
        lastZoom = y;
    }
}

void reshape(int w, int h) {
    camera.resize(w, h);
}

int main(int argc, char **argv) {
    if (argc > 2) {
        printUsage();
        exit(EXIT_FAILURE);
    }

    Settings::selected_preset = Settings::Presets::PHASE_2_SOFT_SHADOWS;
    Settings::selected_renderer = 2;
    Settings::applySelectedPreset();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(Settings::SCREEN_WIDTH, Settings::SCREEN_HEIGHT);
    window = glutCreateWindow("gMini");

    init();
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    key('?', 0, 0);

    camera.move(0., 0., -3.1);
    renderers.resize(5);
    renderers[0].setup_single_sphere();
    renderers[1].setup_single_square();
    renderers[2].setup_cornell_box();
    renderers[3].setup_single_mesh();
    renderers[4].setup_refraction_test();

    glutMainLoop();
    return EXIT_SUCCESS;
}
