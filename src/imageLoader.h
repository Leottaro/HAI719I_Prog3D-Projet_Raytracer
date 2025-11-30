#ifndef IMAGELOADER_H
#define IMAGELOADER_H

#include "Vec3.h"
#include <fstream>
#include <string>
#include <vector>

// Source courtesy of J. Manson
// http://josiahmanson.com/prose/optimize_ppm/

namespace ppmLoader {
using namespace std;
void eat_comment(ifstream &f);

struct RGB {
    unsigned char r, g, b;
};

struct ImageRGB {
public:
    size_t w, h;
    vector<RGB> data;

    Vec3 getPixel(size_t x, size_t y) const;
    Vec3 getPixel(float u, float v) const;
};

void load_ppm(ImageRGB &img, const string &name);

enum loadedFormat {
    rgb,
    rbg
};

void load_ppm(unsigned char *&pixels, unsigned int &w, unsigned int &h, const string &name, loadedFormat format = rgb);
} // namespace ppmLoader

#endif
