// raytracing.cpp - simple raytracing
#include <iostream>
#include <fstream>
#include <cmath>


struct Vector {
    float x;
    float y;
    float z;
};

struct Sphere {
    Vector c;  // sphere centre
    float r;  // sphere radius

    bool ray_intersect(Vector p, Vector d, float& di1) {
        Vector v = {c.x - p.x, c.y - p.y, c.z - p.z};
        float vdotd = d.x * v.x + d.y * v.y + d.z + v.z;
        // sphere is infront of ray
        if (vdotd > 0) {
            float m = vdotd / std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);  // |d|
            Vector pc = {v.x * m + v.y * m + v.z * m};  // project of c on the line
            // no intersection
            if (std::sqrt((v.x - pc.x) * (v.x - pc.x) + (v.y - pc.y) * (v.y - pc.y) + (v.z - pc.z) * (v.z - pc.z)) > r) {
                return false;
            }
            else {
                // pythagorean theorem (distance from pc to il)
                float dist = std::sqrt(r * r - ((v.x - pc.x) * (v.x - pc.x) + (v.y - pc.y) * (v.y - pc.y) + (v.z - pc.z) * (v.z - pc.z)));
                di1 = (pc.x - p.x) * (pc.x - p.x) + (pc.y - p.y) * (pc.y - p.y) + (pc.z - p.z) * (pc.z - p.z) - std::sqrt(dist);
                return true;
            }
        }
        else {
            return false;
        }
    }
};

/*
 * render an image
 */
void render(Sphere sphere) {
    const int width = 1024;
    const int height = 768;

    float* r = new float[width * height];
    float* g = new float[width * height];
    float* b = new float[width * height];

    // gradient image colour to black
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            r[i * width + j] = 0;
            g[i * width + j] = 0;
            b[i * width + j] = 0;
        }
    }

//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            float x = (2 * (j + 0.5) / (float)width - 1) * std::tan(1.5708/2) * width / (float) height;
//            float y = -(2 * (i + 0.5) / (float) height - 1) * std::tan(1.5708 / 2);
//            Vector d = {x, y, -1};
//            Vector p = {0, 0, 0};
//            float di1 = 0;
//            if (sphere.ray_intersect(p, d, di1)) {
//                std::cout << "Here" << std::endl;
//                r[i * width + j] = 1;
//                g[i * width + j] = 1;
//                b[i * width + j] = 1;
//            }
//        }
//    }

    // save image to file
    std::ofstream file;
    file.open("out.ppm");
    file << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height; i++) {
        file << (char)(255 * 1);
        file << (char)(255 * 1);
        file << (char)(255 * 1);
    }
    file.close();

    // clear memory
    delete[] r;
    delete[] g;
    delete[] b;
}


int main() {
    Vector c = {0, 0, 2};  // sphere centre
    Sphere sphere = {c, 2};  // radius 10

    Vector p = {0, 0, 0};
    Vector d = {0, 0, 1};
    render(sphere);
}

// TODO make sure we are adjusting everythign by d
