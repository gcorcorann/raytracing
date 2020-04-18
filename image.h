#pragma once
#include <fstream>
#include "vector.h"

class Image {
    public:
        int width;
        int height;
        Vector* framebuffer = nullptr;
        Image(int w, int h) {
            width = w;
            height = h;
            framebuffer = new Vector [width * height];  // rgb values
        }
        ~Image() { delete[] framebuffer; };
        void writeImage() {
            std::ofstream file;
            file.open("out.ppm", std::ios::binary);
            file << "P6\n" << width << " " << height << "\n255\n";
            for (int j = 0; j < height; j++) {
                for (int i = 0; i < width; i++) {
                    file << (char)(255 * framebuffer[(height-1-j)*width+i].x);
                    file << (char)(255 * framebuffer[(height-1-j)*width+i].y);
                    file << (char)(255 * framebuffer[(height-1-j)*width+i].z);
                }
            }
            file.close();
        }
        void normalizeImage() {
            float max = 0;
            for (int j = 0; j < height; j++) {
                for (int i = 0; i < width; i++) {
                    float m = (framebuffer[j*width+i].x > framebuffer[j*width+i].y) ? framebuffer[j*width+i].x : framebuffer[j*width+i].y;
                    m = (m > framebuffer[j*width+i].z) ? m : framebuffer[j*width+i].z;
                    max = (max > m) ? max : m;
                }
            }
            if (max != 0) {  // TODO this might be changed to (max > 1)
                for (int j = 0; j < height; j++) {
                    for (int i = 0; i < width; i++) {
                        framebuffer[j*width+i] = {framebuffer[j*width+i].x / max, framebuffer[j*width+i].y / max, framebuffer[j*width+i].z / max};
                    }
                }
            }
        }
        Vector* operator [] (int i) { return &framebuffer[i]; };
};
