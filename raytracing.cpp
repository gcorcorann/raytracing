#include <fstream>
#include <iostream>
#include <cmath>

#include "vector.h"
#include "sphere.h"
#include "ray.h"

struct Light {
    Vector p;  // position
    Vector i;  // intensity
};

void pixelImagePlane(int i, int j, int nx, int ny, float l, float r, float t, float b, float& u, float& v) {
    u = l + (r - l) * ((float) i + 0.5) / (float) nx;
    v = b + (t - b) * ((float) j + 0.5) / (float) ny;
}

void perspectiveProjection(float u, float v, float f, Ray& r) {
    r.o = {0, 0, 0};
    r.d = {u, v, -f};
}

void orthographicProjection(float u, float v, float f, Ray& r) {
    r.o = {u, v, 0};
    r.d = {0, 0, -1};
}

void computeViewingRay(int i, int j, int nx, int ny, float l, float r, float t, float b, float f, Ray& ray) {
    float u;
    float v;
    pixelImagePlane(i, j, nx, ny, l, r, t, b, u, v);
    orthographicProjection(u, v, f, ray);
}

Vector lambertianShading(Vector n, Light light, Vector p, Vector s) {
    Vector l = light.p - p;
    l = l.scale(1 / l.norm());
    float r = (n.dot(l) > 0) ? n.dot(l) : 0;  // illumination proportion
    Vector dc = {s.x * light.i.x * r, s.y * light.i.y * r, s.z * light.i.z * r};  // diffuse component
    return (s * light.i).scale(r);
}

/*
 * n surface normal
 * p surface hit position
 * s surface colour
 */
Vector blinnPhongShading(Vector n, Light light, Vector p, Vector s, Ray ray) {
    Vector v = (ray.o - p).unit();
    Vector l = (light.p - p).unit();
    Vector h = (v + l).unit();
    float r = (n.dot(h) > 0) ? std::pow(n.dot(h), 25) : 0;  // illumination proportion
    return (s * light.i).scale(r);
}

/*
 * s surface ambient colour
 * i ambient light intensity
 */
Vector ambientShading(Vector s, Vector i) {
    return (s * i).scale(0.25);  // TODO scaling for now
}

void writeImage(Vector* img, int width, int height) {
    std::ofstream file;
    file.open("out.ppm", std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            file << (char)(255 * img[(height-1-j)*width+i].x);
            file << (char)(255 * img[(height-1-j)*width+i].y);
            file << (char)(255 * img[(height-1-j)*width+i].z);
        }
    }
    file.close();
}

void normalizeImage(Vector* img, int width, int height) {
    float max = 0;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            float m = (img[j*width+i].x > img[j*width+i].y) ? img[j*width+i].x : img[j*width+i].y;
            m = (m > img[j*width+i].z) ? m : img[j*width+i].z;
            max = (max > m) ? max : m;
        }
    }
    if (max != 0) {  // TODO this might be changed to (max > 1)
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                img[j*width+i] = {img[j*width+i].x / max, img[j*width+i].y / max, img[j*width+i].z / max};
            }
        }
    }
}

int main() {
    const int width = 1280;
    const int height = 720;
    float focal_length = 200;  // only for perspective projection
    // image plane location
    float r = 200;
    float t = 150;
    float l = -r;
    float b = -t;
    Vector* img = new Vector [width * height];  // rgb image
    Ray ray;  // viewing ray
    // centre, radius, colour
    Sphere spheres [] = {{30, 30, -300, 30, 1, 0, 0},
                         {50, 50, -400, 20, 1, 1, 1},
                         {-40, -40, -300, 40, 0, 0, 1}};
    // loc, colour
    Light lights [] = {{320, 300, 0, 0.5, 0.5, 0.5},
                       {-1020, -310, 0, 0.5, 0.5, 0.5}};
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            computeViewingRay(i, j, width, height, l, r, t, b, focal_length, ray);
            Vector n;  // surface normal
            Vector p;  // surface hit point
            float t;  // intersection
            float min_t = -1;
            Sphere hit_sphere;
            Vector hit_position;
            Vector hit_norm;
            bool hit = false;
            for (int k = 0; k < sizeof(spheres) / sizeof(*spheres); k++) {
                if (spheres[k].hit(ray, n, p, t) && (!hit || t < min_t)) {
                    hit = true;
                    min_t = t;
                    hit_sphere = spheres[k];
                    hit_position = p;
                    hit_norm = n;
                }
            }
            if (hit) {
                // diffuse + specular
                Vector L = {0, 0, 0};
                for (int k = 0; k < sizeof(lights) / sizeof(*lights); k++) {
                    L = L + lambertianShading(hit_norm, lights[k], hit_position, hit_sphere.s) + blinnPhongShading(hit_norm, lights[k], hit_position, hit_sphere.s, ray);
                }
                // use first light as ambient
                img[j*width+i] = L + ambientShading(hit_sphere.s, lights[0].i);
            }
            else {
                float bg = (float)(height - j) / (height * 3);
                img[j*width+i] = {bg, bg, bg};
            }
        }
    }
    normalizeImage(img, width, height);
    writeImage(img, width, height);
    delete[] img;
    return 0;
}
