#include <fstream>
#include <iostream>
#include <cmath>

#include "vector.h"
#include "sphere.h"
#include "triangle.h"
#include "ray.h"
#include "image.h"
#include "light.h"

void pixelImagePlane(int i, int j, int nx, int ny, float l, float r, float t, float b, float& u, float& v) {
    u = l + (r - l) * ((float) i + 0.5) / (float) nx;
    v = b + (t - b) * ((float) j + 0.5) / (float) ny;
}

void perspectiveProjection(float u, float v, float f, Ray& r) {
    r.o = {0, 30, 0};
    r.d = {u, v, -f};
}

void orthographicProjection(float u, float v, float f, Ray& r) {
    r.o = {u, v, 0};
    r.d = {0, 0, -1};
}

/*
 * i pixel location (x)
 * j pixel location (y)
 * nx number of pixels (x)
 * ny number of pixels (y)
 * l left image plane location
 * r right image plane location
 * t top image plane location
 * b bottom image plane location
 * f focal length
 * ray viewing ray
 */
void computeViewingRay(int i, int j, int nx, int ny, float l, float r, float t, float b, float f, Ray& ray) {
    float u;
    float v;
    pixelImagePlane(i, j, nx, ny, l, r, t, b, u, v);
    perspectiveProjection(u, v, f, ray);
}

/*
 * n surface normal
 * p surface hit position
 * s surface colour
 */
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

/*
 * ns number of objects
 * n hit surface normal
 * p hit surface hit position
 * t hit distance
 * c hit surface colour
 */
bool checkSurfaces(Surface* surfaces[], int ns, Ray ray, Vector& n, Vector& p, Vector& c, float& t) {
    bool hit = false;
    float min_t = -1;
    Vector hit_norm;
    Vector hit_position;
    for (int k = 0; k < ns; k++) {
        if (surfaces[k]->hit(ray, n, p, t) && (!hit || (t < min_t && t > 0))) {
            hit = true;
            min_t = t;
            c = surfaces[k]->s;
            hit_norm = n;
            hit_position = p;
        }
    }
    n = hit_norm;
    p = hit_position;
    t = min_t;
    return hit;
}

int main() {
    const int width = 1280;
    const int height = 720;
    float focal_length = 200;  // only for perspective projection
    float r = 200 / 1.5;  // image place locations
    float t = 150 / 1.5;
    float l = -r;
    float b = -t;
    Image img = {width, height};  // rgb image
    Ray ray;  // viewing ray
    Sphere spheres [] = {{{30, 30, -300}, 30, {1, 0, 0}},  // centre, radius, colour
                         {{80, 20, -400}, 20, {1, 1, 1}},
                         {{-40, 40, -300}, 40, {0, 0, 1}}};
    Triangle triangle = {{-20000, 0, -20000}, {20000, 0, -20000}, {0, 0, 20000}, {0.5, 0.5, 0.5}};  // v1, v2, v3, colour
    Surface* surfaces [] = {&triangle, &spheres[0], &spheres[1], &spheres[2]};
    int ns = sizeof(surfaces) / sizeof(*surfaces);
    Light lights [] = {{320, 300, 0, 0.5, 0.5, 0.5}, {-320, 300, 0, 0.5, 0.5, 0.5}, {0, 100, 0, 0.5, 0.5, 0.5}};  // loc, colour
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            computeViewingRay(i, j, width, height, l, r, t, b, focal_length, ray);
            Vector hit_n;  // hit surface normal
            Vector hit_p;  // hit surface point
            Vector hit_c;  // hit surface colour
            float hit_t;  //  hit intersection
            bool hit = checkSurfaces(surfaces, ns, ray, hit_n, hit_p, hit_c, hit_t);
            if (hit) {
                Vector L = {0, 0, 0};
                for (int k = 0; k < sizeof(lights) / sizeof(*lights); k++) {
                    // shadow ray
                    Ray sray = {hit_p, (lights[k].p - hit_p).unit()};
                    sray.o = sray.o + sray.d.scale(1);  // avoid ray hitting surface it's on
                    Vector hit_sn;
                    Vector hit_sp;
                    Vector hit_sc;
                    float hit_st;
                    if (!checkSurfaces(surfaces, ns, sray, hit_sn, hit_sp, hit_sc, hit_st)) {
                        L = L + lambertianShading(hit_n, lights[k], hit_p, hit_c) + blinnPhongShading(hit_n, lights[k], hit_p, hit_c, ray);
                    }
                }
                *img[j*width+i] = L + ambientShading(hit_c, lights[0].i);  // use first light as ambient
            }
            else {
                *img[j*width+i] = {0, 0, 0};
            }
        }
    }
    img.normalizeImage();
    img.writeImage();
    return 0;
}
