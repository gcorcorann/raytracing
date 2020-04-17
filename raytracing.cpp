#include <fstream>
#include <iostream>
#include <cmath>

#include "vector.h"

struct Sphere {
    Vector c;  // centre
    float r;  // radius
    Vector s;  // surface colour
};

struct Ray {
    Vector o;  // origin
    Vector d;  // diriection
};

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

bool sphereIntersection(Sphere s, Ray r, Vector& n, Vector& p, float& t) {
    float d = sqrt(r.o.sub(s.c).dot(r.d) * r.o.sub(s.c).dot(r.d) - r.d.dot(r.d) * (r.o.sub(s.c).dot(r.o.sub(s.c)) - s.r * s.r));
    if (d >= 0) {
        float tp = (-r.o.sub(s.c).dot(r.d) + d) / r.d.dot(r.d);
        float tn = (-r.o.sub(s.c).dot(r.d) - d) / r.d.dot(r.d);
        t = (tp < tn) ? tp : tn;
        p = r.o.add(r.d.scale(t));
        n = p.sub(s.c).unit();
        return true;
    }
    else 
        return false;
}

Vector lambertianShading(Vector n, Light light, Vector p, Vector s) {
    Vector l = light.p.sub(p);
    l = l.scale(1 / l.norm());
    float r = (n.dot(l) > 0) ? n.dot(l) : 0;  // illumination proportion
    Vector dc = {s.x * light.i.x * r, s.y * light.i.y * r, s.z * light.i.z * r};  // diffuse component
    return s.mul(light.i).scale(r);
}

/*
 * n surface normal
 * p surface hit position
 * s surface colour
 */
Vector blinnPhongShading(Vector n, Light light, Vector p, Vector s, Ray ray) {
    Vector v = ray.o.sub(p).unit();
    Vector l = light.p.sub(p).unit();
    Vector h = v.add(l).unit();
    float r = (n.dot(h) > 0) ? std::pow(n.dot(h), 25) : 0;  // illumination proportion
    return s.mul(light.i).scale(r);
}

/*
 * s surface ambient colour
 * i ambient light intensity
 */
Vector ambientShading(Vector s, Vector i) {
    return s.mul(i).scale(0.25);  // TODO scaling for now
}

int main() {
    const int width = 800;
    const int height = 640;
    float focal_length = 200;  // only for perspective projection
    // image plane location
    float r = 100;
    float l = -r;
    float t = r;
    float b = l;
    Vector img [height][width];  // RGB image
    Ray ray;  // viewing ray
    // centre, radius, colour
    Sphere spheres [] = {{30, 30, -300, 30, 0.4, 0.4, 0.4},
                         {50, 50, -400, 20, 0.3, 0.3, 0.3},
                         {-40, -40, -300, 40, 0.25, 0.25, 0.25}};
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
                if (sphereIntersection(spheres[k], ray, n, p, t) && (!hit || t < min_t)) {
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
                    L = L.add((lambertianShading(hit_norm, lights[k], hit_position, hit_sphere.s).add(blinnPhongShading(hit_norm, lights[k], hit_position, hit_sphere.s, ray))));
                }
                // use first light as ambient
                img[j][i] = L.add(ambientShading(hit_sphere.s, lights[0].i));
            }
            else {
                float bg = (float)(height - j) / (height * 3);
                img[j][i] = {bg, bg, bg};
            }
        }
    }
    // normalize image
    float max = 0;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            float m = (img[height-1-j][i].x > img[height-1-j][i].y) ? img[height-1-j][i].x : img[height-1-j][i].y;
            m = (m > img[height-1-j][i].z) ? m : img[height-1-j][i].z;
            if (m > max) {
                max = m;
            }
        }
    }
    if (max != 0) {
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                img[height-1-j][i] = {img[height-1-j][i].x / max, img[height-1-j][i].y / max, img[height-1-j][i].z / max};
            }
        }
    }
    std::ofstream file;
    file.open("out.ppm", std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            file << (char)(255 * img[height - 1 - j][i].x);
            file << (char)(255 * img[height - 1 - j][i].y);
            file << (char)(255 * img[height - 1 - j][i].z);
        }
    }
    file.close();
    return 0;
}
