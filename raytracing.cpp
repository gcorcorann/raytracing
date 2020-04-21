#include <fstream>
#include <iostream>
#include <cmath>

#include "vector.h"
#include "sphere.h"
#include "ray.h"
//#include "triangle.h"
#include "image.h"
#include "light.h"
#include "material.h"

void pixelImagePlane(int i, int j, int nx, int ny, float l, float r, float t, float b, float& u, float& v) {
    u = l + (r - l) * ((float) i + 0.5) / (float) nx;
    v = b + (t - b) * ((float) j + 0.5) / (float) ny;
}

Ray perspectiveProjection(float u, float v, float f) {
    return {{0, 0, 0}, {u, v, -f}};
}

Ray orthographicProjection(float u, float v, float f, Ray& r) {
    return {{u, v, 0}, {0, 0, -1}};
}
//
///*
// * i pixel location (x)
// * j pixel location (y)
// * nx number of pixels (x)
// * ny number of pixels (y)
// * l left image plane location
// * r right image plane location
// * t top image plane location
// * b bottom image plane location
// * f focal length
// * ray viewing ray
// */
//void computeViewingRay(int i, int j, int nx, int ny, float l, float r, float t, float b, float f, Ray& ray) {
//    float u;
//    float v;
//    pixelImagePlane(i, j, nx, ny, l, r, t, b, u, v);
//    perspectiveProjection(u, v, f, ray);
//}

Ray computeViewingRayFOV(int i, int j, int nx, int ny, float fov, float f) {
    // pixel plane location
    float u; 
    float v;
    // image plane boundaries
    float r = tan(fov / 2);
    float l = -r;
    float t = r * ny / nx;
    float b = -t;
    pixelImagePlane(i, j, nx, ny, l, r, t, b, u, v);
    return perspectiveProjection(u, v, f);
}

/*
 * n surface normal
 * p surface hit position
 */
float lambertianShading(Vector n, Light light, Vector p) {
    Vector l = (light.position - p).unit();
    float r = (n.dot(l) > 0) ? n.dot(l) : 0;  // illumination proportion
    return r * light.intensity;
}

// TODO make r norm
float phongReflect(Ray ray, Vector n, Light light, Vector p, float specular_exponent) {
    Vector l = (light.position - p).unit();
    n = n.unit();
    Vector r = n.scale(n.dot(l) * 2) - l;
    Vector v = (ray.e - p).unit();
    float c = (v.dot(r) > 0) ? std::pow(v.dot(r), specular_exponent) : 0;  // illumination proportion
    return c * light.intensity;
}

Vector reflect(Vector d, Vector n) {
    d = d.unit();
    n = n.unit();
    Vector r = d - n.scale(d.dot(n) * 2);
    return r.unit();
}

/*
 * n surface normal
 * p surface hit position
 */
float blinnPhongShading(Vector n, Light light, Vector p, Ray ray, float specular_exponent) {
    Vector v = (ray.e - p).unit();
    Vector l = (light.position - p).unit();
    Vector h = (v + l).unit();
    float r = (n.dot(h) > 0) ? std::pow(n.dot(h), specular_exponent) : 0;  // illumination proportion
    return r * light.intensity;
}
//
///*
// * s surface ambient colour
// * i ambient light intensity
// */
//Vector ambientShading(Vector s, Vector i) {
//    return (s * i).scale(0.25);  // TODO scaling for now
//}
//
///*
// * colour hit surface colour
// */
//Vector mirrorReflection(Vector color){
//    return {0, 0, 0};
//}
//
///*
// * ns number of objects
// * n hit surface normal
// * p hit surface hit position
// * t hit distance
// * c hit surface colour
// */
//bool checkSurfaces(Surface* surfaces[], int ns, Ray ray, Vector& n, Vector& p, Vector& c, float& t) {
//    bool hit = false;
//    float min_t = -1;
//    Vector hit_norm;
//    Vector hit_position;
//    for (int k = 0; k < ns; k++) {
//        if (surfaces[k]->hit(ray, n, p, t) && (!hit || (t < min_t && t > 0))) {
//            hit = true;
//            min_t = t;
//            c = surfaces[k]->s;
//            hit_norm = n;
//            hit_position = p;
//        }
//    }
//    n = hit_norm;
//    p = hit_position;
//    t = min_t;
//    return hit;
//}
//
///*
// * p hit surface
// */
//bool checkShadow(Surface* surfaces[], int ns, Vector p, Light light) {
//    // shadow ray
//    Ray sray = {p, (light.p - p).unit()};
//    sray.o = sray.o + sray.d.scale(1);  // avoid ray hitting surface it's on
//    Vector hit_sn;
//    Vector hit_sp;
//    Vector hit_sc;
//    float hit_st;
//    return checkSurfaces(surfaces, ns, sray, hit_sn, hit_sp, hit_sc, hit_st);
//}
//
///*
// * colour surface colour of hit object
// */
//bool checkReflective(Surface* surfaces[], int ns, Vector p, Vector n, Ray ray, Vector& color) {
//    n = n.unit();
//    Vector d = ray.d.unit();
//    Vector r = d - n.scale(2 * (d.dot(n)));
//    Ray rr = {p, r};  // reflective ray
//    rr.o = rr.o + rr.d.scale(1);  // avoid ray hitting surface it's on
//    Vector hit_n;
//    Vector hit_p;
//    Vector colour;
//    float hit_t;
//    if (checkSurfaces(surfaces, ns, rr, hit_n, hit_p, colour, hit_t)) {
//        return true;
//    }
//    else {
//        return false;
//    }
//}

bool sceneIntersect(Ray ray, Sphere spheres [], int ns, Vector &n, Vector &p, Material &m) {
    bool hit = false;
    float min_t = -1;
    float t;
    Vector hit_n;
    Vector hit_p;
    for (int k = 0; k < ns; k++) {
        if (spheres[k].hit(ray, n, p, t) && (!hit || (t > 0 && t < min_t))) {
            hit = true;
            min_t = t;
            m = spheres[k].material;
            hit_n = n;
            hit_p = p;
        }
    }
    n = hit_n;
    p = hit_p;
    return hit;
}

//Vector castRay(Ray ray, Sphere spheres[], int ns, Light lights[], int nl, int depth = 0) {
Vector castRay(Ray ray, Sphere spheres[], int ns, Light lights[], int nl) {
    Vector n;  // surface normal
    Vector p;  // surface point
    Material m;
//    if (depth > 4 || !sceneIntersect(ray, spheres, ns, n, p, m)) {
    if (!sceneIntersect(ray, spheres, ns, n, p, m)) {
        return {0.2, 0.7, 0.8};  // background
    }
//    Ray rray;  // reflection ray
//    rray.d = reflect(ray.d, n);
//    rray.e = p + rray.d.scale(0.1);
//    Vector reflect_colour = castRay(rray, spheres, ns, lights, nl, depth+1);
//    Vector reflect_colour = {0, 0, 0};

    float diffuse_light_intensity = 0;
    float specular_light_intensity = 0;
    for (int k = 0; k < nl; k++) {
        Ray sray;
        sray.d = (lights[k].position - p).unit();
        sray.e = p + sray.d.scale(0.1);  // make sure point doesn't intersect itself
        Vector shadow_norm;
        Vector shadow_p;
        Material shadow_m;
        float light_dist = (lights[k].position - p).norm();  // check if shadow ray intersects object infront of light source
        if (sceneIntersect(sray, spheres, ns, shadow_norm, shadow_p, shadow_m) && (shadow_p - sray.e).norm() < light_dist) {
            continue;
        }
        diffuse_light_intensity += lambertianShading(n, lights[k], p);
        specular_light_intensity += phongReflect(ray, n, lights[k], p, m.specular_exponent);
    }
//    return m.diffuse_colour.scale(diffuse_light_intensity).scale(m.albedo.x) + Vector(1, 1, 1).scale(specular_light_intensity).scale(m.albedo.y) + reflect_colour.scale(m.albedo.z);
    return m.diffuse_colour.scale(diffuse_light_intensity).scale(m.albedo.x) + Vector(1, 1, 1).scale(specular_light_intensity).scale(m.albedo.y);
}

void render(Sphere spheres[], int ns, Light lights[], int nl, int width, int height, float fov, float focal_length) {
    Image img {width, height};
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            Ray ray = computeViewingRayFOV(i, j, width, height, fov, focal_length);
            *img[j * width + i] = castRay(ray, spheres, ns, lights, nl);
        }
    }
    img.normalizeImage();
    img.writeImage();
}

int main() {
    const int width = 1024;
    const int height = 768;
    float focal_length = 1;  // only for perspective projection
    float fov = M_PI / 2;
    Material ivory {{0.6, 0.3, 0.1}, {0.4, 0.4, 0.3}, 50.};
    Material red_rubber {{0.9, 0.1, 0.0}, {0.3, 0.1, 0.1}, 10.};
    Material mirror {{0.0, 10.0, 0.8}, {1.0, 1.0, 1.0}, 1425.};
    Sphere spheres [] {{{-3, 0, -16}, 2, ivory},  // centre, radius, material
                       {{-1, -1.5, -12}, 2, red_rubber},
                       {{1.5, -0.5, -18}, 3, ivory},
                       {{7, 5, -18}, 4, mirror}};
    Light lights [] {{{-20, 20, 20}, 1.5},
                     {{30, 50, -25}, 1.8},
                     {{30, 20, 30}, 1.7}};
    int ns = sizeof(spheres) / sizeof(*spheres);
    int nl = sizeof(lights) / sizeof(*lights);
    render(spheres, ns, lights, nl, width, height, fov, focal_length);
//    for (int j = 0; j < height; j++) {
//        for (int i = 0; i < width; i++) {
//            Ray ray = computeViewingRayFOV(i, j, width, height, fov, focal_length);
//            *img[j * width + 1] = castRay(
//
//            Vector hit_n;  // hit surface normal
//            Vector hit_p;  // hit surface point
//            Vector hit_c;  // hit surface colour
//            float hit_t;  //  hit intersection
//            bool hit = checkSurfaces(surfaces, ns, ray, hit_n, hit_p, hit_c, hit_t);
//            if (hit) {
//                Vector L = {0, 0, 0};
//                for (int k = 0; k < sizeof(lights) / sizeof(*lights); k++) {
//                    L = lambertianShading(hit_n, lights[k], hit_p, hit_c);
//                    if (!checkShadow(surfaces, ns, hit_p, lights[k])) {
//                        L = L + lambertianShading(hit_n, lights[k], hit_p, hit_c) + blinnPhongShading(hit_n, lights[k], hit_p, hit_c, ray);
//                        if (checkReflective(surfaces, ns, hit_p, hit_n, ray, hit_c)) {
//                            L = L + mirrorReflection(hit_c);
//                        }
//                    }
//                }
//                *img[j*width+i] = L + ambientShading(hit_c, lights[0].i);  // use first light as ambient
//                *img[j*width+i] = L;
//            }
//            else {
//                *img[j*width+i] = {0.2, 0.7, 0.8};  // background colour
//            }
//        }
//    }
//    img.normalizeImage();
//    img.writeImage();
    return 0;
}
