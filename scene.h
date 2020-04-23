#pragma once
#include "vector.h"
#include "sphere.h"
#include "image.h"
#include "light.h"
#include "ray.h"


class Scene {
private:
    int m_width;
    int m_height;
    float m_fov;
    float m_focal_length;
    int m_depth;
    Sphere* m_spheres;
    int m_ns;
    Light* m_lights;
    int m_nl;
public:
    Scene (const int width, const int height, float fov) {
        m_width = width;
        m_height = height;
        m_fov = fov;
        m_focal_length = 0.f;
        m_depth = 0;
    }
    Scene (const int width, const int height, float fov, float focal_length) {
        m_width = width;
        m_height = height;
        m_fov = fov;
        m_focal_length = focal_length;
        m_depth = 0;
    }
    Scene (const int width, const int height, float fov, float focal_length, int depth) {
        m_width = width;
        m_height = height;
        m_fov = fov;
        m_focal_length = focal_length;
        m_depth = depth;
    }

    void pixelImagePlane(int i, int j, float& u, float& v) {
        // image plane boundaries
        float r = tan(m_fov / 2);
        float l = -r;
        float t = r * m_height / m_width;
        float b = -t;
        // convert pixel to image plane location
        u = l + (r - l) * ((float) i + 0.5) / (float) m_width;
        v = b + (t - b) * ((float) j + 0.5) / (float) m_height;
    }
    Ray computeRay (int i, int j) {
        float u, v;  // pixel image plane location
        pixelImagePlane(i, j, u, v);
        if (m_focal_length == 0.f) {
            orthographicProjection(u, v);
        }
        return perspectiveProjection(u, v);
    }
    Ray perspectiveProjection(float u, float v) {
        return {{0, 0, 0}, {u, v, -m_focal_length}};  // origin, direction
    }
    Ray orthographicProjection(float u, float v) {
        return {{u, v, 0}, {0, 0, -1}};
    }
    void addSurfaces(Sphere* spheres, int ns) {
        m_spheres = spheres;
        m_ns = ns;
    }
    void addLights(Light* lights, int nl) {
        m_lights = lights;
        m_nl = nl;
    }
    bool intersect(Ray r, Vector& n, Vector& p, Material& m) {
        Vector tmp_n, tmp_p;
        float t;
        bool hit = false;
        float min_t = -1;
        for (int k = 0; k < m_ns; k++) {
            if (m_spheres[k].hit(r, tmp_n, tmp_p, t) && (!hit || (t > 0 && t < min_t))) {
                hit = true;
                min_t = t;
                m = m_spheres[k].material;
                n = tmp_n;
                p = tmp_p;
            }
        }
        return hit;
    }
    Vector castRay(Ray ray, int depth=0) {
        Material hit_mat;
        Vector hit_norm, hit_pt;
        if (depth > m_depth || !intersect(ray, hit_norm, hit_pt, hit_mat)) {
            return {0.2f, 0.7f, 0.8f};  // background colour
        }
        Ray reflect_ray {{hit_pt}, {reflect(ray.d, hit_norm)}};
        Vector reflect_colour = castRay(reflect_ray, ++depth).scale(hit_mat.albedo.z);
        Vector diffuse_colour, specular_colour;
        float diffuse_intensity = 0.f;
        float specular_intensity = 0.f;
        shading(ray, hit_norm, hit_pt, hit_mat.specular_exponent, diffuse_intensity, specular_intensity);
        diffuse_colour = hit_mat.diffuse_colour.scale(diffuse_intensity).scale(hit_mat.albedo.x);
        // white specular colour
        specular_colour = Vector(1.f, 1.f, 1.f).scale(specular_intensity).scale(hit_mat.albedo.y);
        return diffuse_colour + specular_colour + reflect_colour;
    }
    void shading(Ray r, Vector n, Vector p, float e, float& di, float& si) {
        Vector l;
        Ray shadow_ray;
        Vector s_hit_pt, s_hit_norm;
        Material s_hit_mat;
        for (int k = 0; k < m_nl; k++) {
            l = (m_lights[k].position - p).unit();
            shadow_ray = {p, l};
            if (intersect(shadow_ray, s_hit_norm, s_hit_pt, s_hit_mat)) {
                continue;
            }
            di += lambertianShading(l, m_lights[k].intensity, n, p);
            si += phongShading(r, l, m_lights[k].intensity, n, p, e);
        }
    }
    float lambertianShading(Vector l, float i, Vector n, Vector p) {
        float c = n.dot(l);
        c = (c > 0) ? c : 0;
        return c * i;
    }
    float phongShading(Ray ray, Vector l, float i, Vector n, Vector p, float e) {
        Vector r = n.scale(n.dot(l) * 2) - l;
        Vector v = (ray.e - p).unit();
        float c = v.dot(r);
        c = (c > 0) ? std::pow(c, e) : 0;
        return c * i;
    }
    float blinnPhongShading(Ray ray, Vector l, float i, Vector n, Vector p, float e) {
        Vector v = (ray.e - p).unit();
        Vector h = (v + l).unit();
        float c = n.dot(h);
        c = (c > 0) ? std::pow(c, e) : 0;
        return c * i;
    }
    Vector reflect(Vector d, Vector n) {
        return (d - n.scale(d.dot(n) * 2)).unit();
    }
    void render() {
        Image img {m_width, m_height};
        for (int j = 0; j < m_height; j++) {
            for (int i = 0; i < m_width; i++) {
                Ray ray = computeRay(i, j);
                *img[j * m_width + i] = castRay(ray);
            }
        }
        img.normalize();
        img.write();
    }
};
