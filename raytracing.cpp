#include <fstream>
#include <iostream>
#include <cmath>

struct Vector {
    float x;
    float y;
    float z;
};

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

float dot(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector sub(Vector a, Vector b) {
    Vector r = {a.x - b.x, a.y - b.y, a.z - b.z};
    return r;
}

Vector add(Vector a, Vector b) {
    Vector r = {a.x + b.x, a.y + b.y, a.z + b.z};
    return r;
}

Vector scale(Vector a, float t) {
    Vector r = {a.x * t, a.y * t, a.z * t};
    return r;
}

float norm(Vector a) {
    return sqrt(dot(a, a));
}

Vector unit(Vector a) {
    return scale(a, 1 / norm(a));
}

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
    float d = sqrt(dot(r.d, sub(r.o, s.c)) * dot(r.d, sub(r.o, s.c)) - dot(r.d, r.d) * (dot(sub(r.o, s.c), sub(r.o, s.c)) - s.r * s.r));
    if (d >= 0) {
        float tp = (-dot(r.d, sub(r.o, s.c)) + d) / dot(r.d, r.d);
        float tn = (-dot(r.d, sub(r.o, s.c)) - d) / dot(r.d, r.d);
        t = (tp < tn) ? tp : tn;
        p = add(r.o, scale(r.d, t));
        n = unit(sub(p, s.c));
        return true;
    }
    else 
        return false;
}

Vector lambertianShading(Vector n, Light light, Vector p, Vector s) {
    Vector l = sub(light.p, p);
    l = scale(l, 1 / norm(l));
    float r = (dot(n, l) > 0) ? dot(n, l) : 0;  // illumination proportion
    Vector d = {s.x * light.i.x * r, s.y * light.i.y * r, s.z * light.i.z * r};
    return d;
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
    Sphere spheres [] = {{30, 30, -300, 30, 0.5, 0.5, 0.5},
                         {50, 50, -400, 30, 0.5, 0, 0},
                         {-40, -40, -300, 40, 0, 1, 0}};
    

    Light light = {0, 10, 0, 1, 1, 1};  // loc, colour
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
            if (hit)
                img[j][i] = lambertianShading(hit_norm, light, hit_position, hit_sphere.s);  // diffuse colour
            else
                img[j][i] = {0, 0, 0};
        }
    }
    std::ofstream file;
    file.open("out.ppm");
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
