#pragma once
#include "vector.h"
#include "ray.h"

class Triangle {
    public:
        Vector v1;  // vertex position
        Vector v2;
        Vector v3;
        Vector s;  // surface colour

        /*
         * r viewing ray
         * n surface normal
         * p surface hit location
         * t distance along ray of hit location
         */
        bool hit(Ray r, Vector& n, Vector& p, float& t) {
            float a = v1.x - v2.x;
            float b = v1.y - v2.y;
            float c = v1.z - v2.z;
            float d = v1.x - v3.x;
            float e = v1.y - v3.y;
            float f = v1.z - v3.z;
            float g = r.d.x;
            float h = r.d.y;
            float i = r.d.z;
            float j = v1.x - r.o.x;
            float k = v1.y - r.o.y;
            float l = v1.z - r.o.z;
            float M = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g);
            float beta = j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g);
            beta /= M;
            float gamma = i * (a * k - j * b) + h * (j * c - a * l) + g * (b * l - k * c);
            gamma /= M;
            t = -(f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c));
            t /= M;
            if (t > 0 && (0 <= gamma <= 1) && (0 <= beta <= 1 - gamma)) {
                p = r.o + r.d.scale(t);
                n = {0, 1, 0};
                return true;
            }
            else {
                return false;
            }
        }
};
