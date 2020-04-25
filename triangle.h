#pragma once
#include "surface.h"
#include "vector.h"
#include "ray.h"

class Triangle : public Surface {
public:
    Vector m_a;
    Vector m_b;
    Vector m_c;
    Triangle(Vector a, Vector b, Vector c, Material m) {
        m_a = a;
        m_b = b;
        m_c = c;
        m_material = m;
    }
    bool hit(Ray r, Vector& n, Vector& p, float& t) override {
        float a = m_a.x - m_b.x;
        float b = m_a.y - m_b.y;
        float c = m_a.z - m_b.z;
        float d = m_a.x - m_c.x;
        float e = m_a.y - m_c.y;
        float f = m_a.z - m_c.z;
        float g = r.d.x;
        float h = r.d.y;
        float i = r.d.z;
        float j = m_a.x - r.e.x;
        float k = m_a.y - r.e.y;
        float l = m_a.z - r.e.z;
        float M = a*(e*i - h*f) + b*(g*f - d*i) + c*(d*h - e*g);
        float beta = (j*(e*i - h*f) + k*(g*f - d*i) + l*(d*h - e*g)) / M;
        float gamma = (i*(a*k - j*b) + h*(j*c - a*l) + g*(b*l - k*c)) / M;
        t = -(f*(a*k - j*b) + e*(j*c - a*l) + d*(b*l - k*c)) / M;
        if (t >= 0 && gamma >= 0 && gamma <= 1 && beta >= 0 && beta <= 1 - gamma) {
            p = r.e + r.d.scale(t);
            n = {0, 1, 0};
            return true;
        }
        else {
            return false;
        }
    }
};
