#pragma once
#include "surface.h"
#include "vector.h"
#include "ray.h"
#include "material.h"

class Sphere : public Surface {
public:
    Vector m_centre;
    float m_radius;
    Sphere(Vector c, float r, Material m) {
        m_centre = c;
        m_radius = r;
        m_material = m;
    }
    bool hit(Ray r, Vector& n, Vector& p, float& t) {
        float d = sqrt((r.e - m_centre).dot(r.d) * (r.e - m_centre).dot(r.d) - r.d.dot(r.d) * ((r.e - m_centre).dot(r.e - m_centre) - m_radius * m_radius));
        if (d >= 0) {
            float tp = (-(r.e - m_centre).dot(r.d) + d) / r.d.dot(r.d);
            float tn = (-(r.e - m_centre).dot(r.d) - d) / r.d.dot(r.d);
            t = (tp < tn) ? tp : tn;
            p = r.e + r.d.scale(t);
            n = (p - m_centre).unit();
            if (t >= 0) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
};
