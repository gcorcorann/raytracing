#pragma once
#include "vector.h"
#include "ray.h"
#include "material.h"

class Sphere {
    public:
        Vector centre;  // centre
        float radius;   // radius
        Material material;  // object material
        Sphere(Vector c, float r, Material m) : material(m) {
            centre = c;
            radius = r;
        }
        /*
         * r viewing ray
         * n surface normal
         * p surface hit location
         * t distance along ray of hit location
         */
        bool hit(Ray r, Vector& n, Vector& p, float& t) {
            float d = sqrt((r.e - centre).dot(r.d) * (r.e - centre).dot(r.d) - r.d.dot(r.d) * ((r.e - centre).dot(r.e - centre) - radius * radius));
            if (d >= 0) {
                float tp = (-(r.e - centre).dot(r.d) + d) / r.d.dot(r.d);
                float tn = (-(r.e - centre).dot(r.d) - d) / r.d.dot(r.d);
                t = (tp < tn) ? tp : tn;
                p = r.e + r.d.scale(t);
                n = (p - centre).unit();
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
