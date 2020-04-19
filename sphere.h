#pragma once
#include "vector.h"
#include "ray.h"
#include "surface.h"

class Sphere : public Surface {
    public:
        Vector c;  // centre
        float r;   // radius
        Sphere(Vector centre, float radius, Vector sc) {
            c = centre;
            r = radius;
            s = sc;
        }
        /*
         * r viewing ray
         * n surface normal
         * p surface hit location
         * t distance along ray of hit location
         */
        bool hit(Ray r, Vector& n, Vector& p, float& t) {
            float d = sqrt((r.o - c).dot(r.d) * (r.o - c).dot(r.d) - r.d.dot(r.d) * ((r.o - c).dot(r.o - c) - this->r * this->r));
            if (d >= 0) {
                float tp = (-(r.o - c).dot(r.d) + d) / r.d.dot(r.d);
                float tn = (-(r.o - c).dot(r.d) - d) / r.d.dot(r.d);
                t = (tp < tn) ? tp : tn;
                p = r.o + r.d.scale(t);
                n = (p - c).unit();
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
