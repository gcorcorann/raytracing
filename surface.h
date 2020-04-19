#pragma once
#include "vector.h"
#include "ray.h"

class Surface {
    public:
        Vector s;  // surface colour
        virtual bool hit(Ray r, Vector& n, Vector& p, float& t) { return false; };
};
