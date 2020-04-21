#pragma once
#include "vector.h"
#include "ray.h"
#include "material.h"

class Surface {
    public:
        Material material;
        virtual bool hit(Ray r, Vector& n, Vector& p, float& t) { return false; };
};
