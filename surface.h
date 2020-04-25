#pragma once
#include "vector.h"
#include "ray.h"
#include "material.h"

class Surface {
    public:
        Material m_material;
        virtual bool hit(Ray r, Vector& n, Vector& p, float& t) = 0;
};
