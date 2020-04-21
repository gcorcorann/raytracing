#pragma once
#include "vector.h"

class Material {
    public:
        Vector albedo;
        Vector diffuse_colour;
        float specular_exponent;
        Material() {}
        Material(Vector a, Vector colour, float spec) {
            albedo = a;
            diffuse_colour = colour;
            specular_exponent = spec;
        }
};
