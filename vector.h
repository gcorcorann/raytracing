#pragma once
#include <cmath>
#include <iostream>

class Vector {
    public:
        float x;
        float y;
        float z;
        Vector () { };
        Vector(float a, float b, float c) {
            x = a;
            y = b;
            z = c;
        }
        float dot(Vector o) { return x * o.x + y * o.y + z * o.z; }
        Vector scale(float t) {
            Vector r = {x * t, y * t, z * t};
            return r;
        }
        float norm() { return sqrt(dot(*this)); }
        Vector unit() { return scale(1 / norm()); }
        Vector operator+ (Vector o) {
            Vector r = {x + o.x, y + o.y, z + o.z};
            return r;
        }
        Vector operator- (Vector o) {
            Vector r = {x - o.x, y - o.y, z - o.z};
            return r;
        }
        Vector operator* (Vector o) {
            Vector r = {x * o.x, y * o.y, z * o.z};
            return r;
        }
};
