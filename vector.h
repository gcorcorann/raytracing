#pragma once
#include <cmath>

class Vector {
    public:
        float x;
        float y;
        float z;
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
