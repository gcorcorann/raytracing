#pragma once

class Vector {
    public:
        float x;
        float y;
        float z;
        float dot(Vector o) { return x * o.x + y * o.y + z * o.z; }
        Vector mul(Vector o) {
            Vector r = {o.x * x, o.y * y, o.z * z};
            return r;
        }
        Vector sub(Vector o) {
            Vector r = {x - o.x, y - o.y, z - o.z};
            return r;
        }
        Vector add(Vector o) {
            Vector r = {o.x + x, o.y + y, o.z + z};
            return r;
        }
        Vector scale(float t) {
            Vector r = {x * t, y * t, z * t};
            return r;
        }
        float norm() { return sqrt(dot(*this)); }
        Vector unit() { return scale(1 / norm()); }
};
