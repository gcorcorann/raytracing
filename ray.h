#pragma once
#include "vector.h"

struct Ray {
    Vector e;  // origin
    Vector d;  // direction

    Ray() {};
    Ray(Vector origin, Vector direction) {
        e = origin;
        d = direction.unit();
        e = step();
    }
    Vector step() {
        return e + d.scale(1e-5f);  // avoid numeric impresisions
    }
};
