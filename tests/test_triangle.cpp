#include <iostream>
#include <cassert>

#include "../triangle.h"
#include "../ray.h"
#include "../vector.h"

int main() {
    Vector a = {-100, 0, -10};
    Vector b = {100, 0, -10};
    Vector c = {0, 100, -10};
    Vector colour = {1, 0, 0};
    Triangle tri = {a, b, c, colour};

    Vector e = {0, 10, 0};
    Vector d = {0, 0, -1};
    Ray ray = {e, d};

    Vector n = {0, 0, 0};
    Vector p = {0, 0, 0};
    float t;
    bool hit = tri.hit(ray, n, p, t);
    assert (hit == true);
    assert (t == 10);
    return 0;
}
