#include "scene.h"
#include "sphere.h"
#include "light.h"
#include "material.h"

int main() {
    const int width = 1024;
    const int height = 768;
    float focal_length = 1;  // only for perspective projection
    float fov = M_PI / 2;
    int max_depth = 4;
    Material ivory {{0.6, 0.3, 0.1}, {0.4, 0.4, 0.3}, 50.};
    Material red_rubber {{0.9, 0.1, 0.0}, {0.3, 0.1, 0.1}, 10.};
    Material mirror {{0.0, 3.0, 0.8}, {1.0, 1.0, 1.0}, 1425.};
    Sphere spheres [] {{{-3, 0, -16}, 2, ivory},  // centre, radius, material
                       {{-1, -1.5, -12}, 2, mirror},
                       {{1.5, -0.5, -18}, 3, red_rubber},
                       {{7, 5, -18}, 4, mirror}};
    int ns = sizeof(spheres) / sizeof(*spheres);
    Light lights [] {{{-20, 20, 20}, 1.5},
                     {{30, 50, -25}, 1.8},
                     {{30, 20, 30}, 1.7}};
    int nl = sizeof(lights) / sizeof(*lights);
    Scene scene (width, height, fov, focal_length, max_depth);
    scene.addSurfaces(spheres, ns);
    scene.addLights(lights, nl);
    scene.render();
    return 0;
}
