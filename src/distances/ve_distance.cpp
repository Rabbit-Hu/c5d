#include <iostream>

#include "../gjk.hpp"
#include "distance_utils.hpp"

namespace c5d {

Scalar ve_distance(const Vector3S &x0, const Vector3S &x1, const Vector3S &x2) {
    Vector3S a = x1 - x0;
    Vector3S b = x2 - x0;
    Vector3S c = x1 - x2;
    Vector3S e = a.cross(b);
    Scalar E = e.dot(e);
    Scalar F = c.dot(c);
    Scalar D = E / F;
    return D;
}

}  // namespace abd
