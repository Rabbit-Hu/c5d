#include "vt_distance.hpp"

#include "distance_utils.hpp"

namespace c5d {

Scalar vt_distance(const Vector3S &x0, const Vector3S &x1, const Vector3S &x2,
                   const Vector3S &x3) {
    Vector3S a = x0 - x3;
    Vector3S b = x1 - x3;
    Vector3S c = x2 - x3;
    Vector3S e = b.cross(c);
    Scalar E = e.dot(e);
    Scalar F = a.dot(e);
    Scalar G = F * F;
    Scalar D = G / E;
    return D;
}

}  // namespace abd
