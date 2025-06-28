#include "vv_distance.hpp"
#include "distance_utils.hpp"

namespace c5d {

Scalar vv_distance(const Vector3S &x0, const Vector3S &x1) {
    return (x1 - x0).squaredNorm();
}

}  // namespace abd