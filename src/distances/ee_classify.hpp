#pragma once

#include "../gjk.hpp"

namespace c5d {

namespace EEType {
enum {
    VV02 = 0,
    VV03 = 1,
    VV12 = 2,
    VV13 = 3,
    VE023 = 4,
    VE123 = 5,
    VE201 = 6,
    VE301 = 7,
    EE = 8,
};
}  // namespace EEType

int ee_classify(const Vector3S &x0, const Vector3S &x1, const Vector3S &x2,
                const Vector3S &x3, Scalar &t0, Scalar &t1,
                Scalar parallel_sin_sq_tol = 1e-10);

}  // namespace abd
