#pragma once

#include "../gjk.hpp"

namespace c5d {

namespace VTType {
enum {
    VV01 = 0,
    VV02 = 1,
    VV03 = 2,
    VE012 = 3,
    VE023 = 4,
    VE031 = 5,
    VT = 6,
};
}  // namespace EEType

int vt_classify(const Vector3S &x0, const Vector3S &x1, const Vector3S &x2,
                const Vector3S &x3, Scalar &t2, Scalar &t3);

}  // namespace abd
