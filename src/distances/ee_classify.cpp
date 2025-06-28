#include "ee_classify.hpp"

namespace c5d {

int ee_classify(const Vector3S &x0, const Vector3S &x1, const Vector3S &x2,
                const Vector3S &x3, Scalar &t0, Scalar &t1,
                Scalar parallel_sin_sq_tol) {
    auto e0 = x1 - x0;
    auto e1 = x3 - x2;
    auto x02 = x2 - x0;

    // Find closest point on the two lines:
    //   h0 = x0 + t0 * e0
    //   h1 = x2 + t1 * e1

    auto A = e0.squaredNorm();
    auto B = -e0.dot(e1);
    auto C = e1.squaredNorm();

    auto D = x02.dot(e0);
    auto E = -x02.dot(e1);

    // [A B; B C] @ [t0; t1] = [D; E]
    auto D0 = A * C - B * B;    // det(A, B; B, C), always >= 0
    auto N0 = (C * D - B * E);  // t0 * D, "Nominator 0"
    auto N1 = (A * E - B * D);  // t1 * E, "Nominator 1"
    auto D1 = D0;               // Denominator 1

    int type = EEType::EE;
    bool no_EE = e0.cross(e1).squaredNorm() < parallel_sin_sq_tol * A * C;

    if (N0 <= 0.0 ||
        (no_EE && N0 < D0 * 0.5))  // the t0=0 edge of the square is visible
    {
        N1 = E;
        D1 = C;
        t0 = 0.0;
        type = EEType::VE023;  // if 0<t1<1 then return type = VE023
    } else if (N0 >= D0 ||
               (no_EE &&
                N0 >= D0 * 0.5))  // the t0=1 edge of the square is visible
    {
        N1 = E - B;
        D1 = C;
        t0 = 1.0;
        type = EEType::VE123;  // if 0<t1<1 then return type = VE123
    } else {
        t0 = N0 / D0;  // type = EEType::EE
    }

    if (N1 <= 0.0)  // determined t1=0, i.e. vertex 2 active
    {
        t1 = 0.0;
        if (D <= 0.0) {
            t0 = 0.0;
            return EEType::VV02;
        } else if (D >= A) {
            t0 = 1.0;
            return EEType::VV12;
        } else {
            t0 = D / A;
            return EEType::VE201;
        }
    } else if (N1 >= D1)  // determined t1=1, i.e. vertex 3 active
    {
        t1 = 1.0;
        if (D - B <= 0.0) {
            t0 = 0.0;
            return EEType::VV03;
        } else if (D - B >= A) {
            t0 = 1.0;
            return EEType::VV13;
        } else {
            t0 = (D - B) / A;
            return EEType::VE301;
        }
    }

    t1 = N1 / D1;
    return type;
}

}  // namespace abd