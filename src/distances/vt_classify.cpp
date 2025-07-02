#include "vt_classify.hpp"

namespace c5d {

int vt_classify(const Vector3S &x0, const Vector3S &x1, const Vector3S &x2,
                const Vector3S &x3, Scalar &t2, Scalar &t3) {
    Vector3S x01 = x1 - x0;
    Vector3S x12 = x2 - x1;
    Vector3S x13 = x3 - x1;

    if ((x12.cross(x13)).squaredNorm() < 1e-18) {
        // Check for 3 points overlapping
        if (x12.squaredNorm() < 1e-18 && x13.squaredNorm() < 1e-18) {
            t2 = 0.0;
            t3 = 0.0;
            return VTType::VV01;
        }

        // Fall back to VE distance
        if (x12.dot(x13) <= 0.0) {
            // VE023, t1 = 0
            Vector3S e = x3 - x2;
            t3 = (x0 - x3).dot(e) / e.squaredNorm();
            if (t3 >= 1.0) {
                t3 = 1.0;
                t2 = 0.0;
                return VTType::VV03;
            } else if (t3 <= 0.0) {
                t3 = 0.0;
                t2 = 1.0;
                return VTType::VV02;
            } else {
                t2 = 1.0 - t3;
                return VTType::VE023;
            }
        }
        else if ((x1 - x2).dot(x3 - x2) <= 0.0) {
            // VE031, t2 = 0
            t2 = 0.0;
            Vector3S e = x3 - x1;
            t3 = (x0 - x1).dot(e) / e.squaredNorm();
            if (t3 >= 1.0) {
                t3 = 1.0;
                return VTType::VV03;
            } else if (t3 <= 0.0) {
                t3 = 0.0;
                return VTType::VV01;
            } else {
                return VTType::VE031;
            }
        }
        else {
            // VE012, t3 = 0
            t3 = 0.0;
            Vector3S e = x2 - x1;
            t2 = (x0 - x1).dot(e) / e.squaredNorm();
            if (t2 >= 1.0) {
                t2 = 1.0;
                return VTType::VV02;
            } else if (t2 <= 0.0) {
                t2 = 0.0;
                return VTType::VV01;
            } else {
                return VTType::VE012;
            }
        }
    }

    Vector3S n = x12.cross(x13).normalized();
    Vector3S x0_proj = x0 - n.dot(x01) * n;
    Vector3S x10 = x0_proj - x1;
    Vector3S x20 = x0_proj - x2;
    Vector3S x30 = x0_proj - x3;

    Scalar a00 = x12.dot(x12);
    Scalar a01 = x12.dot(x13);
    Scalar a11 = x13.dot(x13);
    Scalar b0 = x10.dot(x12);
    Scalar b1 = x10.dot(x13);
    Scalar det_A = a00 * a11 - a01 * a01;
    t2 = b0 * a11 - b1 * a01;
    t3 = a00 * b1 - b0 * a01;
    Scalar t1 = det_A - t2 - t3;

    int n_negative_bary = (t1 < 0.0) + (t2 < 0.0) + (t3 < 0.0);

    if (n_negative_bary == 0) {
        // x0_proj inside triangle, distance type is "point-triangle".
        t2 /= det_A;
        t3 /= det_A;
        return VTType::VT;
    } else if (n_negative_bary == 1.0) {
        // Edge region
        // Project x0 onto the corresponding edge,
        //     if on the edge, distance type is "point-edge";
        //     if on the side of one vertex, distance type is "point-point".
        if (t1 < 0.0) {  // Edge x2x3
            Vector3S x23 = x3 - x2;
            Scalar dot_res = x23.dot(x20);
            if (dot_res < 0.0) {  // On the side of x2
                t2 = 1.0;
                t3 = 0.0;
                return VTType::VV02;
            } else if (dot_res > x23.dot(x23)) {  // On the side of x3
                t2 = 0.0;
                t3 = 1.0;
                return VTType::VV03;
            } else {  // On edge x2x3
                t3 = dot_res / x23.dot(x23);
                t2 = 1.0 - t3;
                return VTType::VE023;
            }
        } else if (t2 < 0.0) {  // Edge x3x1
            Vector3S x31 = x1 - x3;
            Scalar dot_res = x31.dot(x30);
            t2 = 0.0;
            if (dot_res < 0.0) {  // On the side of x3
                t3 = 1.0;
                return VTType::VV03;
            } else if (dot_res > x31.dot(x31)) {  // On the side of x1
                t3 = 0.0;
                return VTType::VV01;
            } else {  // On edge x3x1
                t3 = 1.0 - dot_res / x31.dot(x31);
                return VTType::VE031;
            }
        } else {  // Edge x1x2
            Scalar dot_res = x12.dot(x10);
            t3 = 0.0;
            if (dot_res < 0.0) {  // On the side of x1
                t2 = 0.0;
                return VTType::VV01;
            } else if (dot_res > x12.dot(x12)) {  // On the side of x2
                t2 = 1.0;
                return VTType::VV02;
            } else {  // On edge x1x2
                t2 = dot_res / x12.dot(x12);
                return VTType::VE012;
            }
        }
    } else {
        if (t1 >= 0.0) {
            Scalar dot12 = x12.dot(x10);
            Scalar dot13 = x13.dot(x10);
            if (dot12 > 0.0 && dot13 <= 0.0) {  // On edge x1x2
                t3 = 0.0;
                Scalar x12_sq = x12.dot(x12);
                if (dot12 > x12_sq) {
                    t2 = 1.0;
                    return VTType::VV02;
                } else {
                    t2 = dot12 / x12_sq;
                    return VTType::VE012;
                }
            } else if (dot13 > 0.0 && dot12 <= 0.0) {  // On edge x1x3
                t2 = 0.0;
                Scalar x13_sq = x13.dot(x13);
                if (dot13 > x13_sq) {
                    t3 = 1.0;
                    return VTType::VV03;
                } else {
                    t3 = dot13 / x13_sq;
                    return VTType::VE031;
                }
            } else {
                t2 = 0.0;
                t3 = 0.0;
                return VTType::VV01;
            }
        } else if (t2 >= 0.0) {
            Vector3S x23 = x3 - x2;
            Vector3S x21 = x1 - x2;
            Scalar dot23 = x23.dot(x20);
            Scalar dot21 = x21.dot(x20);
            if (dot23 > 0.0 && dot21 <= 0.0) {  // On edge x2x3
                Scalar x23_sq = x23.dot(x23);
                if (dot23 > x23_sq) {
                    t2 = 0.0;
                    t3 = 1.0;
                    return VTType::VV03;
                } else {
                    t3 = dot23 / x23_sq;
                    t2 = 1.0 - t3;
                    return VTType::VE023;
                }
            } else if (dot21 > 0.0 && dot23 <= 0.0) {  // On edge x2x1
                t3 = 0.0;
                Scalar x21_sq = x21.dot(x21);
                if (dot21 > x21_sq) {
                    t2 = 0.0;
                    return VTType::VV01;
                } else {
                    t2 = 1.0 - dot21 / x21_sq;
                    return VTType::VE012;
                }
            } else {
                t2 = 1.0;
                t3 = 0.0;
                return VTType::VV02;
            }
        } else {
            Vector3S x31 = x1 - x3;
            Vector3S x32 = x2 - x3;
            Scalar dot31 = x31.dot(x30);
            Scalar dot32 = x32.dot(x30);
            if (dot31 > 0.0 && dot32 <= 0.0) {  // On edge x3x1
                t2 = 0.0;
                Scalar x31_sq = x31.dot(x31);
                if (dot31 > x31_sq) {
                    t3 = 0.0;
                    return VTType::VV01;
                } else {
                    t3 = 1.0 - dot31 / x31_sq;
                    return VTType::VE031;
                }
            } else if (dot32 > 0.0 && dot31 <= 0.0) {  // On edge x3x2
                Scalar x32_sq = x32.dot(x32);
                if (dot32 > x32_sq) {
                    t2 = 1.0;
                    t3 = 0.0;
                    return VTType::VV02;
                } else {
                    t2 = dot32 / x32_sq;
                    t3 = 1.0 - t2;
                    return VTType::VE023;
                }
            } else {
                t2 = 0.0;
                t3 = 1.0;
                return VTType::VV03;
            }
        }
    }
}

}  // namespace abd