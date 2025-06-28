#include "ACCD.hpp"

// #include <omp.h>
#include <spdlog/spdlog.h>

#include "distances/distances.hpp"

namespace c5d {

namespace exp {

static Scalar VTAdditiveCCD(Vector3S x0, Vector3S x1, Vector3S x2, Vector3S x3,
                            Vector3S v0, Vector3S v1, Vector3S v2, Vector3S v3,
                            Scalar slackness, Scalar thickness, int max_iter,
                            int &num_iter) {
    Vector3S v_avg = (v0 + v1 + v2 + v3) / 4.0;
    v0 -= v_avg;
    v1 -= v_avg;
    v2 -= v_avg;
    v3 -= v_avg;

    Scalar l_v =
        v0.norm() + std::max(v1.norm(), std::max(v2.norm(), v3.norm()));

    Scalar t2, t3;
    int vt_type = vt_classify(x0, x1, x2, x3, t2, t3);
    Vector3S h = x2 * t2 + x3 * t3 + x1 * (1.0 - t2 - t3);
    Scalar d_sqr = (x0 - h).squaredNorm();
    Scalar d_minus_thickness =
        (d_sqr - thickness * thickness) / (sqrt(d_sqr) + thickness);

    if (d_minus_thickness <= 0.0) {
        spdlog::warn("VTCollision CCD d_minus_thickness = {} <= 0.0",
                     d_minus_thickness);
        return 0.0;
    }

    Scalar g = (1.0 - slackness) * d_minus_thickness;
    Scalar t = 0.0;
    Scalar t_l = slackness * d_minus_thickness / l_v;

    int i_iter = 0;
    while (i_iter < max_iter) {
        i_iter++;

        x0 += t_l * v0;
        x1 += t_l * v1;
        x2 += t_l * v2;
        x3 += t_l * v3;

        vt_type = vt_classify(x0, x1, x2, x3, t2, t3);
        h = x2 * t2 + x3 * t3 + x1 * (1.0 - t2 - t3);
        d_sqr = (x0 - h).squaredNorm();
        d_minus_thickness =
            (d_sqr - thickness * thickness) / (sqrt(d_sqr) + thickness);
        if (i_iter > 1 && d_minus_thickness <= g) {
            break;
        }

        t += t_l;
        if (t >= 1.0) {
            break;
        }
        t_l = 0.9 * d_minus_thickness / l_v;
    }

    num_iter = i_iter;

    if (i_iter == max_iter) {
        spdlog::warn("VTCollision CCD did not converge");
    }

    return t;
}

static Scalar EEAdditiveCCD(Vector3S x0, Vector3S x1, Vector3S x2, Vector3S x3,
                            Vector3S v0, Vector3S v1, Vector3S v2, Vector3S v3,
                            Scalar slackness, Scalar thickness, int max_iter,
                            int &num_iter) {
    Vector3S v_avg = (v0 + v1 + v2 + v3) / 4.0;
    v0 -= v_avg;
    v1 -= v_avg;
    v2 -= v_avg;
    v3 -= v_avg;

    Scalar l_v =
        std::max(v0.norm(), v1.norm()) + std::max(v2.norm(), v3.norm());

    Scalar t0, t1;
    int ee_type = ee_classify(x0, x1, x2, x3, t0, t1);
    Vector3S h0 = x0 + t0 * (x1 - x0);
    Vector3S h1 = x2 + t1 * (x3 - x2);
    Scalar d_sqr = (h0 - h1).squaredNorm();
    Scalar d_minus_thickness =
        (d_sqr - thickness * thickness) / (sqrt(d_sqr) + thickness);

    if (d_minus_thickness <= 0.0) {
        spdlog::warn("EECollision CCD d_minus_thickness = {} <= 0.0",
                     d_minus_thickness);
        return 0.0;
    }

    Scalar g = (1.0 - slackness) * d_minus_thickness;
    Scalar t = 0.0;
    Scalar t_l = slackness * d_minus_thickness / l_v;

    int i_iter = 0;
    while (i_iter < max_iter) {
        i_iter++;

        x0 += t_l * v0;
        x1 += t_l * v1;
        x2 += t_l * v2;
        x3 += t_l * v3;

        ee_type = ee_classify(x0, x1, x2, x3, t0, t1);
        h0 = x0 + t0 * (x1 - x0);
        h1 = x2 + t1 * (x3 - x2);
        d_sqr = (h0 - h1).squaredNorm();
        d_minus_thickness =
            (d_sqr - thickness * thickness) / (sqrt(d_sqr) + thickness);
        if (i_iter > 1 && d_minus_thickness <= g) {
            break;
        }

        t += t_l;
        if (t >= 1.0) {
            break;
        }
        t_l = 0.9 * d_minus_thickness / l_v;
    }

    num_iter = i_iter;

    if (i_iter == max_iter) {
        spdlog::warn("EECollision CCD did not converge");
    }

    return t;
}

Scalar additiveCCD(const MatrixXS &V0, const MatrixXS &V1, const MatrixXi &F0,
                   const MatrixXi &F1, const Matrix3S &P0, const Matrix3S &P1,
                   const Vector3S &x0, const Vector3S &x1, const Matrix3S &A0,
                   const Matrix3S &A1, const Vector3S &v0, const Vector3S &v1,
                   Scalar slackness, Scalar thickness, int max_iter,
                   double &num_iter) {
    Scalar step_size = 1.0;

    MatrixXi E0, E1;
    igl::edges(F0, E0);
    igl::edges(F1, E1);

    double sum_num_iter = 0.0;
    double max_num_iter = 0.0;
    int num_iter_i = 0;

    for (int v_i = 0; v_i < V0.rows(); v_i++) {            // vertex id in V0
        Vector3S x_0 = P0 * V0.row(v_i).transpose() + x0;  // vertex position
        Vector3S v_0 = A0 * V0.row(v_i).transpose() + v0;  // vertex velocity

        for (int t_i = 0; t_i < F1.rows(); t_i++) {  // triangle id in F1
            Vector3S x_1 = P1 * V1.row(F1(t_i, 0)).transpose() + x1;
            Vector3S v_1 = A1 * V1.row(F1(t_i, 0)).transpose() + v1;
            Vector3S x_2 = P1 * V1.row(F1(t_i, 1)).transpose() + x1;
            Vector3S v_2 = A1 * V1.row(F1(t_i, 1)).transpose() + v1;
            Vector3S x_3 = P1 * V1.row(F1(t_i, 2)).transpose() + x1;
            Vector3S v_3 = A1 * V1.row(F1(t_i, 2)).transpose() + v1;

            Scalar new_step_size =
                VTAdditiveCCD(x_0, x_1, x_2, x_3, v_0, v_1, v_2, v_3, slackness,
                              thickness, max_iter, num_iter_i);
            step_size = std::min(step_size, new_step_size);
            sum_num_iter += num_iter_i;
            max_num_iter = std::max(max_num_iter, (double)num_iter_i);
        }
    }

    for (int v_i = 0; v_i < V1.rows(); v_i++) {            // vertex id in V1
        Vector3S x_0 = P1 * V1.row(v_i).transpose() + x1;  // vertex position
        Vector3S v_0 = A1 * V1.row(v_i).transpose() + v1;  // vertex velocity

        for (int t_i = 0; t_i < F0.rows(); t_i++) {  // triangle id in F0
            Vector3S x_1 = P0 * V0.row(F0(t_i, 0)).transpose() + x0;
            Vector3S v_1 = A0 * V0.row(F0(t_i, 0)).transpose() + v0;
            Vector3S x_2 = P0 * V0.row(F0(t_i, 1)).transpose() + x0;
            Vector3S v_2 = A0 * V0.row(F0(t_i, 1)).transpose() + v0;
            Vector3S x_3 = P0 * V0.row(F0(t_i, 2)).transpose() + x0;
            Vector3S v_3 = A0 * V0.row(F0(t_i, 2)).transpose() + v0;

            Scalar new_step_size =
                VTAdditiveCCD(x_0, x_1, x_2, x_3, v_0, v_1, v_2, v_3, slackness,
                              thickness, max_iter, num_iter_i);
            step_size = std::min(step_size, new_step_size);
            sum_num_iter += num_iter_i;
            max_num_iter = std::max(max_num_iter, (double)num_iter_i);
        }
    }

    for (int e_i = 0; e_i < E0.rows(); e_i++) {
        Vector3S x_0 = P0 * V0.row(E0(e_i, 0)).transpose() + x0;
        Vector3S v_0 = A0 * V0.row(E0(e_i, 0)).transpose() + v0;
        Vector3S x_1 = P0 * V0.row(E0(e_i, 1)).transpose() + x0;
        Vector3S v_1 = A0 * V0.row(E0(e_i, 1)).transpose() + v0;

        for (int e_j = 0; e_j < E1.rows(); e_j++) {
            Vector3S x_2 = P1 * V1.row(E1(e_j, 0)).transpose() + x1;
            Vector3S v_2 = A1 * V1.row(E1(e_j, 0)).transpose() + v1;
            Vector3S x_3 = P1 * V1.row(E1(e_j, 1)).transpose() + x1;
            Vector3S v_3 = A1 * V1.row(E1(e_j, 1)).transpose() + v1;

            Scalar new_step_size =
                EEAdditiveCCD(x_0, x_1, x_2, x_3, v_0, v_1, v_2, v_3, slackness,
                              thickness, max_iter, num_iter_i);
            step_size = std::min(step_size, new_step_size);
            sum_num_iter += num_iter_i;
            max_num_iter = std::max(max_num_iter, (double)num_iter_i);
        }
    }

    // num_iter = sum_num_iter / (V0.rows() * F1.rows() + V1.rows() * F0.rows()
    // +
    //    E0.rows() * E1.rows());
    num_iter = max_num_iter;

    return step_size;
}

// Scalar additiveCCD_parallel(const MatrixXS &V0, const MatrixXS &V1,
//                             const MatrixXi &F0, const MatrixXi &F1,
//                             const Matrix3S &P0, const Matrix3S &P1,
//                             const Vector3S &x0, const Vector3S &x1,
//                             const Matrix3S &A0, const Matrix3S &A1,
//                             const Vector3S &v0, const Vector3S &v1,
//                             Scalar slackness, Scalar thickness, int max_iter,
//                             double &num_iter) {
//     Scalar step_size = 1.0;

//     MatrixXi E0, E1;
//     igl::edges(F0, E0);
//     igl::edges(F1, E1);

//     double sum_num_iter = 0.0;
//     double max_num_iter = 0.0;

//     const int num_vt_threads = V0.rows() * F1.rows();
//     const int num_tv_threads = V1.rows() * F0.rows();
//     const int num_ee_threads = E0.rows() * E1.rows();

// #pragma omp parallel for reduction(min : step_size) \
//     reduction(+ : sum_num_iter) reduction(max : max_num_iter)
//     for (int vt_i = 0; vt_i < num_vt_threads; vt_i++) {
//         int v_i = vt_i / F1.rows();
//         int t_i = vt_i % F1.rows();

//         Vector3S x_0 = P0 * V0.row(v_i).transpose() + x0;  // vertex position
//         Vector3S v_0 = A0 * V0.row(v_i).transpose() + v0;  // vertex velocity
//         Vector3S x_1 = P1 * V1.row(F1(t_i, 0)).transpose() + x1;
//         Vector3S v_1 = A1 * V1.row(F1(t_i, 0)).transpose() + v1;
//         Vector3S x_2 = P1 * V1.row(F1(t_i, 1)).transpose() + x1;
//         Vector3S v_2 = A1 * V1.row(F1(t_i, 1)).transpose() + v1;
//         Vector3S x_3 = P1 * V1.row(F1(t_i, 2)).transpose() + x1;
//         Vector3S v_3 = A1 * V1.row(F1(t_i, 2)).transpose() + v1;

//         int num_iter_i = 0;
//         Scalar new_step_size =
//             VTAdditiveCCD(x_0, x_1, x_2, x_3, v_0, v_1, v_2, v_3, slackness,
//                           thickness, max_iter, num_iter_i);

//         step_size = step_size < new_step_size ? step_size : new_step_size;
//         sum_num_iter += (double)num_iter_i;
//         max_num_iter = max_num_iter > (double)num_iter_i ? max_num_iter
//                                                          : (double)num_iter_i;

//         // #pragma omp critical
//         //         {
//         //             step_size = step_size < new_step_size ? step_size :
//         //             new_step_size; sum_num_iter += (double)num_iter_i;
//         //             max_num_iter = max_num_iter > (double)num_iter_i
//         //                                ? max_num_iter
//         //                                : (double)num_iter_i;
//         //         }
//     }

// #pragma omp parallel for reduction(min : step_size) \
//     reduction(+ : sum_num_iter) reduction(max : max_num_iter)
//     for (int vt_i = 0; vt_i < num_tv_threads; vt_i++) {
//         int v_i = vt_i / F0.rows();
//         int t_i = vt_i % F0.rows();

//         Vector3S x_0 = P1 * V1.row(v_i).transpose() + x1;  // vertex position
//         Vector3S v_0 = A1 * V1.row(v_i).transpose() + v1;  // vertex velocity
//         Vector3S x_1 = P0 * V0.row(F0(t_i, 0)).transpose() + x0;
//         Vector3S v_1 = A0 * V0.row(F0(t_i, 0)).transpose() + v0;
//         Vector3S x_2 = P0 * V0.row(F0(t_i, 1)).transpose() + x0;
//         Vector3S v_2 = A0 * V0.row(F0(t_i, 1)).transpose() + v0;
//         Vector3S x_3 = P0 * V0.row(F0(t_i, 2)).transpose() + x0;
//         Vector3S v_3 = A0 * V0.row(F0(t_i, 2)).transpose() + v0;

//         int num_iter_i = 0;
//         Scalar new_step_size =
//             VTAdditiveCCD(x_0, x_1, x_2, x_3, v_0, v_1, v_2, v_3, slackness,
//                           thickness, max_iter, num_iter_i);

//         step_size = step_size < new_step_size ? step_size : new_step_size;
//         sum_num_iter += (double)num_iter_i;
//         max_num_iter = max_num_iter > (double)num_iter_i ? max_num_iter
//                                                          : (double)num_iter_i;

//         // #pragma omp critical
//         //         {
//         //             step_size = step_size < new_step_size ? step_size :
//         //             new_step_size; sum_num_iter += (double)num_iter_i;
//         //             max_num_iter = max_num_iter > (double)num_iter_i
//         //                                ? max_num_iter
//         //                                : (double)num_iter_i;
//         //         }
//     }

// #pragma omp parallel for reduction(min : step_size) \
//     reduction(+ : sum_num_iter) reduction(max : max_num_iter)
//     for (int ee_i = 0; ee_i < num_ee_threads; ee_i++) {
//         int e_i = ee_i / E1.rows();
//         int e_j = ee_i % E1.rows();

//         Vector3S x_0 = P0 * V0.row(E0(e_i, 0)).transpose() + x0;
//         Vector3S v_0 = A0 * V0.row(E0(e_i, 0)).transpose() + v0;
//         Vector3S x_1 = P0 * V0.row(E0(e_i, 1)).transpose() + x0;
//         Vector3S v_1 = A0 * V0.row(E0(e_i, 1)).transpose() + v0;
//         Vector3S x_2 = P1 * V1.row(E1(e_j, 0)).transpose() + x1;
//         Vector3S v_2 = A1 * V1.row(E1(e_j, 0)).transpose() + v1;
//         Vector3S x_3 = P1 * V1.row(E1(e_j, 1)).transpose() + x1;
//         Vector3S v_3 = A1 * V1.row(E1(e_j, 1)).transpose() + v1;

//         int num_iter_i = 0;
//         Scalar new_step_size =
//             EEAdditiveCCD(x_0, x_1, x_2, x_3, v_0, v_1, v_2, v_3, slackness,
//                           thickness, max_iter, num_iter_i);

//         step_size = step_size < new_step_size ? step_size : new_step_size;
//         sum_num_iter += (double)num_iter_i;
//         max_num_iter = max_num_iter > (double)num_iter_i ? max_num_iter
//                                                          : (double)num_iter_i;

//         // #pragma omp critical
//         //         {
//         //             step_size = step_size < new_step_size ? step_size :
//         //             new_step_size; sum_num_iter += (double)num_iter_i;
//         //             max_num_iter = max_num_iter > (double)num_iter_i
//         //                                ? max_num_iter
//         //                                : (double)num_iter_i;
//         //         }
//     }

//     num_iter = max_num_iter;

//     return step_size;
// }

}  // namespace exp

}  // namespace c5d