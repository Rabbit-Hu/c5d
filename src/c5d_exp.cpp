#include "c5d_exp.hpp"

#include <igl/writePLY.h>
#include <spdlog/spdlog.h>

#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>

namespace c5d {

namespace exp {

Scalar c5d_full_gjk_linear(const Polytope *B0, const Polytope *B1,
                           const Matrix3S &P0, const Matrix3S &P1,
                           const Vector3S &x0, const Vector3S &x1,
                           const Matrix3S &A0, const Matrix3S &A1,
                           const Vector3S &v0, const Vector3S &v1,
                           Scalar slackness, Scalar thickness, int max_iter,
                           int &num_iter, Scalar *t_history) {
    Scalar d, nu, dt, d_minus_thickness;
    Vector3S n, ds0, ds1, u;

    d = gjk_transformed(B0, B1, P0, P1, x0, x1, n);
    d_minus_thickness = (d - thickness * thickness) / (sqrt(d) + thickness);
    assert(d_minus_thickness > 0.0);
    Scalar g = (1.0 - slackness) * d_minus_thickness;

    Scalar t = 0.0;
    int k = 0;
    while (k < max_iter) {
        k++;

        d = gjk_transformed(B0, B1, P0 + t * A0, P1 + t * A1, x0 + t * v0,
                            x1 + t * v1, n);
        n = -n;

        d_minus_thickness = (d - thickness * thickness) / (sqrt(d) + thickness);
        if (k > 0 && d_minus_thickness <= g) {
            break;
        }

        ds0 = support_transformed(B0, A0, v0, n);
        ds1 = support_transformed(B1, A1, v1, -n);

        u = ds1 - ds0;
        nu = n.dot(u);
        if (nu >= 0.0) {
            t = 1.0;
            break;
        } else {
            dt = 0.9 * d_minus_thickness / (-nu / sqrt(d));

            t += dt;
            t = std::min(t, (Scalar)1.0);

            if (t >= 1.0) {
                break;
            }
        }

        if (t_history != nullptr) {
            t_history[k] = t;
        }
    }

    num_iter = k;

    if (k == max_iter) {
        spdlog::warn("c5d_full_gjk_linear did not converge after {} iterations",
                     max_iter);
    }

    return t;
}

static inline Scalar spectral_norm(const Matrix3S &A) {
    // return largest singular value
    return A.jacobiSvd().singularValues()(0);
}

Scalar c5d_full_gjk_quad(const Polytope *B0, const Polytope *B1,
                         const Matrix3S &P0, const Matrix3S &P1,
                         const Vector3S &x0, const Vector3S &x1,
                         const Matrix3S &A0, const Matrix3S &A1,
                         const Vector3S &v0, const Vector3S &v1,
                         Scalar slackness, Scalar thickness, Scalar K,
                         int max_iter, int &num_iter, Scalar *t_history) {
    Matrix3S P0_inv = P0.inverse();
    Matrix3S P1_inv = P1.inverse();
    Matrix3S A0_new = A0 * P0_inv;       // A0'
    Matrix3S A1_new = A1 * P1_inv;       // A1'
    Vector3S v0_new = v0 - A0_new * x0;  // v0'
    Vector3S v1_new = v1 - A1_new * x1;  // v1'

    // B' = B(0) = P B + x
    // B(t) = (I + A' t) B' + v' t
    //      = (I + A' t) (P B + x) + v' t
    //      = (P + A t) B + (x + v t)

    // To swap every iter
    const Polytope *B0_ptr = B0, *B1_ptr = B1;
    Matrix3S *A0_ptr = &A0_new, *A1_ptr = &A1_new;
    const Matrix3S *P0_ptr = &P0, *P1_ptr = &P1;
    Vector3S *v0_ptr = &v0_new, *v1_ptr = &v1_new;
    const Vector3S *x0_ptr = &x0, *x1_ptr = &x1;

    Scalar d, t_clip, rho0, rho1, R, r0, dt, d_minus_thickness, v_max,
        delta;
    Vector3S n;
    Matrix3S P0_t, P1_t, P0_t_inv, P1_t_inv, A0_t, A1_t;
    Vector3S x0_t, x1_t, v0_t, v1_t, ds1_rel;
    auto I = Matrix3S::Identity();

    d = gjk_transformed(B0, B1, P0, P1, x0, x1, n);
    d_minus_thickness = (d - thickness * thickness) / (sqrt(d) + thickness);
    assert(d_minus_thickness > 0.0);
    Scalar g = (1.0 - slackness) * d_minus_thickness;

    Scalar t = 0.0;
    int k = 0;
    while (k < max_iter) {
        k++;

        std::swap(B0_ptr, B1_ptr);
        std::swap(A0_ptr, A1_ptr);
        std::swap(P0_ptr, P1_ptr);
        std::swap(v0_ptr, v1_ptr);
        std::swap(x0_ptr, x1_ptr);

        P0_t = I + t * (*A0_ptr);
        P1_t = I + t * (*A1_ptr);
        x0_t = t * (*v0_ptr);
        x1_t = t * (*v1_ptr);

        // B(t) = (I + A' t) B' + v' t
        //      = P_t B' + x_t
        //      = P_t (P B + x) + x_t
        //      = (P_t P) B + (P_t x + x_t)
        //      = (I + t A') P B + (I + t A') x + x_t
        //      = (P + t A) B + (x + t (A' x + v'))
        //      = (P + t A) B + (x + t v)

        d = gjk_transformed(B0_ptr, B1_ptr, P0_t * *P0_ptr, P1_t * *P1_ptr,
                            P0_t * *x0_ptr + x0_t, P1_t * *x1_ptr + x1_t, n);
        d_minus_thickness = (d - thickness * thickness) / (sqrt(d) + thickness);
        if (d_minus_thickness <= g) {
            break;
        }
        n = -n;

        // B(t + dt) = (I + A' (t + dt)) B' + v' (t + dt)
        //           = (P_t + A' dt) B' + (x_t + v' dt)
        //           = (I + (A' P_t^-1) dt) (P_t B') + (x_t + v' dt)
        //           = (I + (A' P_t^-1) dt) (P_t B' + x_t)
        //             - (I + (A' P_t^-1) dt) x_t + (x_t + v' dt)
        //           = (I + (A' P_t^-1) dt) (P_t B' + x_t)
        //             + (v' - (A' P_t^-1) x_t) dt                ......(*)
        //           = (I + (A' P_t^-1) dt) (P_t B' + x_t) + (P_t^-1 v') dt
        // (*): v' - (A' P_t^-1) x_t = v' - t A' (I + t A')^-1 v'
        //                           = (I - t A' (I + t A')^-1) v'
        //                           = ((I + t A')(I + t A')^-1
        //                             - t A' (I + t A')^-1) v'
        //                           = (I + t A')^-1 v'

        P0_t_inv = P0_t.inverse();
        P1_t_inv = P1_t.inverse();
        A0_t = *A0_ptr * P0_t_inv;
        A1_t = *A1_ptr * P1_t_inv;
        v0_t = P0_t_inv * *v0_ptr;
        v1_t = P1_t_inv * *v1_ptr;

        rho0 = spectral_norm(A0_t);
        rho1 = spectral_norm(A1_t);
        t_clip = (K - 1.0) / (K * rho0 + EPSILON_ABS);

        // s0 = support_transformed(B0_ptr, P0_t * *P0_ptr, P0_t * *x0_ptr +
        // x0_t,
        //                          n);
        // dt = solve_pointwise_toi_quad(
        //     B1_ptr, P1_t * *P1_ptr, P1_t * *x1_ptr + x1_t, A1_t - A0_t,
        //     v1_t - v0_t, s0, -n, K * rho0 * (rho0 + rho1),
        //     K * rho0 * (v1_t - v0_t).norm(), thickness);
        // dt = std::min(dt, t_clip) * 0.9;

        r0 = 0.0;
        for (int i = 0, N = B0_ptr->numVertices(); i < N; i++) {
            r0 = std::max(
                r0, (P0_t * (*P0_ptr * B0_ptr->getVertex(i) + *x0_ptr) + x0_t)
                         .norm());
        }
        R = K * rho0 * ((rho0 + rho1) * r0 + (v1_t - v0_t).norm());
        ds1_rel = support_transformed(B1_ptr, (A1_t - A0_t) * P1_t * *P1_ptr,
                                      (A1_t - A0_t) * P1_t * *x1_ptr +
                                          (A1_t - A0_t) * x1_t + (v1_t - v0_t),
                                      -n);

        v_max = -n.normalized().dot(ds1_rel);
        delta = v_max * v_max + 4.0 * R * d_minus_thickness;
        if (abs(R) < 1e-9) {
            dt = d_minus_thickness / v_max;
        } else {
            dt = (-v_max + sqrt(delta)) / (2 * R);
        }
        dt = std::min(dt, t_clip) * 0.9;

        assert(spectral_norm((Matrix3S::Identity() + dt * A0_t).inverse()) <=
               K);

        t = std::min(t + dt, (Scalar)1.0);

        if (k > 9990) {
            spdlog::info(
                "r0 = {}, R = {}, v_max = {}, t_clip = {}, dt = {}, t "
                "= {}",
                r0, R, v_max, t_clip, dt, t);
        }

        if (t_history != nullptr) {
            t_history[k] = t;
        }

        if (t >= 1.0) {
            break;
        }
    }

    num_iter = k;

    if (k == max_iter) {
        spdlog::warn("c5d_full_gjk_quad did not converge after {} iterations",
                     max_iter);
    }

    assert(t >= 0.0);

    d = gjk_transformed(B0, B1, P0 + t * A0, P1 + t * A1, x0 + t * v0,
                        x1 + t * v1, n);
    if (d < thickness * thickness * 0.9 - 1e-16) {
        spdlog::warn(
            "c5d_full_gjk_quad_pw: d = {0:.6e} < thickness^2 = {1:.6e}", d,
            thickness * thickness);
        spdlog::info("k = {0}, t = {1:.6e}", k, t);
        exit(1);
    }

    return t;
}

static inline Scalar solve_pointwise_toi_quad(
    const Polytope *B, const Matrix3S &P, const Vector3S &x, const Matrix3S &A,
    const Vector3S &v, const Vector3S &s0, const Vector3S &n, Scalar R1,
    Scalar R0, Scalar thickness) {
    Vector3S n_unit = n.normalized();
    Vector3S A_T_n_unit = A.transpose() * n_unit;
    Scalar v_dot_n_unit = v.dot(n_unit);
    Scalar s0_dot_n_unit = s0.dot(n_unit) - thickness;
    Scalar min_dt = 1.0;
    Vector3S b;
    Scalar qa, qb, qc, delta, dt;
    for (int i = 0, N = B->numVertices(); i < N; i++) {
        b = P * B->getVertex(i) + x;
        qa = R1 * b.norm() + R0;
        qb = A_T_n_unit.dot(b) + v_dot_n_unit;
        qc = n_unit.dot(b) - s0_dot_n_unit;
        delta = qb * qb - 4.0 * qa * qc;

        if (qc > 0.0) {
            if (qc > 1e-15) {
                spdlog::warn("qc = {} > 0.0", qc);
                spdlog::info("qa = {}, qb = {}", qa, qb);
                spdlog::info("d (sqrt) = {}, (b - s0).dot(n) = {}", n.norm(),
                             (b - s0).dot(n));
                exit(1);
            }
            return 0.0;
        }

        assert(delta >= 0.0);
        assert(qa >= 0.0);
        assert(qc <= 0.0);

        const Scalar eps = EPSILON_REL;
        if (qa < eps) {
            dt = std::abs(qb) < eps ? 1.0 : -qc / qb;
        } else {
            dt = (-qb + sqrt(delta)) / (2.0 * qa);
        }
        if (dt >= 0.0 && dt < min_dt) {
            min_dt = dt;
        }
    }
    return min_dt;
}

Scalar c5d_full_gjk_quad_pw(const Polytope *B0, const Polytope *B1,
                            const Matrix3S &P0, const Matrix3S &P1,
                            const Vector3S &x0, const Vector3S &x1,
                            const Matrix3S &A0, const Matrix3S &A1,
                            const Vector3S &v0, const Vector3S &v1,
                            Scalar slackness, Scalar thickness, Scalar K,
                            int max_iter, int &num_iter, Scalar *t_history) {
    Matrix3S P0_inv = P0.inverse();
    Matrix3S P1_inv = P1.inverse();
    Matrix3S A0_new = A0 * P0_inv;       // A0'
    Matrix3S A1_new = A1 * P1_inv;       // A1'
    Vector3S v0_new = v0 - A0_new * x0;  // v0'
    Vector3S v1_new = v1 - A1_new * x1;  // v1'

    // B' = B(0) = P B + x
    // B(t) = (I + A' t) B' + v' t
    //      = (I + A' t) (P B + x) + v' t
    //      = (P + A t) B + (x + v t)

    // To swap every iter
    const Polytope *B0_ptr = B0, *B1_ptr = B1;
    Matrix3S *A0_ptr = &A0_new, *A1_ptr = &A1_new;
    const Matrix3S *P0_ptr = &P0, *P1_ptr = &P1;
    Vector3S *v0_ptr = &v0_new, *v1_ptr = &v1_new;
    const Vector3S *x0_ptr = &x0, *x1_ptr = &x1;

    Scalar d, t_clip, rho0, rho1, dt, d_minus_thickness;
    Vector3S n;
    Matrix3S P0_t, P1_t, P0_t_inv, P1_t_inv, A0_t, A1_t;
    Vector3S x0_t, x1_t, v0_t, v1_t, s0;
    auto I = Matrix3S::Identity();

    d = gjk_transformed(B0, B1, P0, P1, x0, x1, n);
    d_minus_thickness = (d - thickness * thickness) / (sqrt(d) + thickness);
    assert(d_minus_thickness > 0.0);
    Scalar g = (1.0 - slackness) * d_minus_thickness;

    Scalar t = 0.0;
    int k = 0;
    while (k < max_iter) {
        k++;
        
        std::swap(B0_ptr, B1_ptr);
        std::swap(A0_ptr, A1_ptr);
        std::swap(P0_ptr, P1_ptr);
        std::swap(v0_ptr, v1_ptr);
        std::swap(x0_ptr, x1_ptr);

        P0_t = I + t * (*A0_ptr);
        P1_t = I + t * (*A1_ptr);
        x0_t = t * (*v0_ptr);
        x1_t = t * (*v1_ptr);

        // B(t) = (I + A' t) B' + v' t
        //      = P_t B' + x_t
        //      = P_t (P B + x) + x_t
        //      = (P_t P) B + (P_t x + x_t)
        //      = (I + t A') P B + (I + t A') x + x_t
        //      = (P + t A) B + (x + t (A' x + v'))
        //      = (P + t A) B + (x + t v)

        d = gjk_transformed(B0_ptr, B1_ptr, P0_t * *P0_ptr, P1_t * *P1_ptr,
                            P0_t * *x0_ptr + x0_t, P1_t * *x1_ptr + x1_t, n);
        d_minus_thickness = (d - thickness * thickness) / (sqrt(d) + thickness);
        if (d_minus_thickness <= g) {
            break;
        }
        n = -n;

        // B(t + dt) = (I + A' (t + dt)) B' + v' (t + dt)
        //           = (P_t + A' dt) B' + (x_t + v' dt)
        //           = (I + (A' P_t^-1) dt) (P_t B') + (x_t + v' dt)
        //           = (I + (A' P_t^-1) dt) (P_t B' + x_t)
        //             - (I + (A' P_t^-1) dt) x_t + (x_t + v' dt)
        //           = (I + (A' P_t^-1) dt) (P_t B' + x_t)
        //             + (v' - (A' P_t^-1) x_t) dt                ......(*)
        //           = (I + (A' P_t^-1) dt) (P_t B' + x_t) + (P_t^-1 v') dt
        // (*): v' - (A' P_t^-1) x_t = v' - t A' (I + t A')^-1 v'
        //                           = (I - t A' (I + t A')^-1) v'
        //                           = ((I + t A')(I + t A')^-1
        //                             - t A' (I + t A')^-1) v'
        //                           = (I + t A')^-1 v'

        P0_t_inv = P0_t.inverse();
        P1_t_inv = P1_t.inverse();
        A0_t = *A0_ptr * P0_t_inv;
        A1_t = *A1_ptr * P1_t_inv;
        v0_t = P0_t_inv * *v0_ptr;
        v1_t = P1_t_inv * *v1_ptr;

        rho0 = spectral_norm(A0_t);
        rho1 = spectral_norm(A1_t);
        t_clip = (K - 1.0) / (K * rho0 + EPSILON_ABS);
        s0 = support_transformed(B0_ptr, P0_t * *P0_ptr, P0_t * *x0_ptr + x0_t,
                                 n);

        dt = solve_pointwise_toi_quad(
            B1_ptr, P1_t * *P1_ptr, P1_t * *x1_ptr + x1_t, A1_t - A0_t,
            v1_t - v0_t, s0, -n, K * rho0 * (rho0 + rho1),
            K * rho0 * (v1_t - v0_t).norm(), thickness);
        dt = std::min(dt, t_clip) * 0.9;

        assert(spectral_norm((Matrix3S::Identity() + dt * A0_t).inverse()) <=
               K);

        t = std::min(t + dt, (Scalar)1.0);

        if (t_history != nullptr) {
            t_history[k] = t;
        }

        if (t >= 1.0) {
            break;
        }
    }

    num_iter = k;

    if (k == max_iter) {
        spdlog::warn(
            "c5d_full_gjk_quad_pw did not converge after {} iterations",
            max_iter);
    }

    assert(t >= 0.0);

    d = gjk_transformed(B0, B1, P0 + t * A0, P1 + t * A1, x0 + t * v0,
                        x1 + t * v1, n);
    if (d < thickness * thickness * 0.9 - 1e-16) {
        spdlog::warn(
            "c5d_full_gjk_quad_pw: d = {0:.6e} < thickness^2 = {1:.6e}", d,
            thickness * thickness);
        spdlog::info("k = {0}, t = {1:.6e}", k, t);
        exit(1);
    }

    return t;
}

}  // namespace exp

}  // namespace c5d