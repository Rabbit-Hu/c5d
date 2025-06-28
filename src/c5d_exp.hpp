#pragma once

#include "gjk.hpp"

namespace c5d {

namespace exp {

Scalar c5d_full_gjk_linear(const Polytope *B0, const Polytope *B1,
                           const Matrix3S &P0, const Matrix3S &P1,
                           const Vector3S &x0, const Vector3S &x1,
                           const Matrix3S &A0, const Matrix3S &A1,
                           const Vector3S &v0, const Vector3S &v1,
                           Scalar slackness, Scalar thickness,
                           int max_iter, int &num_iter, Scalar *t_history);

Scalar c5d_full_gjk_quad(const Polytope *B0, const Polytope *B1,
                            const Matrix3S &P0, const Matrix3S &P1,
                            const Vector3S &x0, const Vector3S &x1,
                            const Matrix3S &A0, const Matrix3S &A1,
                            const Vector3S &v0, const Vector3S &v1,
                            Scalar slackness, Scalar thickness, Scalar K,
                            int max_iter, int &num_iter, Scalar *t_history);
                            
Scalar c5d_full_gjk_quad_pw(const Polytope *B0, const Polytope *B1,
                            const Matrix3S &P0, const Matrix3S &P1,
                            const Vector3S &x0, const Vector3S &x1,
                            const Matrix3S &A0, const Matrix3S &A1,
                            const Vector3S &v0, const Vector3S &v1,
                            Scalar slackness, Scalar thickness, Scalar K,
                            int max_iter, int &num_iter, Scalar *t_history);


}  // namespace exp

}  // namespace c5d