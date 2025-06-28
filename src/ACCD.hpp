#pragma once

#include <igl/edges.h>

#include <Eigen/Dense>
#include "gjk.hpp"

namespace c5d {

namespace exp {

using MatrixXS = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixXi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

Scalar additiveCCD(const MatrixXS &V0, const MatrixXS &V1, const MatrixXi &F0,
                   const MatrixXi &F1, const Matrix3S &P0, const Matrix3S &P1,
                   const Vector3S &x0, const Vector3S &x1, const Matrix3S &A0,
                   const Matrix3S &A1, const Vector3S &v0, const Vector3S &v1,
                   Scalar slackness, Scalar thickness, int max_iter,
                   double &num_iter);

// Scalar additiveCCD_parallel(const MatrixXS &V0, const MatrixXS &V1,
//                             const MatrixXi &F0, const MatrixXi &F1,
//                             const Matrix3S &P0, const Matrix3S &P1,
//                             const Vector3S &x0, const Vector3S &x1,
//                             const Matrix3S &A0, const Matrix3S &A1,
//                             const Vector3S &v0, const Vector3S &v1,
//                             Scalar slackness, Scalar thickness, int max_iter,
//                             double &num_iter);

}  // namespace exp

}  // namespace c5d
