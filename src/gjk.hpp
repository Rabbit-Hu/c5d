#pragma once

#include <Eigen/Dense>
#include <limits>

namespace c5d {

#ifdef C5D_DOUBLE_PRECISION
using Scalar = double;
#else
using Scalar = float;
#endif

const Scalar EPSILON = std::numeric_limits<Scalar>::epsilon();
const Scalar EPSILON_REL = EPSILON * 1e4;
const Scalar EPSILON_ABS = EPSILON * 1e2;
// const Scalar EPSILON_ABS_SQR = EPSILON_ABS * EPSILON_ABS;

using Vector3S = Eigen::Matrix<Scalar, 3, 1>;
using Matrix3S = Eigen::Matrix<Scalar, 3, 3>;
using Matrix34S = Eigen::Matrix<Scalar, 3, 4>;
using Matrix3XS = Eigen::Matrix<Scalar, 3, Eigen::Dynamic>;

class Polytope {
  public:
    virtual Vector3S support(const Vector3S &dir) const = 0;
    /// @brief Return the number of vertices.
    virtual inline int numVertices() const = 0;
    // /// @brief Return true if can call getVertex(i).
    // virtual inline bool hasVertex(int i) const = 0;
    /// @brief Get the i-th vertex of the polytope.
    virtual inline const Eigen::Ref<const Vector3S> getVertex(int i) const = 0;
    /// @brief For debugging: export the polytope as mesh.
    virtual void getMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) const = 0;
};

class PointSet : public Polytope {
  public:
    Matrix3XS vertices;
    Vector3S support(const Vector3S &dir) const override;
    inline int numVertices() const override { return vertices.cols(); }
    // inline bool hasVertex(int i) const override;
    inline const Eigen::Ref<const Vector3S> getVertex(int i) const override { return vertices.col(i); }
    void getMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) const override;
};

class Simplex {
  public:
    int n_vertices;
    Vector3S vertices[4];
};

Scalar gjk(Polytope *B0, Polytope *B1, Vector3S &v, int max_iter = 100);

Vector3S support_transformed(const Polytope *body, const Matrix3S &P,
                             const Vector3S &x, const Vector3S &v);

Scalar gjk_transformed(const Polytope *B0, const Polytope *B1,
                       const Matrix3S &A0, const Matrix3S &A1,
                       const Vector3S &x0, const Vector3S &x1, Vector3S &v,
                       int max_iter = 100);

}  // namespace c5d
