#include "gjk.hpp"

#include <spdlog/spdlog.h>

namespace c5d {

static inline Scalar determinant(const Vector3S &p, const Vector3S &q,
                                 const Vector3S &r) {
    return p.dot(q.cross(r));
}

/// @brief Project origin on the line (p, q).
static inline Vector3S projectOnLine(const Vector3S &p, const Vector3S &q) {
    Vector3S pq = q - p;
    Scalar tmp = p.dot(pq) / pq.dot(pq);
    return p - pq * tmp;
}

/// @brief Project v on the plane (p, q, r).
static inline Vector3S projectOnPlane(const Vector3S &p, const Vector3S &q,
                                      const Vector3S &r) {
    Vector3S pq = q - p;
    Vector3S pr = r - p;
    Vector3S n = pq.cross(pr);
    Scalar tmp = p.dot(n) / n.dot(n);
    return n * tmp;
}

/// @brief Determine the projection of the origin is on which side of p on the
/// line (p, q). Return true if origin is on the same side as q.
static inline bool getSide1D(const Vector3S &p, const Vector3S &q) {
    return p.dot(p - q) > 0;  // keep q if true
}

/// @brief Determine the projection of the origin is on which side of line (p,
/// q) on the plane (p, q, r). Return true if origin is on the same side as r.
static inline bool getSide2D(const Vector3S &p, const Vector3S &q,
                             const Vector3S &r) {
    Vector3S pq = q - p;
    Vector3S pr = r - p;

    return p.dot(pq.cross(pq.cross(pr))) >= 0;  // keep r if true
}

/// @brief Determine the origin is on which side of the plane (p, q, r).
/// Return true if origin is on the same side as s.
static inline bool getSide3D(const Vector3S &p, const Vector3S &q,
                             const Vector3S &r, const Vector3S &s) {
    Vector3S pq = q - p;
    Vector3S pr = r - p;
    Vector3S ps = s - p;

    Vector3S n = pq.cross(pr);
    return (p.dot(n) * ps.dot(n)) < 0;  // keep s if true
}

static inline void updateSimplex(Simplex *s, const Vector3S &s0, Vector3S &v) {
    s->n_vertices = 1;
    s->vertices[0] = s0;
    v = s0;
}

static inline void updateSimplex(Simplex *s, const Vector3S &s0,
                                 const Vector3S &s1, Vector3S &v) {
    s->n_vertices = 2;
    s->vertices[0] = s0;
    s->vertices[1] = s1;
    v = projectOnLine(s1, s0);
}

static inline void updateSimplex(Simplex *s, const Vector3S &s0,
                                 const Vector3S &s1, const Vector3S &s2,
                                 Vector3S &v) {
    s->n_vertices = 3;
    s->vertices[0] = s0;
    s->vertices[1] = s1;
    s->vertices[2] = s2;
    v = projectOnPlane(s2, s1, s0);
}

static inline int S1D(Simplex *s, Vector3S &v) {
    Vector3S s0 = s->vertices[0];
    Vector3S s1 = s->vertices[1];  // always keep s1

    if (getSide1D(s1, s0)) {  // keep s0
        // spdlog::debug("keep s0");
        v = projectOnLine(s1, s0);
    } else {  // discard s0
        // spdlog::debug("discard s0");
        updateSimplex(s, s1, v);
    }

    return 0;
}

static inline int S2D(Simplex *s, Vector3S &v) {
    Vector3S s0 = s->vertices[0];
    Vector3S s1 = s->vertices[1];
    Vector3S s2 = s->vertices[2];

    bool side1d_20 = getSide1D(s2, s0);
    bool side1d_21 = getSide1D(s2, s1);

    // V2, E20, E21, or T active
    if (side1d_21) {  // E20, E21, or T active
        bool side2d_210 = getSide2D(s2, s1, s0);
        if (side2d_210) {  // E20 or T active
            bool side2d_201 = getSide2D(s2, s0, s1);
            if (side2d_201) {  // T active
                v = projectOnPlane(s2, s1, s0);
            } else {  // E20 active
                v = projectOnLine(s2, s0);
                s->n_vertices = 2;
                s->vertices[1] = s2;
            }
        } else {  // E21 active
            v = projectOnLine(s2, s1);
            s->n_vertices = 2;
            s->vertices[0] = s1;
            s->vertices[1] = s2;
        }
    } else {              // V2, E20, or T active
        if (side1d_20) {  // E20 or T active
            bool side2d_201 = getSide2D(s2, s0, s1);
            if (side2d_201) {  // T active
                v = projectOnPlane(s2, s1, s0);
            } else {  // E20 active
                v = projectOnLine(s2, s0);
                s->n_vertices = 2;
                s->vertices[1] = s2;
            }
        } else {  // V2 active
            v = s2;
            s->n_vertices = 1;
            s->vertices[0] = s2;
        }
    }

    return 0;
}

static inline int S3D(Simplex *s, Vector3S &v) {
    Vector3S s0 = s->vertices[0];
    Vector3S s1 = s->vertices[1];
    Vector3S s2 = s->vertices[2];
    Vector3S s3 = s->vertices[3];

    // Scalar det = determinant(s1 - s0, s2 - s0, s3 - s0);
    // if (det * det <= EPSILON_ABS) {
    //     spdlog::warn("degenerate tetrahedron");
    //     v = s3;  // degenrate tetrahedron
    //     s->n_vertices = 1;
    //     s->vertices[0] = s3;
    //     return;
    // }

    int side1d[3] = {getSide1D(s3, s0), getSide1D(s3, s1), getSide1D(s3, s2)};
    int sum_side1d = side1d[0] + side1d[1] + side1d[2];

    if (sum_side1d == 0) {
        updateSimplex(s, s3, v);
        return 0;
    }

    int side3d[3] = {getSide3D(s3, s1, s2, s0), getSide3D(s3, s2, s0, s1),
                     getSide3D(s3, s0, s1, s2)};
    int sum_side3d = side3d[0] + side3d[1] + side3d[2];

    int i, j, k;
    Vector3S si, sj, sk;

    if (sum_side3d == 3) {
        v.setZero();  // origin is inside tetrahedron
        return 0;
    } else if (sum_side3d == 2) {
        if (!side3d[2]) {
            // discard s2
        } else if (!side3d[1]) {
            // discard s1
            s->vertices[1] = s2;
        } else {
            // discard s0
            s->vertices[0] = s1;  // TODO: necessary?
            s->vertices[1] = s2;
        }
        s->n_vertices = 3;
        s->vertices[2] = s3;
        return S2D(s, v);
    } else if (sum_side3d == 1) {
        // Consider two triangles! Reorder so that the shared edge is k
        if (side3d[0]) {
            // Consider triangles 320, 301
            k = 0, i = 2, j = 1;
        } else if (side3d[1]) {
            // Consider triangles 301, 312
            k = 1, i = 0, j = 2;
        } else {
            // Consider triangles 312, 320
            k = 2, i = 1, j = 0;
        }
        si = s->vertices[i];
        sj = s->vertices[j];
        sk = s->vertices[k];

        if (sum_side1d == 1) {
            if (side1d[k]) {  // k active; E3k or T3ik or T3jk active
                if (getSide2D(s3, sk, si)) {  // T3ik active
                    updateSimplex(s, sk, si, s3, v);
                } else if (getSide2D(s3, sk, sj)) {  // T3jk active
                    updateSimplex(s, sk, sj, s3, v);
                } else {  // E3k active
                    updateSimplex(s, sk, s3, v);
                }
            } else if (side1d[i]) {           // i active; E3i or T3ik active
                if (getSide2D(s3, si, sk)) {  // T3ik active
                    updateSimplex(s, si, sk, s3, v);
                } else {  // E3i active
                    updateSimplex(s, si, s3, v);
                }
            } else {                          // j active; E3j or T3jk active
                if (getSide2D(s3, sj, sk)) {  // T3jk active
                    updateSimplex(s, sj, sk, s3, v);
                } else {  // E3j active
                    updateSimplex(s, sj, s3, v);
                }
            }
        } else if (sum_side1d == 2) {
            if (side1d[i]) {  // should be i and k
                if (getSide2D(s3, sk, si)) {
                    if (getSide2D(s3, si, sk)) {  // T3ik active
                        updateSimplex(s, si, sk, s3, v);
                    } else {  // ~~E3i active~~ E3k active
                        // updateSimplex(s, si, s3, v);
                        updateSimplex(s, sk, s3, v);  // TODO: check this
                    }
                } else {
                    if (getSide2D(s3, sk, sj)) {  // T3jk active
                        updateSimplex(s, sk, sj, s3, v);
                    } else {  // E3k active
                        updateSimplex(s, sk, s3, v);
                    }
                }
            } else if (side1d[j]) {  // should be j and k
                if (getSide2D(s3, sk, sj)) {
                    if (getSide2D(s3, sj, sk)) {  // T3jk active
                        updateSimplex(s, sj, sk, s3, v);
                    } else {  // E3j active
                        updateSimplex(s, sj, s3, v);
                    }
                } else {
                    if (getSide2D(s3, sk, si)) {  // T3ik active
                        updateSimplex(s, sk, si, s3, v);
                    } else {  // E3k active
                        updateSimplex(s, sk, s3, v);
                    }
                }
            } else {
                spdlog::error("unhandled case");
                return 1;
            }
        } else {
            bool side2d_ik = getSide2D(s3, si, sk);
            bool side2d_jk = getSide2D(s3, sj, sk);
            bool side2d_ki = getSide2D(s3, sk, si);
            bool side2d_kj = getSide2D(s3, sk, sj);
            if (side2d_ki && side2d_kj) {
                spdlog::error(
                    "unhandled case: sum_side1d = {}, side2d_ki = 0, side2d_kj "
                    "= 0",
                    sum_side1d);
                return 1;
            } else if (!side2d_ki && !side2d_kj) {
                updateSimplex(s, sk, s3, v);
            } else if (side2d_kj) {  // discard i
                if (side2d_jk) {     // keep k
                    updateSimplex(s, sj, sk, s3, v);
                } else {
                    updateSimplex(s, sj, s3, v);
                }
            } else {              // discard j
                if (side2d_ik) {  // keep k
                    updateSimplex(s, si, sk, s3, v);
                } else {
                    updateSimplex(s, si, s3, v);
                }
            }
        }
        return 0;
    } else {  // sum_side3d == 0
        if (sum_side1d == 1) {
            // select i s.t. side1d[i] == 1
            if (side1d[0]) {
                k = 1, i = 0, j = 2;
            } else if (side1d[1]) {
                k = 2, i = 1, j = 0;
            } else {
                k = 0, i = 2, j = 1;
            }
            si = s->vertices[i];
            sj = s->vertices[j];
            sk = s->vertices[k];

            if (getSide2D(s3, si, sj)) {  // T3ij active
                updateSimplex(s, si, sj, s3, v);
            } else if (getSide2D(s3, si, sk)) {  // T3ik active
                updateSimplex(s, si, sk, s3, v);
            } else {  // E3i active
                updateSimplex(s, si, s3, v);
            }
        } else if (sum_side1d == 2) {
            // select i s.t. side1d[i] == 0
            if (!side1d[0]) {
                k = 1, i = 0, j = 2;
            } else if (!side1d[1]) {
                k = 2, i = 1, j = 0;
            } else {
                k = 0, i = 2, j = 1;
            }
            si = s->vertices[i];
            sj = s->vertices[j];
            sk = s->vertices[k];

            if (getSide2D(s3, sj, sk)) {
                if (getSide2D(s3, sk, sj)) {  // T3jk active
                    updateSimplex(s, sj, sk, s3, v);
                } else if (getSide2D(s3, sk, si)) {  // T3ik active
                    updateSimplex(s, sk, si, s3, v);
                } else {  // E3k active
                    updateSimplex(s, sk, s3, v);
                }
            } else if (getSide2D(s3, sj, si)) {  // T3ij active
                updateSimplex(s, sj, si, s3, v);
            } else {  // E3j active
                updateSimplex(s, sj, s3, v);
            }
        } else {
            spdlog::error("unhandled case: sum_side3d = {}, sum_side1d = {}",
                          sum_side3d, sum_side1d);
            return 1;
        }
        return 0;
    }
    spdlog::error("unhandled end of S3D");
}

static inline int subAlgorithm(Simplex *s, Vector3S &v) {
    if (s->n_vertices == 2) {
        return S1D(s, v);
    } else if (s->n_vertices == 3) {
        return S2D(s, v);
    } else if (s->n_vertices == 4) {
        return S3D(s, v);
    } else {
        spdlog::error("unhandled case: n_vertices = {}", s->n_vertices);
        return 1;
    }
}

Scalar gjk(Polytope *B0, Polytope *B1, Vector3S &v, int max_iter) {
    Simplex s;
    v = B0->support(Vector3S::Zero()) - B1->support(Vector3S::Zero());
    s.n_vertices = 1;
    s.vertices[0] = v;

    for (int i = 0; i < max_iter; i++) {
        Vector3S w = B0->support(-v) - B1->support(v);

        Scalar vw = v.dot(w);
        Scalar v_norm_sqr = v.squaredNorm();
        Scalar dist_w_simplex = v_norm_sqr - vw;

        if (dist_w_simplex <= EPSILON_REL * v_norm_sqr ||
            dist_w_simplex < EPSILON_ABS) {
            break;  // w on the simplex, converged
        }

        s.vertices[s.n_vertices++] = w;
        int sub_algo_ret = subAlgorithm(&s, v);
        if (sub_algo_ret) {
            spdlog::error("subAlgorithm returned {}", sub_algo_ret);
            return std::numeric_limits<Scalar>::infinity();
        }

        if (s.n_vertices == 4 || v.squaredNorm() < EPSILON_ABS) {
            break;  // origin is inside (B0 - B1)
        }

        if (i == max_iter - 1) {
            spdlog::warn("GJK did not converge after {} iterations", max_iter);
        }
    }

    return v.squaredNorm();
}

Vector3S support_transformed(const Polytope *body, const Matrix3S &P,
                             const Vector3S &x, const Vector3S &v) {
    return P * body->support(P.transpose() * v) + x;
}

Scalar gjk_transformed(const Polytope *B0, const Polytope *B1,
                       const Matrix3S &P0, const Matrix3S &P1,
                       const Vector3S &x0, const Vector3S &x1, Vector3S &v,
                       int max_iter) {
    Simplex s;

    v = support_transformed(B0, P0, x0, Vector3S::Zero()) -
        support_transformed(B1, P1, x1, Vector3S::Zero());
    s.n_vertices = 1;
    s.vertices[0] = v;

    for (int i = 0; i < max_iter; i++) {
        Vector3S w = support_transformed(B0, P0, x0, -v) -
                     support_transformed(B1, P1, x1, v);

        Scalar vw = v.dot(w);
        Scalar v_norm_sqr = v.squaredNorm();
        Scalar dist_w_simplex = v_norm_sqr - vw;

        if (dist_w_simplex < EPSILON_REL * v_norm_sqr ||
            dist_w_simplex < EPSILON_ABS) {
            break;  // w on the simplex, converged
        }

        s.vertices[s.n_vertices++] = w;
        subAlgorithm(&s, v);

        if (s.n_vertices == 4 || v.squaredNorm() < EPSILON_ABS) {
            break;  // origin is inside (B0 - B1)
        }

        if (i == max_iter - 1) {
            spdlog::warn("GJK did not converge after {} iterations", max_iter);
        }
    }

    return v.squaredNorm();
}

Vector3S PointSet::support(const Vector3S &dir) const {
    assert(vertices.cols() > 0);

    Scalar max_dot = -std::numeric_limits<Scalar>::infinity();
    int max_vertex_id = -1;
    for (int i = 0; i < vertices.cols(); i++) {
        Scalar dot = vertices.col(i).dot(dir);
        if (dot > max_dot) {
            max_dot = dot;
            max_vertex_id = i;
        }
    }
    assert(max_vertex_id >= 0);
    return vertices.col(max_vertex_id);
}

inline int PointSet::numVertices() const { return vertices.cols(); }

// inline bool PointSet::hasVertex(int i) const {
//     return i < vertices.cols() && i >= 0;
// }

inline const Eigen::Ref<const Vector3S> PointSet::getVertex(int i) const {
    return vertices.col(i);
}

void PointSet::getMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F) const {
    V.resize(3, vertices.cols());
    for (int i = 0; i < vertices.cols(); i++) {
        V.col(i) = vertices.col(i);
    }
    F.resize(3, 0);
}

}  // namespace c5d
