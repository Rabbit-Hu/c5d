#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <random>

#include "ACCD.hpp"
#include "DataConfig.hpp"
#include "DataPoint.hpp"
#include "c5d_exp.hpp"

namespace c5d {

namespace exp {

static void verticesToConvexHull(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

    std::vector<Point_3> points;
    for (int i = 0; i < V.rows(); i++) {
        points.push_back(Point_3(V(i, 0), V(i, 1), V(i, 2)));
    }

    Surface_mesh sm;
    CGAL::convex_hull_3(points.begin(), points.end(), sm);

    // Convert Surface_mesh to Eigen matrices V and F
    V.resize(sm.number_of_vertices(), 3);
    F.resize(sm.number_of_faces(), 3);

    int cnt_V = 0;
    for (auto vi : sm.vertices()) {
        Point_3 p = sm.point(vi);
        V.row(cnt_V++) = Eigen::RowVector3d(p.x(), p.y(), p.z());
    }
    int cnt_F = 0;
    for (auto fi : sm.faces()) {
        CGAL::Vertex_around_face_circulator<Surface_mesh> vc(sm.halfedge(fi),
                                                             sm),
            done(vc);
        int j = 0;
        do {
            F(cnt_F, j++) = *vc++;
        } while (vc != done);
        cnt_F++;
    }
}

static Eigen::Matrix3d generateRotationMatrix(std::mt19937 &rng) {
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    Eigen::Matrix3d R;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R(i, j) = normal_dist(rng);
        }
    }
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    return svd.matrixU() * svd.matrixV().transpose();
}

static Eigen::MatrixXd generateGaussianMatrix(std::mt19937 &rng, int rows,
                                              int cols) {
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    Eigen::MatrixXd M(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            M(i, j) = normal_dist(rng);
        }
    }
    return M;
}

static Eigen::Matrix3d skewSymmetricMatrix(const Eigen::Vector3d &v) {
    Eigen::Matrix3d S;
    S << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    return S;
}

static Eigen::Vector3d sampleUnitBall(std::mt19937 &rng) {
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    // sample direction
    Eigen::Vector3d x(normal_dist(rng), normal_dist(rng), normal_dist(rng));
    // sample norm
    double r = std::cbrt(uniform_dist(rng));
    return r * x / x.norm();
}

static double spectralNorm(const Eigen::Matrix3d &M) {
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    return svd.singularValues()(0);
}

static double bruteForceCCD(
    const PointSet &B0, const PointSet &B1, const Eigen::Matrix3d &P0,
    const Eigen::Matrix3d &P1, const Eigen::Vector3d &x0,
    const Eigen::Vector3d &x1, const Eigen::Matrix3d &A0,
    const Eigen::Matrix3d &A1, const Eigen::Vector3d &v0,
    const Eigen::Vector3d &v1, double step_size) {
    for (double t = 0.0; t <= 1.0; t += step_size) {
        Eigen::Vector3d normal;
        double distance =
            gjk_transformed(&B0, &B1, P0 + t * A0, P1 + t * A1, x0 + t * v0,
                            x1 + t * v1, normal, 100);
        if (distance <= 0.0) {
            return t - step_size;
        }
    }
    return 1e9;
}

static DataPoint generateDataPoint(DataConfig config, std::mt19937 &rng) {
    DataPoint dp;

    // Gaussian distribution
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    // Generate initial poses: P0, P1, x0, x1
    int pose_rejection_count = 0;
    int velocity_rejection_count = 0;
    while (true) {
        // Generate raw vertices: V0, V1
        int n_raw_vertices = config.vertex.n_raw_vertices;
        dp.V0 = generateGaussianMatrix(rng, n_raw_vertices, 3) *
                config.vertex.vertex_std;
        dp.V1 = generateGaussianMatrix(rng, n_raw_vertices, 3) *
                config.vertex.vertex_std;

        // Convert raw vertices to convex hulls: V0, V1, F0, F1
        verticesToConvexHull(dp.V0, dp.F0);
        verticesToConvexHull(dp.V1, dp.F1);

        PointSet B0, B1;
        B0.vertices = dp.V0.transpose().cast<Scalar>();
        B1.vertices = dp.V1.transpose().cast<Scalar>();

        dp.P0 =
            generateRotationMatrix(rng) +
            generateGaussianMatrix(rng, 3, 3) * config.pose.rotation_noise_std;
        dp.P1 =
            generateRotationMatrix(rng) +
            generateGaussianMatrix(rng, 3, 3) * config.pose.rotation_noise_std;
        dp.x0 = sampleUnitBall(rng) * config.pose.translation_max_norm;
        dp.x1 = sampleUnitBall(rng) * config.pose.translation_max_norm;

        // Reject if initially colliding
        Vector3S init_normal;
        double init_distance =
            gjk_transformed(&B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, init_normal);
        if (init_distance <= 1e-9) {
            pose_rejection_count++;
            spdlog::debug("Rejected due to initially colliding");
            continue;
        }

        if (config.velocity.near_skew) {
            dp.A0 = skewSymmetricMatrix(sampleUnitBall(rng));
            dp.A1 = skewSymmetricMatrix(sampleUnitBall(rng));
        } else {
            dp.A0 = generateGaussianMatrix(rng, 3, 3);
            dp.A1 = generateGaussianMatrix(rng, 3, 3);
            dp.A0 = dp.A0 / spectralNorm(dp.A0);
            dp.A1 = dp.A1 / spectralNorm(dp.A1);
        }
        
        dp.A0 = dp.A0 * config.velocity.affine_max_norm * uniform_dist(rng) +
                generateGaussianMatrix(rng, 3, 3) *
                    config.velocity.affine_noise_std;
        dp.A1 = dp.A1 * config.velocity.affine_max_norm * uniform_dist(rng) +
                generateGaussianMatrix(rng, 3, 3) *
                    config.velocity.affine_noise_std;
        dp.v0 = sampleUnitBall(rng) * config.velocity.linear_max_norm;
        dp.v1 = sampleUnitBall(rng) * config.velocity.linear_max_norm;

        Scalar slackness = 1.0 - 1e-6;
        Scalar thickness = 1e-9;

        int max_iter = 10000;
        int num_iter = 0;
        double num_iter_d = 0.0;
        Scalar *t_history = new Scalar[max_iter];

        spdlog::debug("Running c5d_quad_pw:");
        double quad_pw_toi = c5d_full_gjk_quad_pw(
            &B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0, dp.A1, dp.v0, dp.v1,
            slackness, thickness, 1.1, max_iter, num_iter, t_history);

        spdlog::debug("Running c5d_quad:");
        double quad_toi = c5d_full_gjk_quad(
            &B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0, dp.A1, dp.v0, dp.v1,
            slackness, thickness, 1.1, max_iter, num_iter, t_history);

        spdlog::debug("Running c5d_linear:");
        double linear_toi = c5d_full_gjk_linear(
            &B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0, dp.A1, dp.v0, dp.v1,
            slackness, thickness, max_iter, num_iter, t_history);

        spdlog::debug("Running ACCD:");
        double accd_toi = additiveCCD(
            dp.V0, dp.V1, dp.F0, dp.F1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0,
            dp.A1, dp.v0, dp.v1, slackness, thickness, max_iter, num_iter_d);

        dp.toi = accd_toi;

        delete[] t_history;

        spdlog::info("PW ToI: {}, Quad ToI: {}, Linear ToI: {}, ACCD ToI: {}",
                     quad_pw_toi, quad_toi, linear_toi, accd_toi);
        if (abs(quad_pw_toi - accd_toi) > 1e-5 ||
            abs(quad_toi - accd_toi) > 1e-5 ||
            abs(linear_toi - accd_toi) > 1e-5) {
            double brute_toi = bruteForceCCD(B0, B1, dp.P0, dp.P1, dp.x0, dp.x1,
                                             dp.A0, dp.A1, dp.v0, dp.v1, 5e-6);
            spdlog::error(
                "PW ToI: {}, Linear ToI: {}, ACCD ToI: {}, Brute ToI: {}",
                quad_pw_toi, linear_toi, accd_toi, brute_toi);
            // exit(1);
            continue;
        }

        if (dp.toi >= 1.0 || dp.toi < 0.01) {
            velocity_rejection_count++;
            spdlog::debug(
                "Rejected due to no collision or too early collision (toi = "
                "{})",
                dp.toi);
            continue;
        }

        // double brute_toi = bruteForceCCD(B0, B1, dp.P0, dp.P1, dp.x0, dp.x1,
        //                                  dp.A0, dp.A1, dp.v0, dp.v1, 1e-5);

        // spdlog::debug("C5D-PW ToI: {}, ACCD ToI: {}, Brute ToI: {}", dp.toi,
        //               accd_toi, brute_toi);

        break;
    }
    spdlog::debug("ToI: {}", dp.toi);
    spdlog::info("Pose rejection count: {}, velocity rejection count: {}",
                  pose_rejection_count, velocity_rejection_count);

    return dp;
}

struct Dataset {
    std::vector<DataPoint> datapoints;

    void generate_from_config(DataConfig config, bool save,
                              std::string output_root) {
        std::mt19937 rng(config.seed);
        for (int i = 0; i < config.n_data; i++) {
            datapoints.push_back(generateDataPoint(config, rng));
            spdlog::info("Generated data point {}/{}", i + 1, config.n_data);

            if (save) {
                std::string dir = output_root + "/" + std::to_string(i);
                std::filesystem::create_directories(dir);
                datapoints[i].save(dir);
            }
        }
    }

    void load(std::string root) {
        datapoints.clear();
        // walk through all directories in root
        // for (const auto &entry : std::filesystem::directory_iterator(root)) {
        //     if (entry.is_directory()) {
        //         DataPoint dp;
        //         dp.load(entry.path());
        //         datapoints.push_back(dp);
        //     }
        // }
        int i = 0;
        while (true) {
            std::string dir = root + "/" + std::to_string(i);
            if (!std::filesystem::exists(dir)) {
                break;
            }
            DataPoint dp;
            dp.load(dir);
            datapoints.push_back(dp);
            i++;
        }
        spdlog::info("Loaded {} data points", datapoints.size());
    }

    // void save(std::string root) {
    //     for (int i = 0; i < datapoints.size(); i++) {
    //         std::string dir = root + "/" + std::to_string(i);
    //         std::filesystem::create_directories(dir);
    //         datapoints[i].save(dir);
    //     }
    // }
};

}  // namespace exp

}  // namespace c5d
