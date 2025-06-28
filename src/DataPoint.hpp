#pragma once

#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace c5d {

namespace exp {

struct DataPoint {
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi F0, F1;
    Eigen::Matrix3d P0, P1, A0, A1;
    Eigen::Vector3d x0, x1, v0, v1;
    double toi;

    void save(std::string root) {
        // Save meshes (V0, F0) and (V1, F1) to binary .ply files
        // root/mesh_0.ply, root/mesh_1.ply
        igl::writePLY(root + "/mesh_0.ply", V0, F0, igl::FileEncoding::Binary);
        igl::writePLY(root + "/mesh_1.ply", V1, F1, igl::FileEncoding::Binary);

        std::ofstream file(root + "/data.bin",
                           std::ios::out | std::ios::binary);

        // Save P0, x0, P1, x1 to binary file root/data.bin
        file.write(reinterpret_cast<const char *>(P0.data()),
                   9 * sizeof(double));
        file.write(reinterpret_cast<const char *>(P1.data()),
                   9 * sizeof(double));
        file.write(reinterpret_cast<const char *>(x0.data()),
                   3 * sizeof(double));
        file.write(reinterpret_cast<const char *>(x1.data()),
                   3 * sizeof(double));

        // Save A0, A1, v0, v1
        file.write(reinterpret_cast<const char *>(A0.data()),
                   9 * sizeof(double));
        file.write(reinterpret_cast<const char *>(A1.data()),
                   9 * sizeof(double));
        file.write(reinterpret_cast<const char *>(v0.data()),
                   3 * sizeof(double));
        file.write(reinterpret_cast<const char *>(v1.data()),
                   3 * sizeof(double));

        // Save toi
        file.write(reinterpret_cast<const char *>(&toi), sizeof(double));

        file.close();
    }

    void load(std::string root) {
        // Load meshes (V0, F0) and (V1, F1) from binary .ply files
        // root/mesh_0.ply, root/mesh_1.ply
        igl::readPLY(root + "/mesh_0.ply", V0, F0);
        igl::readPLY(root + "/mesh_1.ply", V1, F1);

        std::ifstream file(root + "/data.bin", std::ios::in | std::ios::binary);
        // Load P0, x0, P1, x1 from binary file root/data.bin
        file.read(reinterpret_cast<char *>(P0.data()), 9 * sizeof(double));
        file.read(reinterpret_cast<char *>(P1.data()), 9 * sizeof(double));
        file.read(reinterpret_cast<char *>(x0.data()), 3 * sizeof(double));
        file.read(reinterpret_cast<char *>(x1.data()), 3 * sizeof(double));

        // Load A0, A1, v0, v1
        file.read(reinterpret_cast<char *>(A0.data()), 9 * sizeof(double));
        file.read(reinterpret_cast<char *>(A1.data()), 9 * sizeof(double));
        file.read(reinterpret_cast<char *>(v0.data()), 3 * sizeof(double));
        file.read(reinterpret_cast<char *>(v1.data()), 3 * sizeof(double));

        // Load toi
        file.read(reinterpret_cast<char *>(&toi), sizeof(double));

        file.close();
    }
};

}  // namespace exp

}  // namespace c5d
