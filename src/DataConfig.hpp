#pragma once

#include <iostream>

#include "DataPoint.hpp"

namespace c5d {

namespace exp {

struct DataConfig {
    int seed, n_data;
    struct {
        int n_raw_vertices;
        double vertex_std;

        void loadYaml(const YAML::Node& config) {
            n_raw_vertices = config["n_raw_vertices"].as<int>();
            vertex_std = config["vertex_std"].as<double>();
        }
    } vertex;
    struct {
        double rotation_noise_std, translation_max_norm;

        void loadYaml(const YAML::Node& config) {
            rotation_noise_std = config["rotation_noise_std"].as<double>();
            translation_max_norm = config["translation_max_norm"].as<double>();
        }
    } pose;
    struct {
        bool near_skew;
        double affine_max_norm, affine_noise_std, linear_max_norm;

        void loadYaml(const YAML::Node& config) {
            near_skew = config["near_skew"].as<bool>();
            affine_max_norm = config["affine_max_norm"].as<double>();
            affine_noise_std = config["affine_noise_std"].as<double>();
            linear_max_norm = config["linear_max_norm"].as<double>();
        }
    } velocity;

    void loadYaml(const YAML::Node& config) {
        seed = config["seed"].as<int>();
        n_data = config["n_data"].as<int>();

        vertex.loadYaml(config["vertex"]);
        pose.loadYaml(config["pose"]);
        velocity.loadYaml(config["velocity"]);
    }
};

std::ostream& operator<<(std::ostream& os, const DataConfig& dc) {
    os << "DataConfig:\
    \n  seed: "
       << dc.seed << "\
    \n  n_data: "
       << dc.n_data << "\
    \n  vertex:\
    \n    n_raw_vertices: "
       << dc.vertex.n_raw_vertices << "\
    \n    vertex_std: "
       << dc.vertex.vertex_std << "\
    \n  pose:\
    \n    rotation_noise_std: "
       << dc.pose.rotation_noise_std << "\
    \n    translation_max_norm: "
       << dc.pose.translation_max_norm << "\
    \n  velocity:\
    \n    near_skew: "
       << dc.velocity.near_skew << "\
    \n    affine_max_norm: "
       << dc.velocity.affine_max_norm << "\
    \n    affine_noise_std: "
       << dc.velocity.affine_noise_std << "\
    \n    linear_max_norm: "
       << dc.velocity.linear_max_norm << std::endl;
    return os;
}

}  // namespace exp

}  // namespace c5d