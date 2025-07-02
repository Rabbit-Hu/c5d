#include <boost/timer/timer.hpp>
#include <iostream>
#include <fstream>

#include "Dataset.hpp"

#ifdef C5D_DOUBLE_PRECISION
using Scalar = double;
#else
using Scalar = float;
#endif

int main(int argc, char *argv[]) {
    spdlog::set_level(spdlog::level::info);

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <dataset_dir> <output_dir>"
                  << std::endl;
        return 1;
    }

    spdlog::info("Loading dataset from {}", argv[1]);
    c5d::exp::Dataset dataset;
    dataset.load(argv[1]);

    Scalar slackness = 1.0 - 1e-6;
    Scalar thickness = 1e-9;
    Scalar K = 1.1;
    constexpr int n_key_iters = 3;
    constexpr int max_iter = 10000;
    Scalar key_ratios[n_key_iters] = {0.9, 0.95, 0.99};
    int key_iters[n_key_iters];
    int n_not_converged_per_iter[4][n_key_iters][max_iter + 1];
    memset(n_not_converged_per_iter, 0, sizeof(n_not_converged_per_iter));
    Scalar *t_history = new Scalar[max_iter + 1];

    // Save into one csv file for each method
    // headers: id, toi, num_iter, key_iter_0.9, key_iter_0.95, key_iter_0.99,
    // wall_time

    std::string header_str =
        "id,toi,num_iter,key_iter_0.9,key_iter_0.95,"
        "key_iter_0.99,wall_time\n";
    std::string output_dir(argv[2]);
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
    }
    std::ofstream c5d_linear_file(output_dir + "/c5d_linear.csv");
    std::ofstream c5d_quad_file(output_dir + "/c5d_quad.csv");
    std::ofstream c5d_quad_pw_file(output_dir + "/c5d_quad_pw.csv");
    std::ofstream accd_file(output_dir + "/accd.csv");
    // std::fstream accd_parallel_file(output_dir + "/accd_parallel.csv",
    //                                 std::ios::out | std::ios::trunc);
    c5d_linear_file << header_str;
    c5d_quad_file << header_str;
    c5d_quad_pw_file << header_str;
    accd_file << header_str;
    // accd_parallel_file << header_str;

    Scalar accd_toi, accd_parallel_toi;

    for (int i = 0; i < dataset.datapoints.size(); i++) {
        // for (int i = 0; i < 10; i++) {
        std::cout << "Data point " << i << std::endl;

        const c5d::exp::DataPoint &dp = dataset.datapoints[i];

        c5d::PointSet B0, B1;
        B0.vertices = dp.V0.transpose().cast<Scalar>();
        B1.vertices = dp.V1.transpose().cast<Scalar>();

        // Test C5D-Linear
        {
            boost::timer::auto_cpu_timer timer;
            int num_iter;
            Scalar toi = c5d::exp::c5d_full_gjk_linear(
                &B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0, dp.A1, dp.v0,
                dp.v1, slackness, thickness, max_iter, num_iter, t_history);
            // std::cout << "C5D-Linear toi: " << toi << std::endl;
            timer.stop();

            // Compute num_iters to reach key ratios
            for (int k = 0; k < n_key_iters; k++) {
                key_iters[k] = std::lower_bound(t_history, t_history + num_iter,
                                                key_ratios[k] * toi) -
                               t_history;
                for (int iter = 0; iter < key_iters[k]; iter++) {
                    n_not_converged_per_iter[0][k][iter]++;
                }
            }
            std::string key_iters_str = "";
            for (int k = 0; k < n_key_iters; k++) {
                key_iters_str += std::to_string(key_iters[k]);
                if (k < n_key_iters - 1) {
                    key_iters_str += ",";
                }
            }

            std::cout << "    C5D-Linear num_iter: " << num_iter
                      << "; key_iters: " << key_iters_str
                      << "; wall time: " << timer.elapsed().wall / 1e6 << " ms"
                      << "; toi: " << toi << std::endl;

            c5d_linear_file << i << "," << toi << "," << num_iter << ","
                            << key_iters_str << ","
                            << timer.elapsed().wall / 1e6 << "\n";
        }

        // Test C5D-Quad
        {
            boost::timer::auto_cpu_timer timer;
            int num_iter;
            Scalar toi = c5d::exp::c5d_full_gjk_quad(
                &B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0, dp.A1, dp.v0,
                dp.v1, slackness, thickness, K, max_iter, num_iter, t_history);
            // std::cout << "C5D-Quad toi: " << toi << std::endl;
            timer.stop();

            // Compute num_iters to reach key ratios
            for (int k = 0; k < n_key_iters; k++) {
                key_iters[k] = std::lower_bound(t_history, t_history + num_iter,
                                                key_ratios[k] * toi) -
                               t_history;
                for (int iter = 0; iter < key_iters[k]; iter++) {
                    n_not_converged_per_iter[1][k][iter]++;
                }
            }
            std::string key_iters_str = "";
            for (int k = 0; k < n_key_iters; k++) {
                key_iters_str += std::to_string(key_iters[k]);
                if (k < n_key_iters - 1) {
                    key_iters_str += ",";
                }
            }

            std::cout << "    C5D-Quad num_iter: " << num_iter
                      << "; key_iters: " << key_iters_str
                      << "; wall time: " << timer.elapsed().wall / 1e6 << " ms"
                      << "; toi: " << toi << std::endl;

            c5d_quad_file << i << "," << toi << "," << num_iter << ","
                          << key_iters_str << "," << timer.elapsed().wall / 1e6
                          << "\n";
        }

        // Test C5D-Quad-PW
        {
            boost::timer::auto_cpu_timer timer;
            int num_iter;
            Scalar toi = c5d::exp::c5d_full_gjk_quad_pw(
                &B0, &B1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0, dp.A1, dp.v0,
                dp.v1, slackness, thickness, K, max_iter, num_iter, t_history);
            // std::cout << "C5D-Quad-PW toi: " << toi << std::endl;
            timer.stop();

            // Compute num_iters to reach key ratios
            for (int k = 0; k < n_key_iters; k++) {
                key_iters[k] = std::lower_bound(t_history, t_history + num_iter,
                                                key_ratios[k] * toi) -
                               t_history;
                for (int iter = 0; iter < key_iters[k]; iter++) {
                    n_not_converged_per_iter[2][k][iter]++;
                }
            }
            std::string key_iters_str = "";
            for (int k = 0; k < n_key_iters; k++) {
                key_iters_str += std::to_string(key_iters[k]);
                if (k < n_key_iters - 1) {
                    key_iters_str += ",";
                }
            }

            std::cout << "    C5D-Quad-PW num_iter: " << num_iter
                      << "; key_iters: " << key_iters_str
                      << "; wall time: " << timer.elapsed().wall / 1e6 << " ms"
                      << "; toi: " << toi << std::endl;

            c5d_quad_pw_file << i << "," << toi << "," << num_iter << ","
                             << key_iters_str << ","
                             << timer.elapsed().wall / 1e6 << "\n";
        }

        // Test ACCD
        {
            boost::timer::auto_cpu_timer timer;
            double num_iter = 0.0;
            accd_toi = c5d::exp::additiveCCD(
                dp.V0, dp.V1, dp.F0, dp.F1, dp.P0, dp.P1, dp.x0, dp.x1, dp.A0,
                dp.A1, dp.v0, dp.v1, slackness, thickness, max_iter, num_iter);
            // std::cout << "ACCD toi: " << toi << std::endl;
            timer.stop();
            std::cout << "    ACCD wall time: " << timer.elapsed().wall / 1e6
                      << " ms" << "; toi: " << accd_toi
                      << "; num_iter: " << num_iter << std::endl;

            accd_file << i << "," << accd_toi << "," << num_iter << ",,,,"
                      << timer.elapsed().wall / 1e6 << "\n";
        }
    }

    // For plotting convergence curves
    for (int key_i = 0; key_i < n_key_iters; key_i++) {
        std::string filename = output_dir + "/curve_" +
                              std::to_string(static_cast<int>(key_ratios[key_i] * 100)) + ".csv";
        std::ofstream curve_file(filename);
        curve_file << "iter,c5d_linear,c5d_quad,c5d_quad_pw\n";
        for (int iter = 0; iter < max_iter; iter++) {
            curve_file << iter << ","
                       << n_not_converged_per_iter[0][key_i][iter] << ","
                       << n_not_converged_per_iter[1][key_i][iter] << ","
                       << n_not_converged_per_iter[2][key_i][iter] << "\n";
            if (n_not_converged_per_iter[0][key_i][iter] +
                    n_not_converged_per_iter[1][key_i][iter] +
                    n_not_converged_per_iter[2][key_i][iter] ==
                0) {
                break;
            }
        }
        curve_file.close();
    }

    c5d_linear_file.close();
    c5d_quad_pw_file.close();
    accd_file.close();

    delete[] t_history;

    return 0;
}