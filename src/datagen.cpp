#include "Dataset.hpp"

int main(int argc, char *argv[]) {
    // spdlog::set_level(spdlog::level::debug);
    spdlog::set_level(spdlog::level::info);

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <exp_config_file> <output_dir>"
                  << std::endl;
        return 1;
    }
    spdlog::info("Reading experiment config file: {}", argv[1]);
    YAML::Node config_node = YAML::LoadFile(argv[1]);
    c5d::exp::DataConfig data_config;
    data_config.loadYaml(config_node);
    spdlog::info("Experiment config loaded:");
    std::cout << data_config << std::endl;

    std::string exp_name = argv[1];
    // extract the filename without extension
    exp_name = exp_name.substr(exp_name.find_last_of("/\\") + 1);
    exp_name = exp_name.substr(0, exp_name.find_last_of('.'));

    c5d::exp::Dataset dataset;
    dataset.generate_from_config(data_config, true,
                                 std::string(argv[2]) + "/" + exp_name);
    spdlog::info("Generated dataset with {} data points",
                 dataset.datapoints.size());

    // spdlog::info("Saving dataset to {}/{}", argv[2], exp_name);
    // dataset.save(std::string(argv[2]) + "/" + exp_name);

    std::cout << "P0: " << dataset.datapoints[0].P0 << std::endl;
    std::cout << "P1: " << dataset.datapoints[0].P1 << std::endl;
    std::cout << "x0: " << dataset.datapoints[0].x0 << std::endl;
    std::cout << "x1: " << dataset.datapoints[0].x1 << std::endl;
    std::cout << "A0: " << dataset.datapoints[0].A0 << std::endl;
    std::cout << "A1: " << dataset.datapoints[0].A1 << std::endl;
    std::cout << "v0: " << dataset.datapoints[0].v0 << std::endl;
    std::cout << "v1: " << dataset.datapoints[0].v1 << std::endl;
    std::cout << "toi: " << dataset.datapoints[0].toi << std::endl;

    return 0;
}
