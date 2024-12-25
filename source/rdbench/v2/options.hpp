#pragma once

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <cxxopts.hpp>

namespace rdbench::v2 {

enum class file_layout {
  canonical,  // Standard data layout for visualization
  log,        // Sequential log format for analysis
};

struct options {
  // Tile partitioning settings
  size_t nr_tiles_x = 1;
  size_t nr_tiles_y = 1;
  size_t sz_tile_x = 1024;
  size_t sz_tile_y = 1024;

  // Simulation settings
  size_t steps = 20000;
  size_t interval = 200;
  bool init_output = false;
  bool verbose = false;

  // I/O settings
  std::string output = "out/o_";
  file_layout layout = file_layout::canonical;
  bool collective = false;
  bool nosync = false;
  bool validate = false;
  bool prettify = false;

  // Gray-Scott model parameters
  double param_f = 0.04;
  double param_k = 0.06075;
  double param_dt = 0.2;
  double param_du = 0.05;
  double param_dv = 0.1;

  // NOLINTNEXTLINE
  static auto parse(int argc, const char* argv[]) -> options {
    options opts;

    auto options = cxxopts::Options(
        "rdbench",
        "Benchmark for MPI-IO using Gray-Scott reaction-diffusion system");

    // clang-format off
    options.add_options()
      // Basic options
      ("h,help", "Print usage")
      ("V,version", "Print version")
      ("v,verbose", "Enable verbose output")
      ("prettify", "Enable JSON output formatting")

      // Tile partitioning settings
      ("nr_tiles_x", "Number of tiles in x direction",
       cxxopts::value<size_t>()->default_value("1"))
      ("nr_tiles_y", "Number of tiles in y direction",
       cxxopts::value<size_t>()->default_value("1"))
      ("sz_tile_x", "Size of tile in x direction",
       cxxopts::value<size_t>()->default_value("128"))
      ("sz_tile_y", "Size of tile in y direction",
       cxxopts::value<size_t>()->default_value("128"))

      // Simulation control
      ("s,steps", "Total number of simulation steps",
       cxxopts::value<size_t>()->default_value("20000"))
      ("i,interval", "Output interval steps",
       cxxopts::value<size_t>()->default_value("200"))
      ("init_output", "Enable initial state output")

      // I/O settings
      ("o,output", "Output file prefix",
       cxxopts::value<std::string>()->default_value("out/o_"))
      ("file_layout", "File layout type (canonical/log)",
       cxxopts::value<std::string>()->default_value("canonical"))
      ("collective", "Enable collective I/O")
      ("nosync", "Disable sync after write")
      ("validate", "Enable write validation")

      // Gray-Scott model parameters
      ("param_F", "Feed rate parameter F",
       cxxopts::value<double>()->default_value("0.04"))
      ("param_k", "Kill rate parameter k",
       cxxopts::value<double>()->default_value("0.06075"))
      ("param_dt", "Time step size",
       cxxopts::value<double>()->default_value("0.2"))
      ("param_Du", "Diffusion coefficient of U",
       cxxopts::value<double>()->default_value("0.05"))
      ("param_Dv", "Diffusion coefficient of V",
       cxxopts::value<double>()->default_value("0.1"));
    // clang-format on

    try {
      auto result = options.parse(argc, argv);

      if (result.count("help") != 0U) {
        std::cout << options.help() << std::endl;
        std::exit(0);  // NOLINT
      }

      if (result.count("version") != 0U) {
        std::cout << "rdbench v2.0.0" << std::endl;
        std::exit(0);  // NOLINT
      }

      // Parse basic settings
      opts.verbose = result.count("verbose") != 0U;
      opts.prettify = result.count("prettify") != 0U;

      // Parse tile settings
      opts.nr_tiles_x = result["nr_tiles_x"].as<size_t>();
      opts.nr_tiles_y = result["nr_tiles_y"].as<size_t>();
      opts.sz_tile_x = result["sz_tile_x"].as<size_t>();
      opts.sz_tile_y = result["sz_tile_y"].as<size_t>();

      // Parse simulation settings
      opts.steps = result["steps"].as<size_t>();
      opts.interval = result["interval"].as<size_t>();
      opts.init_output = result.count("init_output") != 0U;

      // Parse I/O settings
      opts.output = result["output"].as<std::string>();
      auto layout = result["file_layout"].as<std::string>();
      if (layout == "canonical") {
        opts.layout = file_layout::canonical;
      } else if (layout == "log") {
        opts.layout = file_layout::log;
      } else {
        throw std::runtime_error("Invalid file layout: " + layout);
      }

      opts.collective = result.count("collective") != 0U;
      opts.nosync = result.count("nosync") != 0U;
      opts.validate = result.count("validate") != 0U;

      // Parse Gray-Scott model parameters
      opts.param_f = result["param_F"].as<double>();
      opts.param_k = result["param_k"].as<double>();
      opts.param_dt = result["param_dt"].as<double>();
      opts.param_du = result["param_Du"].as<double>();
      opts.param_dv = result["param_Dv"].as<double>();

      return opts;

    } catch (const cxxopts::exceptions::exception& e) {
      throw std::runtime_error(std::string("Error parsing options: ")
                               + e.what());
    }
  }

  void validate_parameters() const {
    if (param_f < 0.0 || param_k < 0.0 || param_dt <= 0.0 || param_du <= 0.0
        || param_dv <= 0.0) {
      throw std::runtime_error("Model parameters must be positive");
    }
  }
};

}  // namespace rdbench::v2
