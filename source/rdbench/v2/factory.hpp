#pragma once

#include <memory>

#include <cxxmpi/cart_comm.hpp>
#include <cxxmpi/comm.hpp>
#include <cxxmpi/dims.hpp>
#include <fmt/format.h>

#include "rdbench/v2/domain.hpp"
#include "rdbench/v2/gray_scott.hpp"
#include "rdbench/v2/options.hpp"
#include "rdbench/v2/tile_filler.hpp"

namespace rdbench::v2 {

class gray_scott_factory {
 public:
  static auto create(const options& opts) -> std::unique_ptr<gray_scott> {
    validate(opts);
    auto model = std::make_unique<gray_scott>(create_cart_comm(opts),
                                              gray_scott::parameters{
                                                  opts.param_du,
                                                  opts.param_dv,
                                                  opts.param_f,
                                                  opts.param_k,
                                                  opts.param_dt,
                                              },
                                              opts.sz_tile_x, opts.sz_tile_y);
    model->apply_tile_filler(constant_filler{1.0}, data_type::u);
    // model->apply_tile_filler(center_block_filler{0.9, 12, 12}, data_type::u);
    model->apply_tile_filler(center_block_filler{0.7, 6, 6}, data_type::v);
    return model;
  }

 private:
  static void validate(const options& opts) {
    opts.validate_parameters();
    validate_model_parameters(opts);
  }

  static auto create_cart_comm(const options& opts) -> cxxmpi::cart_comm {
    const auto& world = cxxmpi::comm_world();

    auto dims = cxxmpi::create_dims(
        static_cast<int>(world.size()),
        {static_cast<int>(opts.nr_tiles_y), static_cast<int>(opts.nr_tiles_x)});

    if (static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1])
        != world.size()) {
      throw std::runtime_error(
          fmt::format("Invalid number of processes for {}x{} grid: {}",
                      opts.nr_tiles_x, opts.nr_tiles_y, world.size()));
    }

    return cxxmpi::cart_comm{
        world,
        {static_cast<size_t>(dims[0]), static_cast<size_t>(dims[1])},
        {true, true},
        true};
  }

  static void validate_model_parameters(const rdbench::v2::options& opts) {
    if (opts.param_f < 0.0) {
      throw std::runtime_error("Feed rate F must be positive");
    }
    if (opts.param_k < 0.0) {
      throw std::runtime_error("Kill rate k must be positive");
    }
    if (opts.param_du <= 0.0 || opts.param_dv <= 0.0) {
      throw std::runtime_error("Diffusion coefficients must be positive");
    }

    if (opts.param_f > 0.1 || opts.param_k > 0.1) {
      std::cerr << "Warning: F and k values above 0.1 might not produce "
                   "typical Turing patterns\n";
    }

    if (opts.param_du <= opts.param_dv) {
      std::cerr << "Warning: Du should typically be greater than Dv for pattern "
                   "formation\n";
    }

    // Check numerical stability condition (CFL condition)
    constexpr double dx = 1.0;  // grid spacing
    const double max_diffusion = std::max(opts.param_du, opts.param_dv);
    const double max_stable_dt = 0.5 * dx * dx / max_diffusion;

    if (opts.param_dt <= 0.0) {
      throw std::runtime_error("Time step dt must be positive");
    }
    if (opts.param_dt > max_stable_dt) {
      throw std::runtime_error(fmt::format(
          "Time step {} exceeds maximum stable value {} (CFL condition)",
          opts.param_dt, max_stable_dt));
    }
  }
};

}  // namespace rdbench::v2
