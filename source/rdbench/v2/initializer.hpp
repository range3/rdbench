#pragma once

#include <experimental/mdspan>

#include "rdbench/v2/domain.hpp"

namespace rdbench::v2 {
class tile_initializer {
 public:
  using domain_type = struct domain;

  tile_initializer() = default;
  tile_initializer(const tile_initializer&) = default;
  tile_initializer(tile_initializer&&) = default;
  auto operator=(const tile_initializer&) -> tile_initializer& = default;
  auto operator=(tile_initializer&&) -> tile_initializer& = default;

  virtual ~tile_initializer() = default;
  virtual void init(mdspan_2d tile, const domain_type& domain) const = 0;
};

class center_block_initializer : public tile_initializer {
 public:
  explicit center_block_initializer(double value,
                                    size_t block_size_x,
                                    size_t block_size_y)
      : value_{value}, block_nx_{block_size_x}, block_ny_{block_size_y} {}

  void init(mdspan_2d tile, const domain_type& domain) const override {
    const auto block_nx = std::min(block_nx_, domain.total_nx);
    const auto block_ny = std::min(block_ny_, domain.total_ny);
    const auto start_x = domain.start_x;
    const auto start_y = domain.start_y;
    const auto block_start_x = domain.total_nx / 2 - block_nx / 2;
    const auto block_start_y = domain.total_ny / 2 - block_ny / 2;
    const auto block_end_x = block_start_x + block_ny_;
    const auto block_end_y = block_start_y + block_ny_;

    for (size_t y = 0; y < domain.ny; ++y) {
      for (size_t x = 0; x < domain.nx; ++x) {
        if (block_start_x <= start_x + x && start_x + x < block_end_x
            && block_start_y <= start_y + y && start_y + y < block_end_y) {
          tile(y + 1, x + 1) = value_;
        }
      }
    }
  }

 private:
  double value_{0.0};
  size_t block_nx_{0};
  size_t block_ny_{0};
};

}  // namespace rdbench::v2
