#pragma once

#include <algorithm>

#include <experimental/mdspan>

#include "rdbench/v2/domain.hpp"

namespace rdbench::v2 {
class tile_filler {
 public:
  using domain_type = struct domain;

  tile_filler() = default;
  tile_filler(const tile_filler&) = default;
  tile_filler(tile_filler&&) = default;
  auto operator=(const tile_filler&) -> tile_filler& = default;
  auto operator=(tile_filler&&) -> tile_filler& = default;

  virtual ~tile_filler() = default;
  virtual void apply(mdspan_2d tile, const domain_type& domain) const = 0;
};

class center_block_filler : public tile_filler {
 public:
  explicit center_block_filler(double value,
                               size_t block_size_x,
                               size_t block_size_y)
      : value_{value}, block_nx_{block_size_x}, block_ny_{block_size_y} {}

  void apply(mdspan_2d tile, const domain_type& domain) const override {
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

class constant_filler : public tile_filler {
 public:
  explicit constant_filler(double value) : value_{value} {}

  void apply(mdspan_2d tile, const domain_type& /*domain*/) const override {
    std::fill(tile.data_handle(), tile.data_handle() + tile.size(), value_);
  }

 private:
  double value_{0.0};
};

}  // namespace rdbench::v2
