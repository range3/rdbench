#pragma once

#include <experimental/mdspan>

namespace rdbench::v2 {

using extent_2d = std::experimental::dextents<size_t, 2>;
using mdspan_2d = std::experimental::mdspan<double, extent_2d>;

struct domain {
  size_t total_nx;  // total number of grid points in x
  size_t total_ny;  // total number of grid points in y
  size_t nx;        // number of grid points in x
  size_t ny;        // number of grid points in y
  size_t start_x;   // starting index in x
  size_t start_y;   // starting index in y
  constexpr auto nx_with_halo() const -> size_t { return nx + 2; }
  constexpr auto ny_with_halo() const -> size_t { return ny + 2; }
  constexpr auto size_with_halo() const -> size_t {
    return nx_with_halo() * ny_with_halo();
  }
  constexpr auto size() const -> size_t { return nx * ny; }
  constexpr auto total_size() const -> size_t { return total_nx * total_ny; }
};
}  // namespace rdbench::v2
