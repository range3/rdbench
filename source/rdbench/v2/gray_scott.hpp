#pragma once

#include <cstddef>
// #include <execution>
// #include <ranges>
#include <vector>

#include <cxxmpi/cart_comm.hpp>
#include <cxxmpi/dtype.hpp>
#include <cxxmpi/request.hpp>

#include <experimental/mdspan>

#include "rdbench/v2/domain.hpp"
#include "rdbench/v2/io.hpp"
#include "rdbench/v2/tile_filler.hpp"

namespace rdbench::v2 {

class gray_scott {
 public:
  using parameters_type = struct parameters {
    double du;  // diffusion rate of u
    double dv;  // diffusion rate of v
    double f;   // feed rate
    double k;   // kill rate
    double dt;  // time step
    // dx = dy = 1.0     // grid spacing
  };

  using domain_type = struct domain;

 private:
  cxxmpi::cart_comm comm_;
  parameters_type params_;
  domain_type domain_;
  cxxmpi::neighbors_2d neighbors_;

  cxxmpi::dtype halo_type_;

  std::vector<double> u_buf_1_;  // activator
  std::vector<double> v_buf_1_;  // inhibitor
  std::vector<double> u_buf_2_;
  std::vector<double> v_buf_2_;

  extent_2d extent_;
  mdspan_2d u_;
  mdspan_2d v_;
  mdspan_2d u_next_;
  mdspan_2d v_next_;

 public:
  gray_scott(cxxmpi::cart_comm comm,
             const parameters_type& params,
             size_t sz_tile_x,
             size_t sz_tile_y)
      : comm_{std::move(comm)},
        params_{params},
        domain_{create_domain(comm_, sz_tile_x, sz_tile_y)},
        neighbors_{comm_.neighbors_2d()},
        halo_type_{create_halo_type()},
        u_buf_1_(domain_.size_with_halo(), 0.0),
        v_buf_1_(domain_.size_with_halo(), 0.0),
        u_buf_2_(domain_.size_with_halo()),
        v_buf_2_(domain_.size_with_halo()),
        extent_{domain_.ny_with_halo(), domain_.nx_with_halo()},
        u_{u_buf_1_.data(), extent_},
        v_{v_buf_1_.data(), extent_},
        u_next_{u_buf_2_.data(), extent_},
        v_next_{v_buf_2_.data(), extent_} {}

  // void print(std::ostream& os) const {
  //   os << std::format("domain: {}x{}+{}+{}\n", domain_.total_nx,
  //                     domain_.total_ny, domain_.start_x, domain_.start_y);
  //   os << std::format("u: {}x{}\n", u_.extent(1), u_.extent(0));
  //   os << std::format("v: {}x{}\n", v_.extent(1), v_.extent(0));
  //   os << u_buf_1_.size() << ' ' << v_buf_1_.size() << ' ' << u_buf_2_.size()
  //      << ' ' << v_buf_2_.size() << '\n';
  //   os << u_.data_handle() << ' ' << v_.data_handle() << ' '
  //      << u_next_.data_handle() << ' ' << v_next_.data_handle() << '\n';
  //   os << u_buf_1_.data() << ' ' << v_buf_1_.data() << ' ' << u_buf_2_.data()
  //      << ' ' << v_buf_2_.data() << '\n';
  //   os << std::fixed << std::setprecision(2);
  //   for (size_t y = 0; y < u_.extent(0); ++y) {
  //     for (size_t x = 0; x < u_.extent(1); ++x) {
  //       os << u_(y, x) << ' ';
  //     }
  //     os << '\n';
  //   }
  // }

  auto comm() const -> const cxxmpi::cart_comm& { return comm_; }
  auto domain() const -> domain_type { return domain_; }
  auto u() const -> mdspan_2d { return u_; }
  auto v() const -> mdspan_2d { return v_; }

  auto apply_tile_filler(const tile_filler& filler, data_type type) -> void {
    if (type == data_type::u) {
      filler.apply(u_, domain_);
    } else {
      filler.apply(v_, domain_);
    }
  }

  void step() {
    exchange_halos();
    // compute_next_state_stdpar();
    compute_next_state();
    swap_tiles();
  }

  void ckpt(const io_strategy& io, size_t idx, data_type type) const {
    if (type == data_type::u) {
      io.write(cxxmpi::weak_comm{comm_}, u_, domain_, type, idx);
    } else {
      io.write(cxxmpi::weak_comm{comm_}, v_, domain_, type, idx);
    }
  }

 private:
  static auto create_domain(const cxxmpi::cart_comm& comm,
                            size_t sz_tile_x,
                            size_t sz_tile_y) -> domain_type {
    const auto& dims = comm.dims();
    const auto& coords = comm.coords();
    const auto dims2 = std::array<size_t, 2>{static_cast<size_t>(dims[0]),
                                             static_cast<size_t>(dims[1])};
    const auto coords2 = std::array<size_t, 2>{static_cast<size_t>(coords[0]),
                                               static_cast<size_t>(coords[1])};
    return domain_type{
        .total_nx = dims2[1] * sz_tile_x,
        .total_ny = dims2[0] * sz_tile_y,
        .nx = sz_tile_x,
        .ny = sz_tile_y,
        .start_x = coords2[1] * sz_tile_x,
        .start_y = coords2[0] * sz_tile_y,
    };
  }

  auto create_halo_type() const -> cxxmpi::dtype {
    auto halo_type = cxxmpi::dtype{
        cxxmpi::as_weak_dtype<double>(),
        static_cast<int>(domain_.ny),
        1,
        static_cast<int>(domain_.nx + 2),
    };
    halo_type.commit();
    return halo_type;
  }

  void exchange_halos() {
    cxxmpi::request_group requests;
    exchange_halos_irecv(requests, u_, 100);
    exchange_halos_irecv(requests, v_, 200);
    exchange_halos_send(u_, 100);
    exchange_halos_send(v_, 200);
    requests.wait_all_without_status();
  }

  void exchange_halos_irecv(cxxmpi::request_group& requests,
                            mdspan_2d tile,
                            int tag_ofs) {
    const int tag_up = tag_ofs;
    const int tag_down = tag_ofs + 1;
    const int tag_left = tag_ofs + 2;
    const int tag_right = tag_ofs + 3;

    comm_.irecv(std::span{&tile(0, 1), domain_.nx}, neighbors_.up, tag_up,
                requests.add());
    comm_.irecv(std::span{&tile(domain_.ny_with_halo() - 1, 1), domain_.nx},
                neighbors_.down, tag_down, requests.add());
    comm_.irecv(std::span<double, 1>{&tile(1, 0), 1},
                cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.left, tag_left,
                requests.add());
    comm_.irecv(std::span<double, 1>{&tile(1, domain_.nx_with_halo() - 1), 1},
                cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.right, tag_right,
                requests.add());
  }

  void exchange_halos_send(mdspan_2d tile, int tag_ofs) {
    const int tag_up = tag_ofs;
    const int tag_down = tag_ofs + 1;
    const int tag_left = tag_ofs + 2;
    const int tag_right = tag_ofs + 3;

    comm_.send(std::span<const double>{&tile(1, 1), domain_.nx}, neighbors_.up,
               tag_down);
    comm_.send(std::span<const double>{&tile(domain_.ny_with_halo() - 2, 1),
                                       domain_.nx},
               neighbors_.down, tag_up);
    comm_.send(std::span<const double, 1>{&tile(1, 1), 1},
               cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.left, tag_right);
    comm_.send(
        std::span<const double, 1>{&tile(1, domain_.nx_with_halo() - 2), 1},
        cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.right, tag_left);
  }

  struct compute_kernel {
    double du, dv, f, k, dt;
    mdspan_2d u, v, u_next, v_next;

    void operator()(int y, int x) const {
      const auto current_u = u(y, x);
      const auto current_v = v(y, x);

      // Laplacian
      const auto laplacian_u =
          u(y, x - 1) + u(y, x + 1) + u(y - 1, x) + u(y + 1, x) - 4 * current_u;
      const auto laplacian_v =
          v(y, x - 1) + v(y, x + 1) + v(y - 1, x) + v(y + 1, x) - 4 * current_v;

      // Reaction
      const auto uv2 = current_u * current_v * current_v;
      const auto u_react = -uv2 + f * (1 - current_u);
      const auto v_react = uv2 - (f + k) * current_v;

      // Update
      u_next(y, x) = current_u + dt * (du * laplacian_u + u_react);
      v_next(y, x) = current_v + dt * (dv * laplacian_v + v_react);
    }
  };

  void compute_next_state() {
    const auto nx = static_cast<int>(domain_.nx);
    const auto ny = static_cast<int>(domain_.ny);
    const auto kernel = compute_kernel{
        .du = params_.du,
        .dv = params_.dv,
        .f = params_.f,
        .k = params_.k,
        .dt = params_.dt,
        .u = u_,
        .v = v_,
        .u_next = u_next_,
        .v_next = v_next_,
    };

    for (int y = 1; y < ny + 1; ++y) {
      for (int x = 1; x < nx + 1; ++x) {
        kernel(y, x);
      }
    }
  }

#ifdef __cpp_lib_ranges_cartesian_product
  void compute_next_state_stdpar() {
    const auto nx = static_cast<int>(domain_.nx);
    const auto ny = static_cast<int>(domain_.ny);
    const auto kernel = compute_kernel{
        .du = params_.du,
        .dv = params_.dv,
        .f = params_.f,
        .k = params_.k,
        .dt = params_.dt,
        .u = u_,
        .v = v_,
        .u_next = u_next_,
        .v_next = v_next_,
    };

    // work around for std::views::cartesian_product
    // std::vector<std::pair<int, int>> indices;
    // indices.reserve(static_cast<size_t>(nx * ny));
    // for (auto y : std::views::iota(1, static_cast<int>(ny + 1))) {
    //   for (auto x : std::views::iota(1, static_cast<int>(nx + 1))) {
    //     indices.emplace_back(y, x);
    //   }
    // }

    auto indices = std::views::cartesian_product(
        std::views::iota(1, static_cast<int>(ny + 1)),
        std::views::iota(1, static_cast<int>(nx + 1)));

    std::for_each(std::execution::par, indices.begin(), indices.end(),
                  [kernel](auto idx) {
                    auto [y, x] = idx;
                    kernel(y, x);
                  });
  }
#endif

  void swap_tiles() {
    std::swap(u_, u_next_);
    std::swap(v_, v_next_);
  }
};
}  // namespace rdbench::v2
