#pragma once

#include <cstddef>
#include <vector>

#include <cxxmpi/cart_comm.hpp>
#include <cxxmpi/dtype.hpp>
#include <cxxmpi/request.hpp>

#include <experimental/mdspan>

#include "rdbench/v2/domain.hpp"
#include "rdbench/v2/initializer.hpp"

namespace rdbench::v2 {

class gray_scott {
 public:
  using parameters_type = struct parameters {
    double du = 0.05;    // diffusion rate of u
    double dv = 0.1;     // diffusion rate of v
    double f = 0.04;     // feed rate
    double k = 0.06075;  // kill rate
    double dt = 0.2;     // time step
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
        u_buf_1_(domain_.size_with_halo(), 1.0),
        v_buf_1_(domain_.size_with_halo(), 0.0),
        u_buf_2_(domain_.size_with_halo()),
        v_buf_2_(domain_.size_with_halo()),
        extent_{domain_.ny_with_halo(), domain_.nx_with_halo()},
        u_{u_buf_1_.data(), extent_},
        v_{v_buf_1_.data(), extent_},
        u_next_{u_buf_2_.data(), extent_},
        v_next_{v_buf_2_.data(), extent_} {}

  auto comm() const -> const cxxmpi::cart_comm& { return comm_; }
  auto domain() const -> domain_type { return domain_; }
  auto u() const -> mdspan_2d { return u_; }
  auto v() const -> mdspan_2d { return v_; }

  auto init_field(const field_initializer& init_u,
                  const field_initializer& init_v) -> void {
    init_u.init(u_, domain_);
    init_v.init(v_, domain_);
  }

  void step() {
    exchange_halos();
    compute_next_state();
    swap_fields();
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
                            mdspan_2d view,
                            int tag_ofs) {
    const int tag_up = tag_ofs;
    const int tag_down = tag_ofs + 1;
    const int tag_left = tag_ofs + 2;
    const int tag_right = tag_ofs + 3;

    comm_.irecv(std::span{&view(0, 1), domain_.nx}, neighbors_.up, tag_up,
                requests.add());
    comm_.irecv(std::span{&view(domain_.ny_with_halo() - 1, 1), domain_.nx},
                neighbors_.down, tag_down, requests.add());
    comm_.irecv(std::span<double, 1>{&view(1, 0), 1},
                cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.left, tag_left,
                requests.add());
    comm_.irecv(std::span<double, 1>{&view(1, domain_.nx_with_halo() - 1), 1},
                cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.right, tag_right,
                requests.add());
  }

  void exchange_halos_send(mdspan_2d view, int tag_ofs) {
    const int tag_up = tag_ofs;
    const int tag_down = tag_ofs + 1;
    const int tag_left = tag_ofs + 2;
    const int tag_right = tag_ofs + 3;

    comm_.send(std::span<const double>{&view(1, 1), domain_.nx}, neighbors_.up,
               tag_down);
    comm_.send(std::span<const double>{&view(domain_.ny_with_halo() - 2, 1),
                                       domain_.nx},
               neighbors_.down, tag_up);
    comm_.send(std::span<const double, 1>{&view(1, 1), 1},
               cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.left, tag_right);
    comm_.send(
        std::span<const double, 1>{&view(1, domain_.nx_with_halo() - 2), 1},
        cxxmpi::weak_dtype{halo_type_}, 1, neighbors_.right, tag_left);
  }

  void compute_next_state() {
    const auto nx = domain_.nx;
    const auto ny = domain_.ny;
    const auto du = params_.du;
    const auto dv = params_.dv;
    const auto f = params_.f;
    const auto k = params_.k;
    const auto dt = params_.dt;

    for (size_t y = 1; y < ny + 1; ++y) {
      for (size_t x = 1; x < nx + 1; ++x) {
        const auto u = u_(y, x);
        const auto v = v_(y, x);
        const auto laplacian_u =
            u_(y, x - 1) + u_(y, x + 1) + u_(y - 1, x) + u_(y + 1, x) - 4 * u;
        const auto laplacian_v =
            v_(y, x - 1) + v_(y, x + 1) + v_(y - 1, x) + v_(y + 1, x) - 4 * v;
        const auto uv2 = u * v * v;
        const auto u_react = -uv2 + f * (1 - u);
        const auto v_react = uv2 - (f + k) * v;

        u_next_(y, x) = u + dt * (du * laplacian_u + u_react);
        v_next_(y, x) = v + dt * (dv * laplacian_v + v_react);
      }
    }
  }

  void swap_fields() {
    std::swap(u_, u_next_);
    std::swap(v_, v_next_);
  }
};
}  // namespace rdbench::v2
