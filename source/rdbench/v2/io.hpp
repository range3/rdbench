#pragma once
#include <array>
#include <format>
#include <functional>
#include <string>

#include <cxxmpi/dtype.hpp>
#include <cxxmpi/file.hpp>

#include "rdbench/v2/domain.hpp"
#include "rdbench/v2/options.hpp"

namespace rdbench::v2 {

class io_strategy {
 public:
  explicit io_strategy(const options& opts, const domain& domain)
      : opts_{opts}, tile_dtype_{create_tile_dtype(domain)} {}

  constexpr io_strategy(const io_strategy&) = delete;
  auto operator=(const io_strategy&) -> io_strategy& = delete;
  constexpr io_strategy(io_strategy&&) noexcept = default;
  auto operator=(io_strategy&&) noexcept -> io_strategy& = default;
  virtual ~io_strategy() = default;

  virtual void write(cxxmpi::weak_comm comm,
                     cmdspan_2d data,
                     const domain& domain,
                     data_type type,
                     size_t idx) const = 0;

 protected:
  auto opts() const -> const options& { return opts_.get(); }

  auto weak_tile_type() const -> cxxmpi::weak_dtype {
    return cxxmpi::weak_dtype{tile_dtype_};
  }

 private:
  std::reference_wrapper<const options> opts_;
  cxxmpi::dtype tile_dtype_;

  static auto create_tile_dtype(const domain& domain) -> cxxmpi::dtype {
    auto tile_dtype = cxxmpi::dtype{
        cxxmpi::as_weak_dtype<double>(),
        std::array<int, 2>{static_cast<int>(domain.ny_with_halo()),
                           static_cast<int>(domain.nx_with_halo())},
        std::array<int, 2>{static_cast<int>(domain.ny),
                           static_cast<int>(domain.nx)},
        std::array<int, 2>{static_cast<int>(1), static_cast<int>(1)}};
    tile_dtype.commit();
    return tile_dtype;
  }
};

class file_io_strategy : public io_strategy {
 public:
  explicit file_io_strategy(const options& opts, const domain& domain)
      : io_strategy(opts, domain) {}
  constexpr file_io_strategy(const file_io_strategy&) = delete;
  auto operator=(const file_io_strategy&) -> file_io_strategy& = delete;
  constexpr file_io_strategy(file_io_strategy&&) noexcept = default;
  auto operator=(file_io_strategy&&) noexcept -> file_io_strategy& = default;
  ~file_io_strategy() override = default;

  void write(cxxmpi::weak_comm comm,
             cmdspan_2d data,
             const domain& domain,
             data_type type,
             size_t idx) const override {
    auto file =
        cxxmpi::open(file_path(domain, type, idx), comm,
                     MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN);
    file.set_atomicity(false);

    write_file(comm, cxxmpi::weak_file{file}, data, domain);

    if (!opts().nosync) {
      file.sync();
    }
  }

 protected:
  virtual void write_file(cxxmpi::weak_comm comm,
                          cxxmpi::weak_file file,
                          cmdspan_2d data,
                          const domain& domain) const = 0;

 private:
  auto file_path(const domain& domain,
                 data_type type,
                 size_t idx) const -> std::string {
    return std::format("{}{}", opts().output, filename(domain, type, idx));
  }

  static constexpr auto filename(const domain& domain,
                                 data_type type,
                                 size_t idx) -> std::string {
    const auto* const type_str = type == data_type::u ? "u" : "v";
    return std::format("{}-{}x{}-{}x{}-{:06d}.bin", type_str, domain.total_nx,
                       domain.total_ny, domain.nx, domain.ny, idx);
  }
};

class canonical_io final : public file_io_strategy {
 public:
  explicit canonical_io(const options& opts, const domain& domain)
      : file_io_strategy(opts, domain),
        file_view_type_{create_file_view_type(domain)} {}
  constexpr canonical_io(const canonical_io&) = delete;
  auto operator=(const canonical_io&) -> canonical_io& = delete;
  constexpr canonical_io(canonical_io&&) noexcept = default;
  auto operator=(canonical_io&&) noexcept -> canonical_io& = default;
  ~canonical_io() override = default;

  void write_file(cxxmpi::weak_comm /*comm*/,
                  cxxmpi::weak_file file,
                  cmdspan_2d data,
                  const domain& /*domain*/) const override {
    file.set_view(0, cxxmpi::as_weak_dtype<double>(),
                  cxxmpi::weak_dtype{file_view_type_});

    if (opts().collective) {
      file.write_at_all(0, data.data_handle(), 1, weak_tile_type());
    } else {
      file.write_at(0, data.data_handle(), 1, weak_tile_type());
    }
  }

 private:
  cxxmpi::dtype file_view_type_;

  static auto create_file_view_type(const domain& domain) -> cxxmpi::dtype {
    auto file_view_type =
        cxxmpi::dtype{cxxmpi::as_weak_dtype<double>(),
                      std::array<int, 2>{static_cast<int>(domain.total_ny),
                                         static_cast<int>(domain.total_nx)},
                      std::array<int, 2>{static_cast<int>(domain.ny),
                                         static_cast<int>(domain.nx)},
                      std::array<int, 2>{static_cast<int>(domain.start_y),
                                         static_cast<int>(domain.start_x)}};
    file_view_type.commit();
    return file_view_type;
  }
};

class log_io final : public file_io_strategy {
 public:
  explicit log_io(const options& opts, const domain& domain)
      : file_io_strategy(opts, domain) {}
  constexpr log_io(const log_io&) = delete;
  auto operator=(const log_io&) -> log_io& = delete;
  constexpr log_io(log_io&&) noexcept = default;
  auto operator=(log_io&&) noexcept -> log_io& = default;
  ~log_io() override = default;

  void write_file(cxxmpi::weak_comm comm,
                  cxxmpi::weak_file file,
                  cmdspan_2d data,
                  const domain& domain) const override {
    MPI_Offset const ofs = static_cast<MPI_Offset>(sizeof(double))
                         * static_cast<int64_t>(domain.size()) * comm.rank();

    if (opts().collective) {
      file.write_at_all(ofs, data.data_handle(), 1, weak_tile_type());
    } else {
      file.write_at(ofs, data.data_handle(), 1, weak_tile_type());
    }
  }
};

}  // namespace rdbench::v2
