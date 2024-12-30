#pragma once

#include <array>
#include <format>
#include <functional>
#include <string>

#include <cxxmpi/dtype.hpp>
#include <cxxmpi/file.hpp>

#include "rdbench/v2/gray_scott.hpp"
#include "rdbench/v2/options.hpp"

namespace rdbench::v2 {

class io_strategy {
 public:
  explicit io_strategy(const options& opts,
                       const gray_scott::domain_type& domain)
      : opts_{opts}, domain_type_{create_domain_dtype(domain)} {}
  constexpr io_strategy(const io_strategy&) = delete;
  auto operator=(const io_strategy&) -> io_strategy& = delete;
  constexpr io_strategy(io_strategy&&) noexcept = default;
  auto operator=(io_strategy&&) noexcept -> io_strategy& = default;
  virtual ~io_strategy() = default;

  void write_model_checkpoint(const gray_scott& model, size_t idx) {
    auto file =
        cxxmpi::open(file_path(model, idx), cxxmpi::weak_comm{model.comm()},
                     MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN);
    file.set_atomicity(false);
    write(cxxmpi::weak_file{file}, model);

    if (!opts().nosync) {
      file.sync();
    }
  }

 protected:
  virtual void write(cxxmpi::weak_file, const gray_scott&) const = 0;

  auto opts() const -> const options& { return opts_.get(); }
  auto weak_domain_type() const -> cxxmpi::weak_dtype {
    return cxxmpi::weak_dtype{domain_type_};
  }

 private:
  std::reference_wrapper<const options> opts_;
  cxxmpi::dtype domain_type_;

  auto file_path(const gray_scott& model, size_t idx) const -> std::string {
    return std::format("{}{}", opts().output, filename(model.domain(), idx));
  }

  static constexpr auto filename(const gray_scott::domain_type& domain,
                                 size_t idx) -> std::string {
    return std::format("{}x{}-{}x{}-{:06d}.bin", domain.total_nx,
                       domain.total_ny, domain.nx, domain.ny, idx);
  }

  static auto create_domain_dtype(const gray_scott::domain_type& domain)
      -> cxxmpi::dtype {
    auto data_type = cxxmpi::dtype{
        cxxmpi::as_weak_dtype<double>(),
        std::array<int, 2>{static_cast<int>(domain.ny_with_halo()),
                           static_cast<int>(domain.nx_with_halo())},
        std::array<int, 2>{static_cast<int>(domain.ny),
                           static_cast<int>(domain.nx)},
        std::array<int, 2>{static_cast<int>(1), static_cast<int>(1)}};
    data_type.commit();
    return data_type;
  }
};

class canonical_io final : public io_strategy {
 public:
  explicit canonical_io(const options& opts,
                        const gray_scott::domain_type& domain)
      : io_strategy(opts, domain),
        file_view_type_{create_file_view_type(domain)} {}
  constexpr canonical_io(const canonical_io&) = delete;
  auto operator=(const canonical_io&) -> canonical_io& = delete;
  constexpr canonical_io(canonical_io&&) noexcept = default;
  auto operator=(canonical_io&&) noexcept -> canonical_io& = default;
  ~canonical_io() override = default;

 protected:
  void write(cxxmpi::weak_file file, const gray_scott& model) const override {
    file.set_view(0, cxxmpi::as_weak_dtype<double>(),
                  cxxmpi::weak_dtype{file_view_type_});
    if (opts().collective) {
      file.write_at_all(
          0, std::span<const double>{model.u().data_handle(), model.u().size()},
          weak_domain_type());
    } else {
      file.write_at(
          0, std::span<const double>{model.u().data_handle(), model.u().size()},
          weak_domain_type());
    }
  }

 private:
  cxxmpi::dtype file_view_type_;

  static auto create_file_view_type(const gray_scott::domain_type& domain)
      -> cxxmpi::dtype {
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

class log_io final : public io_strategy {
 public:
  explicit log_io(const options& opts, const gray_scott::domain_type& domain)
      : io_strategy(opts, domain) {}
  constexpr log_io(const log_io&) = delete;
  auto operator=(const log_io&) -> log_io& = delete;
  constexpr log_io(log_io&&) noexcept = default;
  auto operator=(log_io&&) noexcept -> log_io& = default;
  ~log_io() override = default;

 protected:
  void write(cxxmpi::weak_file file, const gray_scott& model) const override {
    MPI_Offset const ofs = static_cast<int64_t>(sizeof(double))
                         * static_cast<int64_t>(model.domain().size())
                         * model.comm().rank();
    if (opts().collective) {
      file.write_at_all(
          ofs,
          std::span<const double>{model.u().data_handle(), model.u().size()},
          weak_domain_type());
    } else {
      file.write_at(
          ofs,
          std::span<const double>{model.u().data_handle(), model.u().size()},
          weak_domain_type());
    }
  }
};

}  // namespace rdbench::v2
