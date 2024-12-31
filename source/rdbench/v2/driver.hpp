#pragma once

#include <iostream>
#include <memory>

#include "rdbench/v2/factory.hpp"
#include "rdbench/v2/gray_scott.hpp"
#include "rdbench/v2/io.hpp"
#include "rdbench/v2/options.hpp"

namespace rdbench::v2 {

class gray_scott_driver {
 public:
  explicit gray_scott_driver(const options& opts)
      : opts_{opts},
        model_{gray_scott_factory::create(opts)},
        io_{create_io_strategy(opts, model_->domain())} {}

  void run() {
    auto ckpt_idx = size_t{0};
    auto const total_steps = opt().steps;
    auto const interval = opt().interval;

    if (opt().init_output) {
      io_->write_model_checkpoint(*model_, ckpt_idx++);
    }

    for (size_t step = 1; step <= total_steps; ++step) {
      model_->step();

      if (interval != 0 && step % interval == 0) {
        io_->write_model_checkpoint(*model_, ckpt_idx++);
        if (opt().verbose && model_->comm().rank() == 0) {
          std::cout << "Checkpoint " << ckpt_idx - 1 << " at step " << step
                    << '\n';
        }
      }
    }
  }

 private:
  std::reference_wrapper<const options> opts_;
  std::unique_ptr<gray_scott> model_;
  std::unique_ptr<io_strategy> io_;

  auto opt() const -> const options& { return opts_.get(); }

  static auto create_io_strategy(const options& opts,
                                 const gray_scott::domain_type& domain)
      -> std::unique_ptr<io_strategy> {
    switch (opts.layout) {
      case file_layout::canonical:
        return std::make_unique<canonical_io>(opts, domain);
      case file_layout::log:
        return std::make_unique<log_io>(opts, domain);
      default:
        throw std::runtime_error("Invalid file layout");
    }
  }
};

}  // namespace rdbench::v2
