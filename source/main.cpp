#include <exception>
#include <iomanip>
#include <iostream>

#include <cxxmpi/comm.hpp>
#include <cxxmpi/error.hpp>
#include <cxxmpi/universe.hpp>

// #include "prof/formatter.hpp"
#include "prof/json_formatter.hpp"
#include "prof/recorder.hpp"
#include "rdbench/v2/driver.hpp"
#include "rdbench/v2/metrics.hpp"
#include "rdbench/v2/options.hpp"
#include "rdbench/v2/profiler.hpp"

auto main(int argc, char** argv) -> int {
  try {
    rdbench::v2::profiler::init();
    const cxxmpi::universe univ(argc, argv);
    auto opt = rdbench::v2::options::parse(argc, argv);
    auto driver = rdbench::v2::gray_scott_driver(opt);
    auto metrics = rdbench::v2::metrics_handler(
        opt, cxxmpi::weak_cart_comm{driver.model().comm()});

    driver.run();

    if (cxxmpi::comm_world().rank() == 0) {
      auto json = metrics.create_output_json(
          prof::recorder::instance().format<prof::json_formatter>());
      if (opt.prettify) {
        std::cout << std::setw(4);
      }
      std::cout << json << '\n';
    }
  } catch (const cxxmpi::mpi_error& e) {
    std::cerr << "MPI error: " << e.what() << '\n';
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
