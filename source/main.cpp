#include <exception>
#include <iostream>

#include <cxxmpi/error.hpp>
#include <cxxmpi/universe.hpp>

#include "rdbench/v2/driver.hpp"
#include "rdbench/v2/options.hpp"

// inline constexpr auto IO_EVENT = prof::event_name{"IO"};

auto main(int argc, char** argv) -> int {
  try {
    // auto& recorder = prof::recorder::instance();

    // recorder.register_event(prof::singleton<prof::event<IO_EVENT>>::instance());
    const cxxmpi::universe univ(argc, argv);
    auto opt = rdbench::v2::options::parse(argc, argv);
    auto driver = rdbench::v2::gray_scott_driver(opt);
    driver.run();
  } catch (const cxxmpi::mpi_error& e) {
    std::cerr << "MPI error: " << e.what() << '\n';
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
