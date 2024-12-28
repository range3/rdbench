#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>

#include <cxxmpi/comm.hpp>
#include <cxxmpi/universe.hpp>
#include <kamping/communicator.hpp>
#include <kamping/environment.hpp>

#include "rdbench/v2/options.hpp"

auto main(int argc, char** argv) -> int {
  try {
    // const kamping::Environment env;
    // auto const& comm = kamping::comm_world();
    const cxxmpi::universe univ(argc, argv);
    auto const& comm = cxxmpi::comm_world();
    std::cout << std::format("Hello from rank {} of {}!\n", comm.rank(),
                             comm.size());
    auto opt = rdbench::v2::options::parse(argc, argv);
    opt.validate_parameters();

    if (comm.size() == opt.nr_tiles_x * opt.nr_tiles_y) {
      throw std::runtime_error(
          std::format("np ({}) must be equal to total number of tiles ({})",
                      comm.size(), opt.nr_tiles_x * opt.nr_tiles_y));
    }

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
