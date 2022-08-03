#pragma once
#include <fmt/core.h>
#include <mpi.h>

#include <stdexcept>

class mpi_error : public std::runtime_error {
public:
  mpi_error() : runtime_error(mpi_error::msg(MPI_ERR_UNKNOWN)) {}
  mpi_error(int ecode) : runtime_error(mpi_error::msg(ecode)), ecode_(ecode) {}

  static std::string msg(int ecode) {
    int eclass, len;
    char estring[MPI_MAX_ERROR_STRING];
    MPI_Error_class(ecode, &eclass);
    MPI_Error_string(ecode, estring, &len);
    return fmt::format("MPI Error: {}: {}", eclass, estring);
  }

private:
  int ecode_ = MPI_ERR_UNKNOWN;
};
