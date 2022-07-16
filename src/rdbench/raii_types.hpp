// SPDX-License-Identifier: Apache-2.0
// Copyright 2022 range3 ( https://github.com/range3 )
#pragma once

#include <mpi.h>

#include <cstddef>
#include <memory>

struct DatatypeWrap {
  DatatypeWrap() = default;
  DatatypeWrap(std::nullptr_t) {}
  explicit DatatypeWrap(MPI_Datatype dt) : dt(dt) {}
  explicit operator bool() const { return dt != MPI_DATATYPE_NULL; }
  friend bool operator==(DatatypeWrap l, DatatypeWrap r) { return l.dt == r.dt; }
  friend bool operator!=(DatatypeWrap l, DatatypeWrap r) { return !(l == r); }
  MPI_Datatype dt = MPI_DATATYPE_NULL;
};

class Datatype {
public:
  Datatype() = default;
  Datatype(MPI_Datatype dt) { dt_ = decltype(dt_)(DatatypeWrap(dt), Deleter()); }
  MPI_Datatype get() { return dt_.get().dt; }

private:
  struct Deleter {
    using pointer = DatatypeWrap;
    void operator()(DatatypeWrap dt) { MPI_Type_free(&dt.dt); }
  };
  std::unique_ptr<void, Deleter> dt_;
};

struct CommWrap {
  CommWrap() = default;
  CommWrap(std::nullptr_t) {}
  explicit CommWrap(MPI_Comm comm) : comm(comm) {}
  explicit operator bool() const { return comm != MPI_COMM_NULL; }
  friend bool operator==(CommWrap l, CommWrap r) { return l.comm == r.comm; }
  friend bool operator!=(CommWrap l, CommWrap r) { return !(l == r); }
  MPI_Comm comm = MPI_COMM_NULL;
};

struct Comm {
public:
  Comm() = default;
  Comm(MPI_Comm comm) { comm_ = decltype(comm_)(CommWrap(comm), Deleter()); }
  MPI_Comm get() { return comm_.get().comm; }

private:
  struct Deleter {
    using pointer = CommWrap;
    void operator()(CommWrap comm) { MPI_Comm_free(&comm.comm); }
  };
  std::unique_ptr<void, Deleter> comm_;
};
