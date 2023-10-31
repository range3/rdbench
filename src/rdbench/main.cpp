// This file is modified from ( https://github.com/kaityo256/sevendayshpc )
// under license CC-BY-4.0.
// https://github.com/kaityo256/sevendayshpc/blob/main/LICENSE

#include <fmt/core.h>
#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <nlohmann/json.hpp>
#include <numeric>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "BitmapPlusPlus.hpp"
#include "error.hpp"
#include "mpix/MPIX_interface_proposal.h"
#include "raii_types.hpp"
#include "rdbench/version.h"
#include "stopwatch.hpp"

using orderd_json = nlohmann::ordered_json;

double F = 0.04;
double k = 0.06075;
double dt = 0.2;
double Du = 0.05;
double Dv = 0.1;

typedef std::vector<double> vd;

enum class IOType { Manual, View, Chunk };
static std::unordered_map<std::string, IOType> iotype_from_string{{"manual", IOType::Manual},
                                                                  {"view", IOType::View},
                                                                  {"chunk", IOType::Chunk}};

enum class RdbenchPhase { Comm, Calc, Write, Init, Finalize };
static std::unordered_map<RdbenchPhase, std::string> phase_to_string{
    {RdbenchPhase::Comm, "comm"},         {RdbenchPhase::Calc, "calc"},
    {RdbenchPhase::Write, "write"},       {RdbenchPhase::Init, "init"},
    {RdbenchPhase::Finalize, "finalize"},
};

struct RdbenchInfo {
  int rank;
  int rank2d;
  int nprocs;
  int xnp;
  int ynp;
  std::vector<int> topo;
  int L;
  int my_grid_x;
  int my_grid_y;
  int chunk_size_x;
  int chunk_size_y;
  int my_begin_x;
  int my_begin_y;
  int my_end_x;
  int my_end_y;
  std::optional<std::string> v_bitmap_path;
  std::optional<std::string> u_bitmap_path;
  std::string output_prefix;
  std::string iot = "view";
  IOType iotype = IOType::View;
  bool collective = false;
  bool sync = true;
  bool validate = true;
  bool initial_output = true;
  bool create_output_dir = true;
  Datatype filetype;
  Datatype memtype;
  Datatype vertical_halo_type;
  Comm comm_2d;
  size_t total_steps;
  size_t interval;
  size_t write_phase_count = 0;
  size_t calc_phase_count = 1;
  bool fixed_x = false;
  bool fixed_y = false;
  bool verbose = false;
  int rank_down;
  int rank_up;
  int rank_left;
  int rank_right;

  static RdbenchInfo create(cxxopts::ParseResult &parsed) {
    RdbenchInfo info;
    MPI_Comm_size(MPI_COMM_WORLD, &info.nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &info.rank);

    {
      info.iot = parsed["iotype"].as<std::string>();
      auto i = iotype_from_string.find(info.iot);
      if (i == iotype_from_string.end()) {
        throw std::invalid_argument(fmt::format("invalid iotype {}", info.iot));
      }
      info.iotype = i->second;
    }

    info.output_prefix = parsed["output"].as<std::string>();
    info.collective = parsed.count("collective") != 0U;
    info.sync = parsed.count("nosync") == 0U;
    info.validate = parsed.count("novalidate") == 0U;
    info.initial_output = parsed.count("disable-initial-output") == 0U;
    info.create_output_dir = parsed.count("nomkdir") == 0U;
    info.fixed_x = parsed["fixed-x"].count() != 0U;
    info.fixed_y = parsed["fixed-y"].count() != 0U;
    info.total_steps = parsed["steps"].as<size_t>();
    info.interval = parsed["interval"].as<size_t>();
    info.L = parsed["L"].as<int>();
    info.verbose = parsed["verbose"].count() != 0U;
    if (parsed.count("init_v_from_bitmap") != 0U) {
      info.v_bitmap_path = parsed["init_v_from_bitmap"].as<std::string>();
    }
    if (parsed.count("init_u_from_bitmap") != 0U) {
      info.u_bitmap_path = parsed["init_u_from_bitmap"].as<std::string>();
    }

    // set Gray-Scott model parameters
    F = parsed["param_F"].as<double>();
    k = parsed["param_k"].as<double>();
    dt = parsed["param_dt"].as<double>();
    Du = parsed["param_Du"].as<double>();
    Dv = parsed["param_Dv"].as<double>();

    if (info.interval != 0) {
      info.write_phase_count = info.total_steps / info.interval + (info.initial_output ? 1 : 0);
      info.calc_phase_count = (info.total_steps + info.interval - 1) / info.interval;
    }

    int periods[] = {info.fixed_y ? 0 : 1, info.fixed_x ? 0 : 1};
    int dims[2] = {parsed["ynp"].as<int>(), parsed["xnp"].as<int>()};
    if (parsed.count("topology") == 0U) {
      MPI_Dims_create(info.nprocs, 2, dims);
      info.xnp = dims[1];
      info.ynp = dims[0];
      info.topo = {info.nprocs};

      MPI_Comm comm_2d;
      MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);
      info.comm_2d = Comm{comm_2d};
    } else {
      int ecode;
      info.topo = parsed["topology"].as<std::vector<int>>();
      std::vector<int> dims_ml(2 * info.topo.size());
      auto nprocs_in_topo = std::accumulate(info.topo.begin(), info.topo.end(), 1,
                                            [](int acc, int cur) { return acc * cur; });
      if (nprocs_in_topo != info.nprocs) {
        throw std::invalid_argument(
            "nprocs specified in --topology does not match nprocs in mpirun");
      }

      if ((ecode = MPIX_Dims_ml_create(info.nprocs, 2, MPIX_WEIGHTS_EQUAL, info.topo.size(),
                                       &info.topo[0], &dims_ml[0]))
          != MPI_SUCCESS) {
        throw mpi_error(ecode);
      }
      if (info.verbose && info.rank == 0) {
        std::string out;
        out = "x: ";
        for (size_t l = 0; l < info.topo.size(); ++l) {
          out += fmt::format("{} ", dims_ml[1 * info.topo.size() + l]);
        }
        out += "\ny: ";
        for (size_t l = 0; l < info.topo.size(); ++l) {
          out += fmt::format("{} ", dims_ml[0 * info.topo.size() + l]);
        }
        out += "\n";
        fmt::print("{}", out);
      }
      MPI_Comm comm_2d;
      MPIX_Cart_ml_create(MPI_COMM_WORLD, 2, periods, info.topo.size(), &dims_ml[0], MPI_INFO_NULL,
                          dims, &comm_2d);
      info.comm_2d = Comm{comm_2d};
      info.xnp = dims[1];
      info.ynp = dims[0];
    }

    MPI_Comm_rank(info.comm_2d.get(), &info.rank2d);

    int coords[2];
    MPI_Cart_coords(info.comm_2d.get(), info.rank2d, 2, coords);
    MPI_Cart_shift(info.comm_2d.get(), 0, 1, &info.rank_up, &info.rank_down);
    MPI_Cart_shift(info.comm_2d.get(), 1, 1, &info.rank_left, &info.rank_right);
    info.my_grid_x = coords[1];
    info.my_grid_y = coords[0];
    info.chunk_size_x = info.L / info.xnp;
    if (info.L % info.xnp) {
      throw std::invalid_argument(fmt::format("L ({}) mod xnp ({}) must be 0", info.L, info.xnp));
    }
    info.chunk_size_y = info.L / info.ynp;
    if (info.L % info.ynp) {
      throw std::invalid_argument(fmt::format("L ({}) mod ynp ({}) must be 0", info.L, info.ynp));
    }
    info.my_begin_x = info.chunk_size_x * info.my_grid_x;
    info.my_begin_y = info.chunk_size_y * info.my_grid_y;
    info.my_end_x = info.my_begin_x + info.chunk_size_x;
    info.my_end_y = info.my_begin_y + info.chunk_size_y;

    MPI_Datatype t;
    int array_shape[] = {info.chunk_size_y + 2, info.chunk_size_x + 2};
    int chunk_shape[] = {info.chunk_size_y, info.chunk_size_x};
    int chunk_start[] = {1, 1};
    MPI_Type_create_subarray(2, array_shape, chunk_shape, chunk_start, MPI_ORDER_C, MPI_DOUBLE, &t);
    MPI_Type_commit(&t);
    info.memtype = Datatype{t};

    MPI_Type_vector(info.chunk_size_y, 1, info.chunk_size_x + 2, MPI_DOUBLE, &t);
    MPI_Type_commit(&t);
    info.vertical_halo_type = Datatype{t};

    if (info.iotype == IOType::View) {
      array_shape[0] = array_shape[1] = info.L;
      chunk_shape[0] = info.chunk_size_y;
      chunk_shape[1] = info.chunk_size_x;
      chunk_start[0] = info.my_begin_y;
      chunk_start[1] = info.my_begin_x;
      MPI_Type_create_subarray(2, array_shape, chunk_shape, chunk_start, MPI_ORDER_C, MPI_DOUBLE,
                               &t);
      MPI_Type_commit(&t);
      info.filetype = Datatype{t};
    }

    return info;
  }

  std::string output_file(const int idx) const {
    return fmt::format("{}{}x{}-{}x{}-{:06}.bin", output_prefix, L, L, chunk_size_x, chunk_size_y,
                       idx);
  }

  size_t chunk_idx(const int iy, const int ix) const {
    return (1 + iy) * (2 + chunk_size_x) + 1 + ix;
  }

  bool is_inside(const int giy, const int gix) const {
    if (gix < my_begin_x) return false;
    if (gix >= my_end_x) return false;
    if (giy < my_begin_y) return false;
    if (giy >= my_end_y) return false;
    return true;
  }

  size_t g2c_idx(int giy, int gix) const {
    const int ix = gix - my_begin_x;
    const int iy = giy - my_begin_y;
    return chunk_idx(iy, ix);
  }
};

bmp::Bitmap crop_bitmap(const bmp::Bitmap &src, int max_width, int max_height) {
  // Calculate the new dimensions by centering the crop area
  int new_width = std::min(src.width(), max_width);
  int new_height = std::min(src.height(), max_height);
  int x_start = (src.width() - new_width) / 2;
  int y_start = (src.height() - new_height) / 2;

  // Create a new bitmap with the desired dimensions
  bmp::Bitmap cropped(new_width, new_height);

  // Copy pixels from the source image to the cropped image
  for (int y = 0; y < new_height; ++y) {
    for (int x = 0; x < new_width; ++x) {
      int src_x = x_start + x;
      int src_y = y_start + y;
      cropped.set(x, y, src.get(src_x, src_y));
    }
  }

  return cropped;
}

std::vector<double> convert_to_grayscale(const bmp::Bitmap &image, int &x_size, int &y_size) {
  // Dimensions
  x_size = image.width();
  y_size = image.height();

  // Vector to hold grayscale values
  std::vector<double> grayscale_values;
  grayscale_values.reserve(x_size * y_size);

  // Convert each pixel to grayscale
  for (int y = 0; y < y_size; ++y) {
    for (int x = 0; x < x_size; ++x) {
      const auto &pixel = image.get(x, y);
      double gray = 0.3 * pixel.r + 0.59 * pixel.g + 0.11 * pixel.b;
      grayscale_values.push_back(1. - gray / 255.);
    }
  }

  return grayscale_values;
}

void print_vd(const vd &v, int width, int height) {
  if (v.size() != width * height) {
    std::cerr << "Error: The size of the vector does not match the specified dimensions."
              << std::endl;
    return;
  }

  std::cout << std::fixed << std::setprecision(3);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      // Print each value in the vector, formatted to show 3 decimal places
      std::cout << v[y * width + x] << " ";
    }
    std::cout << std::endl;
  }
}

void init(vd &u, vd &v, RdbenchInfo &info) {
  int d = 3;
  const int L = info.L;

  if (info.u_bitmap_path) {
    int x_size, y_size;
    auto u_bitmap = crop_bitmap(bmp::Bitmap{*info.u_bitmap_path}, L, L);
    auto u_init = convert_to_grayscale(u_bitmap, x_size, y_size);
    if(info.rank == 0 && info.verbose) {
      std::cout << "initial u" << std::endl;
      print_vd(u_init, x_size, y_size);
    }
    int u_offset_y = L / 2 - y_size / 2;
    int u_offset_x = L / 2 - x_size / 2;
    for (int giy = u_offset_y; giy < L / 2 + y_size / 2; giy++) {
      for (int gix = u_offset_x; gix < L / 2 + x_size / 2; gix++) {
        if (!info.is_inside(giy, gix)) continue;
        // u[info.g2c_idx(giy, gix)] = u_init[(giy - u_offset_y) * x_size + (gix - u_offset_x)];
        u[info.g2c_idx(giy, gix)] = u_init[(giy - u_offset_y) * x_size + (gix - u_offset_x)] > 0 ? 0.7 : 0;
      }
    }
  } else {
    for (int giy = L / 2 - d; giy < L / 2 + d; giy++) {
      for (int gix = L / 2 - d; gix < L / 2 + d; gix++) {
        if (!info.is_inside(giy, gix)) continue;
        u[info.g2c_idx(giy, gix)] = 0.7;
      }
    }
  }
  if (info.v_bitmap_path) {
    int x_size, y_size;
    auto v_bitmap = crop_bitmap(bmp::Bitmap{*info.v_bitmap_path}, L, L);
    auto v_init = convert_to_grayscale(v_bitmap, x_size, y_size);
    if(info.rank == 0 && info.verbose) {
      std::cout << "initial v" << std::endl;
      print_vd(v_init, x_size, y_size);
    }
    int v_offset_y = L / 2 - y_size / 2;
    int v_offset_x = L / 2 - x_size / 2;
    for (int giy = v_offset_y; giy < L / 2 + y_size / 2; giy++) {
      for (int gix = v_offset_x; gix < L / 2 + x_size / 2; gix++) {
        if (!info.is_inside(giy, gix)) continue;
        v[info.g2c_idx(giy, gix)] = v_init[(giy - v_offset_y) * x_size + (gix - v_offset_x)] > 0 ? 0.9 : 0;
      }
    }
  } else {
    d = 6;
    for (int giy = L / 2 - d; giy < L / 2 + d; giy++) {
      for (int gix = L / 2 - d; gix < L / 2 + d; gix++) {
        if (!info.is_inside(giy, gix)) continue;
        v[info.g2c_idx(giy, gix)] = 0.9;
      }
    }
  }
}

double calcU(double tu, double tv) { return tu * tu * tv - (F + k) * tu; }

double calcV(double tu, double tv) { return -tu * tu * tv + F * (1.0 - tv); }

double laplacian(int ix, int iy, vd &s, RdbenchInfo &info) {
  double ts = 0.0;
  const int l = info.chunk_size_x + 2;
  ts += s[ix - 1 + iy * l];
  ts += s[ix + 1 + iy * l];
  ts += s[ix + (iy - 1) * l];
  ts += s[ix + (iy + 1) * l];
  ts -= 4.0 * s[ix + iy * l];
  return ts;
}

void calc(vd &u, vd &v, vd &u2, vd &v2, RdbenchInfo &info) {
  const int lx = info.chunk_size_x + 2;
  const int ly = info.chunk_size_y + 2;
  for (int iy = 1; iy < ly - 1; iy++) {
    for (int ix = 1; ix < lx - 1; ix++) {
      double du = 0;
      double dv = 0;
      const int i = ix + iy * lx;
      du = Du * laplacian(ix, iy, u, info);
      dv = Dv * laplacian(ix, iy, v, info);
      du += calcU(u[i], v[i]);
      dv += calcV(u[i], v[i]);
      u2[i] = u[i] + du * dt;
      v2[i] = v[i] + dv * dt;
    }
  }
}

void sendrecv_halo(vd &local_data, RdbenchInfo &info) {
  const static int iup = 0;
  const static int iright = 1;
  const static int idown = 2;
  const static int ileft = 3;
  const static int tag[] = {100, 101, 102, 103};
  MPI_Request req[4];
  MPI_Status status[4];
  MPI_Comm comm_2d = info.comm_2d.get();
  // up -> down
  MPI_Irecv(&local_data[info.chunk_idx(-1, 0)], info.chunk_size_x, MPI_DOUBLE, info.rank_up,
            tag[iup], comm_2d, &req[iup]);
  // down -> up
  MPI_Irecv(&local_data[info.chunk_idx(info.chunk_size_y, 0)], info.chunk_size_x, MPI_DOUBLE,
            info.rank_down, tag[idown], comm_2d, &req[idown]);
  // left -> right
  MPI_Irecv(&local_data[info.chunk_idx(0, -1)], 1, info.vertical_halo_type.get(), info.rank_left,
            tag[ileft], comm_2d, &req[ileft]);
  // left <- right
  MPI_Irecv(&local_data[info.chunk_idx(0, info.chunk_size_x)], 1, info.vertical_halo_type.get(),
            info.rank_right, tag[iright], comm_2d, &req[iright]);

  // up -> down
  MPI_Send(&local_data[info.chunk_idx(info.chunk_size_y - 1, 0)], info.chunk_size_x, MPI_DOUBLE,
           info.rank_down, tag[iup], comm_2d);
  // down -> up
  MPI_Send(&local_data[info.chunk_idx(0, 0)], info.chunk_size_x, MPI_DOUBLE, info.rank_up,
           tag[idown], comm_2d);
  // left -> right
  MPI_Send(&local_data[info.chunk_idx(0, info.chunk_size_x - 1)], 1, info.vertical_halo_type.get(),
           info.rank_right, tag[ileft], comm_2d);
  // left <- right
  MPI_Send(&local_data[info.chunk_idx(0, 0)], 1, info.vertical_halo_type.get(), info.rank_left,
           tag[iright], comm_2d);

  for (int i = 0; i < 4; ++i) {
    MPI_Wait(&req[i], &status[i]);
  }
}

void write_file(vd &local_data, int index, RdbenchInfo &info) {
  std::string filename = info.output_file(index);

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(info.comm_2d.get(), filename.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh);
  MPI_File_set_atomicity(fh, false);

  int (*write_func)(MPI_File, MPI_Offset, const void *, int, MPI_Datatype, MPI_Status *);
  write_func = info.collective ? MPI_File_write_at_all : MPI_File_write_at;

  if (info.iotype == IOType::View) {
    MPI_File_set_view(fh, 0, MPI_DOUBLE, info.filetype.get(), "native", MPI_INFO_NULL);
    int wcount;
    do {
      write_func(fh, 0, local_data.data(), 1, info.memtype.get(), &status);
      MPI_Get_count(&status, info.memtype.get(), &wcount);
    } while (wcount != 1);
  } else if (info.iotype == IOType::Chunk) {
    int count = info.chunk_size_x * info.chunk_size_y;
    int wcount;
    do {
      write_func(fh, sizeof(vd::value_type) * count * info.rank2d, local_data.data(), 1,
                 info.memtype.get(), &status);
      MPI_Get_count(&status, info.memtype.get(), &wcount);
    } while (wcount != 1);
  } else {
    int wcount, wc;
    for (int iy = 0; iy < info.chunk_size_y; ++iy) {
      int i_begin_local = (iy + 1) * (info.chunk_size_x + 2) + 1;
      int i_begin_global = (iy + info.my_begin_y) * info.L + info.my_begin_x;
      wcount = 0;
      do {
        write_func(fh, (i_begin_global + wcount) * sizeof(double),
                   &local_data[i_begin_local + wcount], info.chunk_size_x - wcount, MPI_DOUBLE,
                   &status);
        MPI_Get_count(&status, MPI_DOUBLE, &wc);
        wcount += wc;
      } while ((info.chunk_size_x - wcount) > 0);
    }
  }

  if (info.sync) {
    MPI_File_sync(fh);
  }

  MPI_File_close(&fh);
}

void read_file(vd &local_data, int index, RdbenchInfo &info) {
  std::string filename = info.output_file(index);

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(info.comm_2d.get(), filename.c_str(), MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN,
                MPI_INFO_NULL, &fh);
  MPI_File_set_atomicity(fh, false);

  int (*read_func)(MPI_File, MPI_Offset, void *, int, MPI_Datatype, MPI_Status *);
  read_func = info.collective ? MPI_File_read_at_all : MPI_File_read_at;

  if (info.iotype == IOType::View) {
    MPI_File_set_view(fh, 0, MPI_DOUBLE, info.filetype.get(), "native", MPI_INFO_NULL);
    int rcount;
    do {
      read_func(fh, 0, &local_data[0], 1, info.memtype.get(), &status);
      MPI_Get_count(&status, info.memtype.get(), &rcount);
    } while (rcount != 1);
  } else if (info.iotype == IOType::Chunk) {
    int count = info.chunk_size_x * info.chunk_size_y;
    int wcount;
    do {
      read_func(fh, sizeof(vd::value_type) * count * info.rank2d, &local_data[0], 1,
                info.memtype.get(), &status);
      MPI_Get_count(&status, info.memtype.get(), &wcount);
    } while (wcount != 1);
  } else {
    int rcount, rc;
    for (int iy = 0; iy < info.chunk_size_y; ++iy) {
      int i_begin_local = (iy + 1) * (info.chunk_size_x + 2) + 1;
      int i_begin_global = (iy + info.my_begin_y) * info.L + info.my_begin_x;
      rcount = 0;
      do {
        read_func(fh, (i_begin_global + rcount) * sizeof(double),
                  &local_data[i_begin_local + rcount], info.chunk_size_x - rcount, MPI_DOUBLE,
                  &status);
        MPI_Get_count(&status, MPI_DOUBLE, &rc);
        rcount += rc;
      } while ((info.chunk_size_x - rcount) > 0);
    }
  }

  MPI_File_close(&fh);
}

bool validate_file_io(vd &u, int index, RdbenchInfo &info) {
  vd validate(u.size(), 0.0);
  write_file(u, index, info);
  read_file(validate, index, info);

  MPI_Barrier(MPI_COMM_WORLD);
  if (info.rank == 0) {
    MPI_File_delete(info.output_file(index).c_str(), MPI_INFO_NULL);
  }

  const size_t stride = info.chunk_size_x + 2;
  auto iu = u.begin() + 1 + stride;
  auto iv = validate.begin() + 1 + stride;

  for (size_t i; i < info.chunk_size_y; ++i) {
    if (!std::equal(iu, iu + info.chunk_size_x, iv, [](double a, double b) {
          return std::fabs(a - b) <= std::numeric_limits<double>::epsilon()
                 || std::fabs(a - b) < std::numeric_limits<double>::min();
        })) {
      return false;
    }
    iu += stride;
    iv += stride;
  }

  return true;
}

orderd_json calc_result(
    const std::vector<std::pair<RdbenchPhase, Stopwatch::duration>> &phase_durations,
    RdbenchInfo &info) {
  size_t nfiles = info.write_phase_count;
  size_t file_size = static_cast<size_t>(info.L) * static_cast<size_t>(info.L) * sizeof(double);
  size_t total_write_size = nfiles * file_size;

  std::vector<std::pair<RdbenchPhase, int64_t>> phase_durations_i8;
  std::transform(phase_durations.begin(), phase_durations.end(),
                 std::back_inserter(phase_durations_i8),
                 [](const auto &p) { return std::make_pair(p.first, p.second.count()); });
  int64_t tc
      = std::accumulate(phase_durations_i8.begin(), phase_durations_i8.end(), 0LL,
                        [](int64_t acc, const decltype(phase_durations_i8)::value_type &cur) {
                          if (cur.first == RdbenchPhase::Calc) {
                            acc += cur.second;
                          }
                          return acc;
                        });

  int64_t tcomm
      = std::accumulate(phase_durations_i8.begin(), phase_durations_i8.end(), 0LL,
                        [](int64_t acc, const decltype(phase_durations_i8)::value_type &cur) {
                          if (cur.first == RdbenchPhase::Comm) {
                            acc += cur.second;
                          }
                          return acc;
                        });

  int64_t tw
      = std::accumulate(phase_durations_i8.begin(), phase_durations_i8.end(), 0LL,
                        [](int64_t acc, const decltype(phase_durations_i8)::value_type &cur) {
                          if (cur.first == RdbenchPhase::Write) {
                            acc += cur.second;
                          }
                          return acc;
                        });

  int64_t max_tc, max_tw, max_tcomm;
  MPI_Reduce(&tc, &max_tc, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&tcomm, &max_tcomm, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&tw, &max_tw, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);

  std::vector<int64_t> durations_i8;
  std::transform(phase_durations_i8.begin(), phase_durations_i8.end(),
                 std::back_inserter(durations_i8), [](const auto &p) { return p.second; });
  std::vector<int64_t> durations_all_i8(phase_durations_i8.size() * info.nprocs);
  MPI_Gather(durations_i8.data(), durations_i8.size(), MPI_INT64_T, durations_all_i8.data(),
             durations_i8.size(), MPI_INT64_T, 0, MPI_COMM_WORLD);

  if (info.rank != 0) {
    return nullptr;
  }

  auto to_sec = [](int64_t t) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(Stopwatch::duration{t})
        .count();
  };

  double calc_time_sec = to_sec(max_tc);
  double comm_time_sec = to_sec(max_tcomm);
  double write_time_sec = to_sec(max_tw);

  orderd_json rdbench_result = {
      {"version", RDBENCH_VERSION},
      {"nprocs", info.nprocs},
      {"topology", info.topo},
      {"xnp", info.xnp},
      {"ynp", info.ynp},
      {"L", info.L},
      {"chunkSizeX", info.chunk_size_x},
      {"chunkSizeY", info.chunk_size_y},
      {"collective", info.collective},
      {"iotype", info.iot},
      {"sync", info.sync},
      {"validate", info.validate},
      {"steps", info.total_steps},
      {"interval", info.interval},
      {"fixedX", info.fixed_x},
      {"fixedY", info.fixed_y},
      {"initialOutput", info.initial_output},
      {"nfiles", nfiles},
      {"fileSize", file_size},
      {"totalWriteSizeByte", total_write_size},
      {"calcTimeSec", calc_time_sec},
      {"commTimeSec", comm_time_sec},
  };

  if (nfiles > 0) {
    rdbench_result["writeTimeSec"] = write_time_sec;
    rdbench_result["writeBandwidthByte"]
        = std::stod(fmt::format("{:.2f}", total_write_size / write_time_sec));
  }

  auto phase_durations_json = nlohmann::ordered_json::array();
  for (size_t i = 0; i < phase_durations.size(); ++i) {
    double max = 0.0;
    for (size_t rank = 0; rank < info.nprocs; ++rank) {
      double cur = to_sec(durations_all_i8[rank * phase_durations.size() + i]);
      if (cur > max) {
        max = cur;
      }
    }

    phase_durations_json.push_back(orderd_json{phase_to_string[phase_durations[i].first], max});
  }
  rdbench_result["initialTimeSec"] = phase_durations_json.front()[1];
  rdbench_result["finalizeTimeSec"] = nullptr;
  rdbench_result["phaseDurationsSec"] = phase_durations_json;

  return rdbench_result;
}

void print_cartesian(RdbenchInfo &info) {
  int rank2d;
  MPI_Comm_rank(info.comm_2d.get(), &rank2d);
  std::vector<int> ranks(info.nprocs);
  MPI_Gather(&info.rank, 1, MPI_INT, &ranks[0], 1, MPI_INT, 0, info.comm_2d.get());

  if (info.rank == 0) {
    for (int y = 0; y < info.ynp; ++y) {
      for (int x = 0; x < info.xnp; ++x) {
        fmt::print("{:>3} ", ranks[y * info.xnp + x]);
      }
      fmt::print("\n");
    }
  }
}

void ensure_output_directory_exists(RdbenchInfo &info) {
  namespace fs = std::filesystem;
  fs::path out = info.output_prefix;
  if (out.has_parent_path()) {
    auto parent = std::regex_replace(out.parent_path().string(), std::regex("^[^:]+:"), "");
    fs::create_directories(parent);
  }
}

int main(int argc, char *argv[]) {
  cxxopts::Options options("rdbench", "MPI/MPI-IO benchmark based on 2D reaction-diffusion system");
  // clang-format off
  options.set_width(100).set_tab_expansion().add_options()
    ("h,help", "Print usage")
    ("V,version", "Print version")
    ("v,verbose", "Verbose output")
    ("t,topology", "Variable level hardware topology (e.g. -t 2,8 2-node 8-core, -t 4,2,4 4-node 2-socket 4-core) default: np (flat topology)", cxxopts::value<std::vector<int>>())
    ("xnp", "Number of processes in x-axis (0 == auto), is ignored if -t is set.", cxxopts::value<int>()->default_value("0"))
    ("ynp", "Number of processes in y-axis (0 == auto), is ignored if -t is set", cxxopts::value<int>()->default_value("0"))
    ("L,length", "Length of a edge of a square region", cxxopts::value<int>()->default_value("128"))
    ("init_v_from_bitmap", "Initialize v from a bitmap file", cxxopts::value<std::string>())
    ("init_u_from_bitmap", "Initialize u from a bitmap file", cxxopts::value<std::string>())
    ("o,output", "Prefix of output files", cxxopts::value<std::string>()->default_value("out/o_"))
    ("iotype",
      "manual / view / chunk\n"
      "\tmanual: 2d block-cyclic, write line-by-line.\n"
      "\tview: 2d block-cyclic, use MPI_File_set_view.\n"
      "\tchunk: block partitioning.\n", cxxopts::value<std::string>()->default_value("view"))
    ("c,collective", "Writing files in collective mode of MPI-IO")
    ("nosync", "MPI_File_sync is no longer called.")
    ("s,steps", "Total steps", cxxopts::value<size_t>()->default_value("20000"))
    ("i,interval", "Write the array into files every i steps (0 == disable file output)", cxxopts::value<size_t>()->default_value("200"))
    ("fixed-x", "Fixed boundary in x-axis")
    ("fixed-y", "Fixed boundary in y-axis")
    ("novalidate", "Disable IO validation feature reading the data written in the file to check if it was written correctly")
    ("disable-initial-output", "Disable file output for initial state")
    ("nomkdir", "Disable output directory craetion using POSIX")
    ("prettify", "Prettify JSON output")
    ("param_F", "Parameter F of the Gray-Scott model", cxxopts::value<double>()->default_value("0.04"))
    ("param_k", "Parameter k of the Gray-Scott model", cxxopts::value<double>()->default_value("0.06075"))
    ("param_dt", "Parameter dt of the Gray-Scott model", cxxopts::value<double>()->default_value("0.2"))
    ("param_Du", "Parameter Du of the Gray-Scott model", cxxopts::value<double>()->default_value("0.05"))
    ("param_Dv", "Parameter Dv of the Gray-Scott model", cxxopts::value<double>()->default_value("0.1"))
  ;
  // clang-format on

  auto parsed = options.parse(argc, argv);
  if (parsed.count("help") != 0U) {
    fmt::print("{}\n", options.help());
    return 0;
  }

  if (parsed.count("version") != 0U) {
    fmt::print("{}\n", RDBENCH_VERSION);
    return 0;
  }

  int ret = 0;
  Stopwatch stopwatch;
  std::vector<std::pair<RdbenchPhase, Stopwatch::duration>> phase_durations;
  orderd_json rdbench_result;
  try {
    stopwatch.reset();
    MPI_Init(&argc, &argv);
    // default error handler for MPI-IO
    MPI_File_set_errhandler(MPI_FILE_NULL, MPI_ERRORS_ARE_FATAL);
    RdbenchInfo info = RdbenchInfo::create(parsed);

    if (info.verbose) {
      print_cartesian(info);
    }

    if (info.create_output_dir && info.rank == 0) {
      ensure_output_directory_exists(info);
    }

    const int V = (info.chunk_size_x + 2) * (info.chunk_size_y + 2);
    vd u(V, 0.0), v(V, 0.0);
    vd u2(V, 0.0), v2(V, 0.0);
    int step;
    int file_idx = 0;
    init(u, v, info);

    if (info.validate) {
      if (!validate_file_io(u, file_idx++, info)) {
        throw std::runtime_error("Read validation failed");
      }
    }

    phase_durations.reserve(1 + info.calc_phase_count + info.write_phase_count + 1);

    phase_durations.emplace_back(RdbenchPhase::Init, stopwatch.get_and_reset());
    if (info.initial_output && info.interval != 0) {
      write_file(u, file_idx++, info);
      phase_durations.emplace_back(RdbenchPhase::Write, stopwatch.get_and_reset());
    }

    Stopwatch::duration acc_comm{};
    Stopwatch::duration acc_calc{};
    for (step = 1; step <= info.total_steps; step++) {
      if (step & 1) {
        sendrecv_halo(u, info);
        sendrecv_halo(v, info);
        acc_comm += stopwatch.get_and_reset();
        calc(u, v, u2, v2, info);
        acc_calc += stopwatch.get_and_reset();
      } else {
        sendrecv_halo(u2, info);
        sendrecv_halo(v2, info);
        acc_comm += stopwatch.get_and_reset();
        calc(u2, v2, u, v, info);
        acc_calc += stopwatch.get_and_reset();
      }
      if (info.interval != 0 && step % info.interval == 0) {
        phase_durations.emplace_back(RdbenchPhase::Comm, acc_comm);
        phase_durations.emplace_back(RdbenchPhase::Calc, acc_calc);
        acc_comm = acc_calc = Stopwatch::duration::zero();
        write_file(step & 1 ? u : u2, file_idx++, info);
        phase_durations.emplace_back(RdbenchPhase::Write, stopwatch.get_and_reset());
      }
    }
    if (info.interval == 0 || info.total_steps % info.interval != 0) {
      phase_durations.emplace_back(RdbenchPhase::Comm, acc_comm);
      phase_durations.emplace_back(RdbenchPhase::Calc, acc_calc);
    }

    if (info.validate) {
      if (!validate_file_io(u, file_idx++, info)) {
        throw std::runtime_error("Read validation failed");
      }
    }

    rdbench_result = calc_result(phase_durations, info);
  } catch (const std::exception &e) {
    fmt::print(stderr, "exception: {}\n", e.what());
    ret = -1;
  }

  MPI_Finalize();
  if (!rdbench_result.is_null()) {
    double finalize_time_sec
        = std::chrono::duration_cast<std::chrono::duration<double>>(stopwatch.get_and_reset())
              .count();
    rdbench_result["finalizeTimeSec"] = finalize_time_sec;
    rdbench_result["phaseDurationsSec"].push_back(
        orderd_json{phase_to_string[RdbenchPhase::Finalize], finalize_time_sec});
    if (parsed.count("prettify") != 0U) {
      std::cout << std::setw(4);
    }
    std::cout << rdbench_result << std::endl;
  }
  return ret;
}
