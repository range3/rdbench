#include <fmt/core.h>
#include <mpi.h>

#include <chrono>
#include <cstdio>
#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "stopwatch.hpp"

using json = nlohmann::json;

const double F = 0.04;
const double k = 0.06075;
const double dt = 0.2;
const double Du = 0.05;
const double Dv = 0.1;

typedef std::vector<double> vd;

struct RdbenchInfo {
  int rank;
  int nprocs;
  int xnp;
  int ynp;
  int L;
  int my_grid_x;
  int my_grid_y;
  int chunk_size_x;
  int chunk_size_y;
  int my_begin_x;
  int my_begin_y;
  int my_end_x;
  int my_end_y;
  std::string output_prefix;
  bool collective;
  bool view;
  bool nosync;
  MPI_Datatype filetype = MPI_DATATYPE_NULL;
  MPI_Datatype memtype = MPI_DATATYPE_NULL;
  MPI_Datatype vertical_halo_type = MPI_DATATYPE_NULL;
  size_t total_steps;
  size_t interval;
  MPI_Comm comm_2d = MPI_COMM_NULL;
  bool fixed_x = false;
  bool fixed_y = false;
  int rank_down;
  int rank_up;
  int rank_left;
  int rank_right;

  static RdbenchInfo create(cxxopts::ParseResult &parsed) {
    RdbenchInfo info;
    MPI_Comm_size(MPI_COMM_WORLD, &info.nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &info.rank);

    info.output_prefix = parsed["output"].as<std::string>();
    info.collective = parsed.count("collective") != 0U;
    info.view = parsed.count("view") != 0U;
    info.nosync = parsed.count("nosync") != 0U;
    int dims[2] = {parsed["ynp"].as<int>(), parsed["xnp"].as<int>()};
    MPI_Dims_create(info.nprocs, 2, dims);
    info.xnp = dims[1];
    info.ynp = dims[0];
    info.fixed_x = parsed["fixed-x"].count() != 0U;
    info.fixed_y = parsed["fixed-y"].count() != 0U;

    int periods[] = {info.fixed_y ? 0 : 1, info.fixed_x ? 0 : 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &info.comm_2d);
    // if (info.rank == 0) {
    //   int c[2];
    //   int rank;
    //   for (c[0] = 0; c[0] < dims[0]; ++c[0]) {
    //     for (c[1] = 0; c[1] < dims[1]; ++c[1]) {
    //       MPI_Cart_rank(info.comm_2d, c, &rank);
    //       fmt::print("({}, {}) = rank {}\n", c[0], c[1], rank);
    //     }
    //   }
    // }

    info.total_steps = parsed["steps"].as<size_t>();
    info.interval = parsed["interval"].as<size_t>();
    info.L = parsed["L"].as<int>();

    int coords[2];
    MPI_Cart_coords(info.comm_2d, info.rank, 2, coords);
    // fmt::print("{}: coords = ({}, {})\n", info.rank, coords[0], coords[1]);
    MPI_Cart_shift(info.comm_2d, 0, 1, &info.rank_up, &info.rank_down);
    MPI_Cart_shift(info.comm_2d, 1, 1, &info.rank_left, &info.rank_right);
    // fmt::print("{}: neighbor = ({}, {}, {}, {})\n", info.rank, info.rank_up, info.rank_right,
    //            info.rank_down, info.rank_left);
    info.my_grid_x = coords[1];
    info.my_grid_y = coords[0];
    info.chunk_size_x = info.L / info.xnp;
    if (info.L % info.xnp) {
      throw std::invalid_argument("L mod xnp must be 0");
    }
    info.chunk_size_y = info.L / info.ynp;
    if (info.L % info.ynp) {
      throw std::invalid_argument("L mod ynp must be 0");
    }
    info.my_begin_x = info.chunk_size_x * info.my_grid_x;
    info.my_begin_y = info.chunk_size_y * info.my_grid_y;
    info.my_end_x = info.my_begin_x + info.chunk_size_x;
    info.my_end_y = info.my_begin_y + info.chunk_size_y;

    if (info.view) {
      int array_shape[] = {info.L, info.L};
      int chunk_shape[] = {info.chunk_size_y, info.chunk_size_x};
      int chunk_start[] = {info.my_begin_y, info.my_begin_x};
      MPI_Type_create_subarray(2, array_shape, chunk_shape, chunk_start, MPI_ORDER_C, MPI_DOUBLE,
                               &info.filetype);
      MPI_Type_commit(&info.filetype);

      array_shape[0] = info.chunk_size_y + 2;
      array_shape[1] = info.chunk_size_x + 2;
      chunk_shape[0] = info.chunk_size_y;
      chunk_shape[1] = info.chunk_size_x;
      chunk_start[0] = chunk_start[1] = 1;
      MPI_Type_create_subarray(2, array_shape, chunk_shape, chunk_start, MPI_ORDER_C, MPI_DOUBLE,
                               &info.memtype);
      MPI_Type_commit(&info.memtype);
    }
    MPI_Type_vector(info.chunk_size_y, 1, info.chunk_size_x + 2, MPI_DOUBLE,
                    &info.vertical_halo_type);
    MPI_Type_commit(&info.vertical_halo_type);

    return info;
  }

  std::string output_file(const int idx) const {
    return fmt::format("{}{:06}.bin", output_prefix, idx);
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

void init(vd &u, vd &v, RdbenchInfo &info) {
  int d = 3;
  const int L = info.L;
  for (int giy = L / 2 - d; giy < L / 2 + d; giy++) {
    for (int gix = L / 2 - d; gix < L / 2 + d; gix++) {
      if (!info.is_inside(giy, gix)) continue;
      u[info.g2c_idx(giy, gix)] = 0.7;
    }
  }
  d = 6;
  for (int giy = L / 2 - d; giy < L / 2 + d; giy++) {
    for (int gix = L / 2 - d; gix < L / 2 + d; gix++) {
      if (!info.is_inside(giy, gix)) continue;
      v[info.g2c_idx(giy, gix)] = 0.9;
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
  // up -> down
  MPI_Irecv(&local_data[info.chunk_idx(-1, 0)], info.chunk_size_x, MPI_DOUBLE, info.rank_up,
            tag[iup], info.comm_2d, &req[iup]);
  // down -> up
  MPI_Irecv(&local_data[info.chunk_idx(info.chunk_size_y, 0)], info.chunk_size_x, MPI_DOUBLE,
            info.rank_down, tag[idown], info.comm_2d, &req[idown]);
  // left -> right
  MPI_Irecv(&local_data[info.chunk_idx(0, -1)], 1, info.vertical_halo_type, info.rank_left,
            tag[ileft], info.comm_2d, &req[ileft]);
  // left <- right
  MPI_Irecv(&local_data[info.chunk_idx(0, info.chunk_size_x)], 1, info.vertical_halo_type,
            info.rank_right, tag[iright], info.comm_2d, &req[iright]);

  // up -> down
  MPI_Send(&local_data[info.chunk_idx(info.chunk_size_y - 1, 0)], info.chunk_size_x, MPI_DOUBLE,
           info.rank_down, tag[iup], info.comm_2d);
  // down -> up
  MPI_Send(&local_data[info.chunk_idx(0, 0)], info.chunk_size_x, MPI_DOUBLE, info.rank_up,
           tag[idown], info.comm_2d);
  // left -> right
  MPI_Send(&local_data[info.chunk_idx(0, info.chunk_size_x - 1)], 1, info.vertical_halo_type,
           info.rank_right, tag[ileft], info.comm_2d);
  // left <- right
  MPI_Send(&local_data[info.chunk_idx(0, 0)], 1, info.vertical_halo_type, info.rank_left,
           tag[iright], info.comm_2d);

  for (int i = 0; i < 4; ++i) {
    MPI_Wait(&req[i], &status[i]);
  }
}

void write_file(vd &local_data, RdbenchInfo &info) {
  static int index = 0;
  std::string filename = info.output_file(index);

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh);
  MPI_File_set_atomicity(fh, false);

  if (info.view) {
    MPI_File_set_view(fh, 0, MPI_DOUBLE, info.filetype, "native", MPI_INFO_NULL);
  }

  int (*write_func)(MPI_File, MPI_Offset, const void *, int, MPI_Datatype, MPI_Status *);
  write_func = info.collective ? MPI_File_write_at_all : MPI_File_write_at;

  if (info.view) {
    int wcount;
    do {
      write_func(fh, 0, local_data.data(), 1, info.memtype, &status);
      MPI_Get_count(&status, info.memtype, &wcount);
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

  if (!info.nosync) {
    MPI_File_sync(fh);
  }

  MPI_File_close(&fh);
  index += 1;
}

void print_result(Stopwatch::duration calc_time, Stopwatch::duration write_time,
                  RdbenchInfo &info) {
  size_t nfiles = 1 + (info.interval > 0 ? info.total_steps / info.interval : 0);
  size_t file_size = info.L * info.L * sizeof(double);
  size_t total_write_size = nfiles * file_size;
  double calc_time_sec
      = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(calc_time).count();
  double write_time_sec
      = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(write_time)
            .count();

  json rdbench_result = {
      {"nprocs", info.nprocs},
      {"xnp", info.xnp},
      {"ynp", info.ynp},
      {"L", info.L},
      {"chunkSizeX", info.chunk_size_x},
      {"chunkSizeY", info.chunk_size_y},
      {"collective", info.collective},
      {"view", info.view},
      {"nosync", info.nosync},
      {"steps", info.total_steps},
      {"interval", info.interval},
      {"fixedX", info.fixed_x},
      {"fixedY", info.fixed_y},
      {"nfiles", nfiles},
      {"fileSize", file_size},
      {"totalWriteSizeByte", total_write_size},
      {"calcTimeSec", calc_time_sec},
      {"writeTimeSec", write_time_sec},
      {"writeBandwidthByte", std::stod(fmt::format("{:.2f}", total_write_size / write_time_sec))}};

  std::cout << rdbench_result << std::endl;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options("rdbench_int",
                           "MPI/MPI-IO benchmark based on 2d reaction-diffusion system");
  // clang-format off
  options.add_options()
    ("h,help", "Print usage")
    ("xnp", "Number of processes in x-axis (0 == auto)", cxxopts::value<int>()->default_value("0"))
    ("ynp", "Number of processes in y-axis (0 == auto)", cxxopts::value<int>()->default_value("0"))
    ("L,length", "Length of a edge of a square region", cxxopts::value<int>()->default_value("128"))
    ("o,output", "Prefix of output files", cxxopts::value<std::string>()->default_value("./out_"))
    ("c,collective", "writing files in collective mode of MPI-IO")
    ("v,view", "Use MPI_File_set_view")
    ("nosync", "MPI_File_sync is no longer called before closing each file.")
    ("s,steps", "Total steps", cxxopts::value<size_t>()->default_value("20000"))
    ("i,interval", "Write the array into files every i steps", cxxopts::value<size_t>()->default_value("200"))
    ("fixed-x", "Fixed boundary in x-axis")
    ("fixed-y", "Fixed boundary in y-axis")
  ;
  // clang-format on

  auto parsed = options.parse(argc, argv);
  if (parsed.count("help") != 0U) {
    fmt::print("{}\n", options.help());
    return 0;
  }

  int ret = 0;
  try {
    MPI_Init(&argc, &argv);
    // default error handler for MPI-IO
    MPI_File_set_errhandler(MPI_FILE_NULL, MPI_ERRORS_ARE_FATAL);
    RdbenchInfo info = RdbenchInfo::create(parsed);
    Stopwatch stopwatch;
    Stopwatch::duration time_calc(0), time_write(0);

    const int V = (info.chunk_size_x + 2) * (info.chunk_size_y + 2);
    vd u(V, 0.0), v(V, 0.0);
    vd u2(V, 0.0), v2(V, 0.0);
    init(u, v, info);

    stopwatch.reset();
    write_file(u, info);
    time_write += stopwatch.get_and_reset();

    for (int step = 1; step <= info.total_steps; step++) {
      if (step & 1) {
        sendrecv_halo(u, info);
        sendrecv_halo(v, info);
        calc(u, v, u2, v2, info);
      } else {
        sendrecv_halo(u2, info);
        sendrecv_halo(v2, info);
        calc(u2, v2, u, v, info);
      }
      time_calc += stopwatch.get_and_reset();
      if (info.interval != 0 && step % info.interval == 0) {
        write_file(step & 1 ? u : u2, info);
        time_write += stopwatch.get_and_reset();
      }
    }

    // FIXME: use RAII
    if (info.filetype != MPI_DATATYPE_NULL) {
      MPI_Type_free(&info.filetype);
    }
    if (info.memtype != MPI_DATATYPE_NULL) {
      MPI_Type_free(&info.memtype);
    }
    if (info.vertical_halo_type != MPI_DATATYPE_NULL) {
      MPI_Type_free(&info.vertical_halo_type);
    }
    if (info.comm_2d != MPI_COMM_NULL) {
      MPI_Comm_free(&info.comm_2d);
    }

    int64_t tc = time_calc.count();
    int64_t tw = time_write.count();
    int64_t max_tc, max_tw;
    MPI_Reduce(&tc, &max_tc, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tw, &max_tw, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);

    if (info.rank == 0) {
      print_result(Stopwatch::duration{max_tc}, Stopwatch::duration{max_tw}, info);
    }
  } catch (const std::exception &e) {
    fmt::print(stderr, "exception: {}\n", e.what());
    ret = -1;
  }

  MPI_Finalize();
  return ret;
}
