#include <fmt/core.h>
#include <mpi.h>

#include <cstdio>
#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <vector>

using json = nlohmann::json;

// const int L = 128;
// const int TOTAL_STEP = 20000;
// const int INTERVAL = 200;
// const double F = 0.04;
// const double k = 0.06075;
// const double dt = 0.2;
// const double Du = 0.05;
// const double Dv = 0.1;

// typedef std::vector<double> vd;
typedef std::vector<int> vd;

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
      MPI_Type_create_subarray(2, array_shape, chunk_shape, chunk_start, MPI_ORDER_C, MPI_INT,
                               &info.filetype);
      MPI_Type_commit(&info.filetype);

      array_shape[0] = info.chunk_size_y + 2;
      array_shape[1] = info.chunk_size_x + 2;
      chunk_shape[0] = info.chunk_size_y;
      chunk_shape[1] = info.chunk_size_x;
      chunk_start[0] = chunk_start[1] = 1;
      MPI_Type_create_subarray(2, array_shape, chunk_shape, chunk_start, MPI_ORDER_C, MPI_INT,
                               &info.memtype);
      MPI_Type_commit(&info.memtype);
    }
    MPI_Type_vector(info.chunk_size_y, 1, info.chunk_size_x + 2, MPI_INT, &info.vertical_halo_type);
    MPI_Type_commit(&info.vertical_halo_type);

    return info;
  }

  std::string output_file(int idx) const { return fmt::format("{}{:06}.bin", output_prefix, idx); }

  size_t chunk_idx(int iy, int ix) const { return (1 + iy) * (2 + chunk_size_x) + 1 + ix; }
};

// void init(vd &u, vd &v, RdbenchInfo& info) {
//   const int L = info.L;
//   int d = 3;
//   for (int i = L / 2 - d; i < L / 2 + d; i++) {
//     for (int j = L / 2 - d; j < L / 2 + d; j++) {
//       u[j + i * L] = 0.7;
//     }
//   }
//   d = 6;
//   for (int i = L / 2 - d; i < L / 2 + d; i++) {
//     for (int j = L / 2 - d; j < L / 2 + d; j++) {
//       v[j + i * L] = 0.9;
//     }
//   }
// }

// double calcU(double tu, double tv) { return tu * tu * tv - (F + k) * tu; }

// double calcV(double tu, double tv) { return -tu * tu * tv + F * (1.0 - tv); }

// double laplacian(int ix, int iy, vd &s) {
//   double ts = 0.0;
//   ts += s[ix - 1 + iy * L];
//   ts += s[ix + 1 + iy * L];
//   ts += s[ix + (iy - 1) * L];
//   ts += s[ix + (iy + 1) * L];
//   ts -= 4.0 * s[ix + iy * L];
//   return ts;
// }

// void calc(vd &u, vd &v, vd &u2, vd &v2) {
//   for (int iy = 1; iy < L - 1; iy++) {
//     for (int ix = 1; ix < L - 1; ix++) {
//       double du = 0;
//       double dv = 0;
//       const int i = ix + iy * L;
//       du = Du * laplacian(ix, iy, u);
//       dv = Dv * laplacian(ix, iy, v);
//       du += calcU(u[i], v[i]);
//       dv += calcV(u[i], v[i]);
//       u2[i] = u[i] + du * dt;
//       v2[i] = v[i] + dv * dt;
//     }
//   }
// }

void init(std::vector<int> &local_data, RdbenchInfo &info) {
  const int offset = info.chunk_size_x * info.chunk_size_y * info.rank;
  for (int iy = 0; iy < info.chunk_size_y; iy++) {
    for (int ix = 0; ix < info.chunk_size_x; ix++) {
      int index = (ix + 1) + (iy + 1) * (info.chunk_size_x + 2);
      int value = ix + iy * info.chunk_size_x + offset;
      local_data[index] = value;
    }
  }
}

void dump_local_sub(std::vector<int> &local_data, RdbenchInfo &info) {
  printf("rank = %d\n", info.rank);
  for (int iy = 0; iy < info.chunk_size_y + 2; iy++) {
    for (int ix = 0; ix < info.chunk_size_x + 2; ix++) {
      unsigned int index = ix + iy * (info.chunk_size_x + 2);
      printf("%03d ", local_data[index]);
    }
    printf("\n");
  }
  printf("\n");
}

void dump_local(std::vector<int> &local_data, RdbenchInfo &info) {
  for (int i = 0; i < info.nprocs; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (i == info.rank) {
      dump_local_sub(local_data, info);
    }
  }
}

void write_file(std::vector<int> &local_data, RdbenchInfo &info) {
  static int index = 0;
  std::string filename = info.output_file(index);

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh);
  MPI_File_set_atomicity(fh, false);

  if (info.view) {
    MPI_File_set_view(fh, 0, MPI_INT, info.filetype, "native", MPI_INFO_NULL);
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
        write_func(fh, (i_begin_global + wcount) * sizeof(int), &local_data[i_begin_local + wcount],
                   info.chunk_size_x - wcount, MPI_INT, &status);
        MPI_Get_count(&status, MPI_INT, &wc);
        wcount += wc;
      } while ((info.chunk_size_x - wcount) > 0);
    }
  }

  MPI_File_close(&fh);
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
  MPI_Irecv(&local_data[info.chunk_idx(-1, 0)], info.chunk_size_x, MPI_INT, info.rank_up, tag[iup],
            info.comm_2d, &req[iup]);
  // down -> up
  MPI_Irecv(&local_data[info.chunk_idx(info.chunk_size_y, 0)], info.chunk_size_x, MPI_INT,
            info.rank_down, tag[idown], info.comm_2d, &req[idown]);
  // left -> right
  MPI_Irecv(&local_data[info.chunk_idx(0, -1)], 1, info.vertical_halo_type, info.rank_left,
            tag[ileft], info.comm_2d, &req[ileft]);
  // left <- right
  MPI_Irecv(&local_data[info.chunk_idx(0, info.chunk_size_x)], 1, info.vertical_halo_type,
            info.rank_right, tag[iright], info.comm_2d, &req[iright]);

  // up -> down
  MPI_Send(&local_data[info.chunk_idx(info.chunk_size_y - 1, 0)], info.chunk_size_x, MPI_INT,
           info.rank_down, tag[iup], info.comm_2d);
  // down -> up
  MPI_Send(&local_data[info.chunk_idx(0, 0)], info.chunk_size_x, MPI_INT, info.rank_up, tag[idown],
           info.comm_2d);
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

int main(int argc, char *argv[]) {
  int ret = 0;
  try {
    MPI_Init(&argc, &argv);

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
      ("s,steps", "Total steps", cxxopts::value<size_t>()->default_value("20000"))
      ("i,interval", "Write the array into files every i steps", cxxopts::value<size_t>()->default_value("200"))
      ("fixed-x", "Fixed boundary in x-axis")
      ("fixed-y", "Fixed boundary in y-axis")
    ;
    // clang-format on

    auto parsed = options.parse(argc, argv);
    if (parsed.count("help") != 0U) {
      fmt::print("{}\n", options.help());
      MPI_Finalize();
      return ret;
    }

    RdbenchInfo info = RdbenchInfo::create(parsed);
    json rdbench_result = {
        {"nprocs", info.nprocs},
        {"xnp", info.xnp},
        {"ynp", info.ynp},
        {"L", info.L},
        {"chunk_size_x", info.chunk_size_x},
        {"chunk_size_y", info.chunk_size_y},
        {"collective", info.collective},
        {"view", info.view},
        {"steps", info.total_steps},
        {"interval", info.interval},
        {"fixedX", info.fixed_x},
        {"fixedY", info.fixed_y},
    };

    std::vector<int> local_data((info.chunk_size_x + 2) * (info.chunk_size_y + 2), 0);
    init(local_data, info);
    sendrecv_halo(local_data, info);
    // if (info.rank == 0) {
    //   fmt::print("--- init ---\n");
    // }
    // dump_local(local_data, info);

    if (info.rank == 0) {
      fmt::print("--- after exchange halo ---\n");
    }
    dump_local(local_data, info);

    // if (info.rank == 0) {
    //   fmt::print("L: {} ({}, {})\n", info.L, info.xnp, info.ynp);
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // for (int i = 0; i < info.nprocs; ++i) {
    //   if (i == info.rank) {
    //     fmt::print("rank: {}, ({}, {}), [{}:{}, {}:{}]\n", info.rank, info.my_grid_x,
    //                info.my_grid_y, info.my_begin_x, info.my_end_x, info.my_begin_y,
    //                info.my_end_y);
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);
    // }

    // write_file(local_data, info);

    if (info.rank == 0) {
      std::cout << rdbench_result << std::endl;
    }

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

  } catch (const std::exception &e) {
    fmt::print(stderr, "exception: {}\n", e.what());
    ret = -1;
  }

  // const int V = L * L;
  // vd u(V, 0.0), v(V, 0.0);
  // vd u2(V, 0.0), v2(V, 0.0);
  // init(u, v);
  // for (int i = 0; i < TOTAL_STEP; i++) {
  //   if (i & 1) {
  //     calc(u2, v2, u, v);
  //   } else {
  //     calc(u, v, u2, v2);
  //   }
  //   if (i % INTERVAL == 0) save_as_dat(u);
  // }
  MPI_Finalize();
  return ret;
}
