# rdbench v2
rdbench is a benchmark developed to evaluate file I/O performance in typical HPC simulation applications.
It simulates the time evolution of a 2D reaction-diffusion system based on the Gray-Scott model, using a five-point stencil computation with distributed parallel execution via MPI. Additionally, it supports periodic checkpointing using MPI-IO.

Furthermore, if a compiler supports ISO C++ `std::execution::par_unseq`, rdbench enables hybrid MPI computation by combining thread-level parallelism with MPI.
For environments where the NVIDIA HPC SDK (nvc++) is available, it also supports offloading computations to GPUs.





![](https://raw.githubusercontent.com/range3/rdbench/master/rdbench-viz/viz.gif)

# Building
## Requirements
- CMake
- C++20 compiler
  - gcc 11 or later
  - nvc++ (optional) for
    - RDBENCH_USE_GPU=ON
    - RDBENCH_USE_MULTICORE=ON
- dependencies
  - fmt
  - cxxmpi
  - cxxopts
  - nlohmann_json
  - kokkos/mdspan
  - MPI
    - OpenMPI
    - Nvidia HPC SDK
 - CUDATookkit for nvtx3 (optional) 

```bash
# basic build
# requires all dependencies to be installed (can find by find_package)
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build/

# use fetch content and vcpkg to install dependencies automatically
## install vcpkg somewhere
git clone https://github.com/microsoft/vcpkg.git
VCPKG_ROOT=$(pwd)/vcpkg
${VCPKG_ROOT}/bootstrap-vcpkg.sh

# basic flat MPI build
cmake -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" \
  -D CMAKE_MODULE_PATH=cmake/fetch_content
cmake --build build/

# (optional) optimal build with Unified Memory GPU support
# (nvc++ -stdpar=gpu -gpu=cc90,mem:unified:nomanagedalloc)
CXX="nvc++" cmake -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" \
  -D CMAKE_MODULE_PATH=cmake/fetch_content \
  -D RDBENCH_USE_GPU=ON \
  -D RDBENCH_USE_UNIFIED_MEMORY_NO_MANAGED_ALLOC=ON
cmake --build build/

# (optional) build with multicore support (nvc++ -stdpar=multicore)
# Hybrid MPI+OpenMP build
# nvc++ -stdpar=multicore uses OpenMP runtime via Thrust’s OpenMP back-end.
CXX="nvc++" cmake -S . -B build \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" \
  -D CMAKE_MODULE_PATH=cmake/fetch_content \
  -D RDBENCH_USE_MULTICORE=ON
cmake --build build/
```

# Usage
see: https://github.com/tsukuba-hpcs/rdbench-eval
```bash
# file output (Flat MPI run)
mkdir out/ # output directory
nnodes=4
ppn=4 # processes per node
np=$((nnodes * ppn))
tpp=1 # threads per process
mpirun \
  --hostfile "${PBS_NODEFILE}" \
  --np $np \
  --map-by "ppr:${ppn}:node:pe=${tpp}" \
build/rdbench \
  --sz_tile_x 128 \ # each process has 128x128 tiles
  --sz_tile_y 128 \
  --nr_tiles_x 0 \ # 0: auto (default) 
  --nr_tiles_y 0 \ # np=16 may become 4x4 cartesian communicator
  --steps 100 \
  --interval 10 \ # output every 10 steps
  --init_output \ # output initial state (step 0)
  --output out/flat_mpi_ \ # output prefix
  --prettify

# no file output, compute only
# Flat MPI run
nnodes=4
ppn=4 # processes per node
np=$((nnodes * ppn))
tpp=1 # threads per process
mpirun \
  --hostfile "${PBS_NODEFILE}" \
  --np $np \
  --map-by "ppr:${ppn}:node:pe=${tpp}" \
build/rdbench \
  -s 100 \
  -i 0 # no output

# MPI + GPU run
nnodes=4
ppn=1 # processes per node
np=$((nnodes * ppn))
tpp=1 # threads per process
mpirun \
  -x OMP_NUM_THREADS=$tpp \ # if building with RDBENCH_USE_MULTICORE=ON
  --hostfile "${PBS_NODEFILE}" \
  --np $np \
  --map-by "ppr:${ppn}:node:pe=${tpp}" \
build/rdbench \
  -s 100 \
  -i 0

# Hybrid MPI run
nnodes=4
ppn=1 # processes per node
np=$((nnodes * ppn))
tpp=4 # threads per process
mpirun \
  -x OMP_NUM_THREADS=$tpp \ # if building with RDBENCH_USE_MULTICORE=ON
  --hostfile "${PBS_NODEFILE}" \
  --np $np \
  --map-by "ppr:${ppn}:node:pe=${tpp}" \
build/rdbench \
  -s 100 \
  -i 0
```

# Other demos
![](https://raw.githubusercontent.com/range3/rdbench/master/rdbench-viz/viz2.gif)
