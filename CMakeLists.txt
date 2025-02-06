cmake_minimum_required(VERSION 3.14...3.31)

include(cmake/prelude.cmake)

project(
    rdbench
    VERSION 2.0.0
    DESCRIPTION "2D reaction-diffusion system benchmark using MPI and MPI-IO."
    HOMEPAGE_URL "https://github.com/range3/rdbench"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)
include(cmake/stdpar.cmake)

# ---- Options ----
option(RDBENCH_USE_GPU "Enable GPU offloading with NVC++" OFF)
option(RDBENCH_USE_UNIFIED_MEMORY "Enable unified memory for GPU execution" OFF)
option(RDBENCH_USE_UNIFIED_MEMORY_NO_MANAGED_ALLOC "Enable to force non-managed allocation for unified memory" OFF)
option(RDBENCH_USE_MULTICORE "Enable multicore execution" OFF)
option(RDBENCH_USE_NVTX "Enable NVIDIA Tools Extension (NVTX) profiling" OFF)
rdbench_setup_stdpar_options()

# ---- Declare library ----

add_library(rdbench_lib INTERFACE)

target_include_directories(
    rdbench_lib ${warning_guard}
    INTERFACE
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_definitions(rdbench_lib INTERFACE MDSPAN_USE_PAREN_OPERATOR=1)
if(RDBENCH_USE_STDPAR)
  target_compile_features(rdbench_lib INTERFACE cxx_std_23)
else()
  target_compile_features(rdbench_lib INTERFACE cxx_std_20)
endif()
rdbench_configure_parallel_execution(rdbench_lib)

find_package(fmt REQUIRED)
find_package(cxxmpi REQUIRED)
find_package(cxxopts REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(mdspan REQUIRED)
find_package(MPI REQUIRED)
target_link_libraries(rdbench_lib
  INTERFACE
    MPI::MPI_CXX
    fmt::fmt
    cxxmpi::cxxmpi
    std::mdspan
)

# NVTX support
if(RDBENCH_USE_NVTX)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(rdbench_lib
        INTERFACE
        CUDA::nvtx3
    )
    target_compile_definitions(rdbench_lib
        INTERFACE
        RDBENCH_USE_NVTX
    )
else()
    target_compile_definitions(rdbench_lib
        INTERFACE
        NVTX_DISABLE
    )
endif()

# ---- Declare executable ----

add_executable(rdbench_exe source/main.cpp)
add_executable(rdbench::exe ALIAS rdbench_exe)

set_property(TARGET rdbench_exe PROPERTY OUTPUT_NAME rdbench)

target_link_libraries(rdbench_exe
  PRIVATE
    rdbench_lib
    cxxopts::cxxopts
    nlohmann_json::nlohmann_json
)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT rdbench_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of rdbench"
  )
endif()

include(cmake/dev-mode.cmake)
