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

# ---- Declare library ----

add_library(
    rdbench_lib OBJECT
    source/lib.cpp
)

target_include_directories(
    rdbench_lib ${warning_guard}
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(rdbench_lib PUBLIC cxx_std_20)

find_package(fmt REQUIRED)
find_package(ByteSpan REQUIRED)
find_package(cxxmpi REQUIRED)
find_package(cxxopts REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(MPI REQUIRED)
find_package(kamping REQUIRED)
target_link_libraries(rdbench_lib
  PUBLIC
    MPI::MPI_CXX
    kamping::kamping
    ByteSpan::ByteSpan
    cxxmpi::cxxmpi
  PRIVATE
    fmt::fmt
)

# ---- Declare executable ----

add_executable(rdbench_exe source/main.cpp)
add_executable(rdbench::exe ALIAS rdbench_exe)

set_property(TARGET rdbench_exe PROPERTY OUTPUT_NAME rdbench)

target_compile_features(rdbench_exe PRIVATE cxx_std_20)

target_link_libraries(rdbench_exe PRIVATE rdbench_lib cxxopts::cxxopts nlohmann_json::nlohmann_json)

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
