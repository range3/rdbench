cmake_minimum_required(VERSION 3.14...3.22)

# ---- Project ----

# Note: update this to your new project's name and version
project(
  Rdbench
  VERSION 0.1.0
  LANGUAGES CXX
)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---- Set default build type ----
# Encourage user to specify a build type (e.g. Release, Debug, etc.), otherwise
# set it to Release.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# ---- Add dependencies via CPM ----
# see https://github.com/TheLartians/CPM.cmake for more info

include(cmake/CPM.cmake)

# PackageProject.cmake will be used to make our target installable
CPMAddPackage("gh:TheLartians/PackageProject.cmake@1.8.0")

CPMAddPackage(
  NAME fmt
  GIT_TAG 8.1.1
  GITHUB_REPOSITORY fmtlib/fmt
  OPTIONS "FMT_INSTALL YES" # create an installable target
  EXCLUDE_FROM_ALL YES
)

CPMAddPackage(
  GITHUB_REPOSITORY jarro2783/cxxopts
  VERSION 3.0.0
  OPTIONS "CXXOPTS_BUILD_EXAMPLES NO" "CXXOPTS_BUILD_TESTS NO" "CXXOPTS_ENABLE_INSTALL YES"
  EXCLUDE_FROM_ALL YES
)

CPMAddPackage(
  NAME nlohmann_json
  VERSION 3.9.1
  GITHUB_REPOSITORY nlohmann/json
  OPTIONS 
    "JSON_BuildTests OFF"
  EXCLUDE_FROM_ALL YES
)

# ---- Other dependencies ----
find_package(MPI REQUIRED)

# ---- Add source files ----

# Note: globbing sources is considered bad practice as CMake's generators may not detect new files
# automatically. Keep that in mind when changing files, or explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/rdbench/*.hpp")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/rdbench/*.cpp")

# ---- Create executable ----
# rdbench
add_executable(${PROJECT_NAME} ${headers} ${sources})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "rdbench")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt cxxopts::cxxopts nlohmann_json::nlohmann_json MPI::MPI_CXX)
install(TARGETS ${PROJECT_NAME} DESTINATION "${CMAKE_INSTALL_BINDIR}")

# rdbench_seq
add_executable(${PROJECT_NAME}_seq "${CMAKE_CURRENT_SOURCE_DIR}/src/rdbench_seq/main.cpp")
set_target_properties(${PROJECT_NAME}_seq PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "rdbench_seq")
target_compile_options(${PROJECT_NAME}_seq PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")
target_link_libraries(${PROJECT_NAME}_seq PRIVATE fmt::fmt cxxopts::cxxopts nlohmann_json::nlohmann_json MPI::MPI_CXX)

# rdbench_int
add_executable(${PROJECT_NAME}_int "${CMAKE_CURRENT_SOURCE_DIR}/src/rdbench_int/main.cpp")
set_target_properties(${PROJECT_NAME}_int PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "rdbench_int")
target_compile_options(${PROJECT_NAME}_int PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")
target_link_libraries(${PROJECT_NAME}_int PRIVATE fmt::fmt cxxopts::cxxopts nlohmann_json::nlohmann_json MPI::MPI_CXX)
