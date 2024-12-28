cmake_policy(SET CMP0169 OLD)
include(FetchContent)
FetchContent_Declare(
    cxxmpi 
    GIT_REPOSITORY https://github.com/range3/cxxmpi.git
    GIT_TAG master
    GIT_REMOTE_UPDATE_STRATEGY REBASE
)

FetchContent_GetProperties(cxxmpi)
if(NOT cxxmpi_POPULATED)
  FetchContent_Populate(cxxmpi)
  add_subdirectory("${cxxmpi_SOURCE_DIR}" "${cxxmpi_BINARY_DIR}" EXCLUDE_FROM_ALL)
endif()

set(cxxmpi_FOUND 1)
