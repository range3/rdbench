
cmake_policy(SET CMP0169 OLD)
include(FetchContent)
FetchContent_Declare(
  kamping
  GIT_REPOSITORY https://github.com/kamping-site/kamping.git
  GIT_TAG v0.1.1
  GIT_SHALLOW YES
)
FetchContent_MakeAvailable(kamping)

# FetchContent_GetProperties(kamping)
# if(NOT kamping_POPULATED)
#   FetchContent_Populate(kamping)
#   add_subdirectory("${kamping_SOURCE_DIR}" "${kamping_BINARY_DIR}" EXCLUDE_FROM_ALL)
# endif()

set(kamping_FOUND 1)
