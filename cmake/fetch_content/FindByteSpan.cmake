cmake_policy(SET CMP0169 OLD)
include(FetchContent)
FetchContent_Declare(
    ByteSpan 
    GIT_REPOSITORY https://github.com/range3/byte-span.git
    GIT_TAG v0.2.0
    GIT_SHALLOW YES
)

FetchContent_GetProperties(ByteSpan)
if(NOT ByteSpan_POPULATED)
  FetchContent_Populate(ByteSpan)
  add_subdirectory("${bytespan_SOURCE_DIR}" "${bytespan_BINARY_DIR}" EXCLUDE_FROM_ALL)
endif()

set(ByteSpan_FOUND 1)
