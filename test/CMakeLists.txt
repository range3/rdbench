# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(rdbenchTests LANGUAGES CXX)

# ---- Dependencies ----

find_package(Catch2 REQUIRED)
include(Catch)

# ---- Tests ----

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/source/*_test.cpp")
add_executable(rdbench_test ${TEST_SOURCES})
target_link_libraries(
    rdbench_test PRIVATE
    rdbench_lib
    Catch2::Catch2WithMain
)
target_compile_features(rdbench_test PRIVATE cxx_std_20)

catch_discover_tests(rdbench_test)

# ---- End-of-file commands ----

add_folders(Test)
