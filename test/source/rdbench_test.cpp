#include <ratio>

#include <catch2/catch_test_macros.hpp>

#include "rdbench/v2/utils/power.hpp"

TEST_CASE("power_basic_operations", "[power]") {
  using namespace rdbench::v2::utils;  // NOLINT

  SECTION("positive exponent") {
    static_assert(power_as_integral<2, 3>() == 8);
    static_assert(power_as_integral<3, 4>() == 81);
    static_assert(power_as_integral<5, 2>() == 25);

    static_assert(power_as_integral<2, 0>() == 1);
    static_assert(power_as_integral<42, 0>() == 1);

    static_assert(power_as_integral<2, 1>() == 2);
    static_assert(power_as_integral<42, 1>() == 42);

    static_assert(power_as_integral<2, 10>() == 1024);
    static_assert(power_as_integral<10, 3>() == 1000);
  }

  SECTION("negative exponent") {
    static_assert(
        std::ratio_equal_v<decltype(power<2, -3>()), std::ratio<1, 8> >);

    static_assert(
        std::ratio_equal_v<decltype(power<3, -2>()), std::ratio<1, 9> >);
  }

  SECTION("edge case") {
    static_assert(power_as_integral<1, 0>() == 1);
    static_assert(power_as_integral<1, 1>() == 1);
    static_assert(power_as_integral<1, 42>() == 1);

    static_assert(power_as_integral<2, 0>() == 1);

    static_assert(
        std::ratio_equal_v<decltype(power<2, 0>()), std::ratio<1, 1> >);
  }
}
