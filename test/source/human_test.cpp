#include <cstdint>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "rdbench/v2/utils/human_readable.hpp"

TEST_CASE("power_basic_operations", "[power]") {
  using namespace rdbench::v2::utils;  // NOLINT

  SECTION("positive exponent") {
    static_assert(power<2, 3>() == 8);
    static_assert(power<3, 4>() == 81);
    static_assert(power<5, 2>() == 25);

    static_assert(power<2, 0>() == 1);
    static_assert(power<42, 0>() == 1);

    static_assert(power<2, 1>() == 2);
    static_assert(power<42, 1>() == 42);

    static_assert(power<2, 10>() == 1024);
    static_assert(power<10, 3>() == 1000);
  }
}

TEST_CASE("format_number basic operations", "[format_number]") {
  using namespace rdbench::v2::utils;  // NOLINT
  CHECK(detail::format_number(123) == "123");
  CHECK(detail::format_number<4>(123) == "123");
  CHECK(detail::format_number(12345) == "12345");
  CHECK(detail::format_number<2>(12345) == "12345");
  CHECK(detail::format_number(123.456F) == "123");
  CHECK(detail::format_number<4>(0.123456789) == "0.123");
  CHECK(detail::format_number(1234.5678F) == "1235");
  CHECK(detail::format_number<6>(1234.5678) == "1234.57");
  CHECK(detail::format_number(-123) == "-123");
  CHECK(detail::format_number<4>(-123) == "-123");
  CHECK(detail::format_number(-12345) == "-12345");
  CHECK(detail::format_number<2>(-12345) == "-12345");
  CHECK(detail::format_number(-123.456F) == "-123");
  CHECK(detail::format_number<4>(-0.123456789) == "-0.123");
  CHECK(detail::format_number(-1234.5678F) == "-1235");
  CHECK(detail::format_number<6>(-1234.5678) == "-1234.57");
  CHECK(detail::format_number(2.) == "2.00");
  CHECK(detail::format_number(static_cast<uint64_t>(10000000000000000000ULL))
        == "10000000000000000000");
}

TEST_CASE("to_human basic operations", "[to_human]") {
  using namespace rdbench::v2::utils;  // NOLINT

  CHECK(to_human(1024) == "1K");
  CHECK(to_human<1000>(1024) == "1.02K");
  CHECK(to_human<1000>(0.00123) == "1.23m");
  CHECK(to_human(static_cast<uint64_t>(1024 * 1024 * 1024)) == "1G");
  CHECK(to_human(123.456) == "123");
  CHECK(to_human(123.456F) == "123");
  CHECK(to_human(0) == "0");
  CHECK(to_human(500) == "500");
  CHECK(to_human(0.5) == "512m");
  CHECK(to_human<1000>(1000000) == "1M");
}

TEST_CASE("from_human basic operations", "[from_human]") {
  using namespace rdbench::v2::utils;  // NOLINT

  CHECK(from_human<int>("1K") == 1024);
  CHECK(from_human<int>("1k") == 1024);
  CHECK(from_human<int, 1000>("1K") == 1000);
  CHECK(from_human<int, 1000>("1k") == 1000);
  CHECK(from_human<float, 1000>("1.23m") == Catch::Approx(0.00123));
  CHECK(from_human<double>("1G") == Catch::Approx(1U << 30U));
  CHECK(from_human<short>("1K") == 1024);
  CHECK(from_human<unsigned short>("32K") == (1U << 15U));
}
