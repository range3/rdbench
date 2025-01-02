#pragma once

#include <cmath>
#include <concepts>
#include <limits>
#include <ratio>

namespace rdbench::v2::utils {

namespace detail {

template <uintmax_t Base, unsigned int Exp>
consteval auto power_positive() -> uintmax_t {
  if constexpr (Exp == 0) {
    return 1;
  }
  if constexpr (Exp == 1) {
    return Base;
  }

  uintmax_t const half = power_positive<Base, Exp / 2>();
  if constexpr (Exp % 2 == 0) {
    return half * half;
  } else {
    return Base * half * half;
  }
}

}  // namespace detail

template <uintmax_t Base, int Exp>
consteval auto power() {
  if constexpr (Exp >= 0) {
    return std::ratio<
        detail::power_positive<Base, static_cast<unsigned int>(Exp)>(), 1>{};
  } else {
    return std::ratio<
        1, detail::power_positive<Base, static_cast<unsigned int>(-Exp)>()>{};
  }
}

template <uintmax_t Base, int Exp, typename FloatT = double>
consteval auto power_as_floating() -> FloatT {
  constexpr auto ratio = power<Base, Exp>();
  return static_cast<FloatT>(ratio.num) / static_cast<FloatT>(ratio.den);
}

template <uintmax_t Base, int Exp>
consteval auto power_as_integral() {
  static_assert(Exp >= 0, "Negative exponent not allowed for integral result");
  constexpr auto ratio = power<Base, Exp>();
  static_assert(ratio.den == 1, "Result must be an integer");
  return ratio.num;
}

}  // namespace rdbench::v2::utils
