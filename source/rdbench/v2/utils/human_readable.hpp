#pragma once

#include <array>
#include <charconv>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>

namespace rdbench::v2::utils {

template <uintmax_t Base, unsigned int Exp>
consteval auto power() -> uintmax_t {
  if constexpr (Exp == 0) {
    return 1;
  }
  if constexpr (Exp == 1) {
    return Base;
  }

  uintmax_t const half = power<Base, Exp / 2>();
  if constexpr (Exp % 2 == 0) {
    return half * half;
  } else {
    return Base * half * half;
  }
}

namespace detail {

template <typename T>
constexpr auto is_approximately_zero(T value) -> bool {
  return std::abs(value) < std::numeric_limits<T>::epsilon() * 100;
}

template <int DecimalFlexLimit = 3, typename T>
  requires std::floating_point<T> || std::integral<T>
auto format_number(T value) -> std::string {
  if constexpr (std::is_integral_v<T>) {
    return std::to_string(value);
  } else {
    std::ostringstream os;
    int integral_digits =
        (is_approximately_zero(value))
            ? 1
            : static_cast<int>(std::log10(std::abs(static_cast<double>(value))))
                  + 1;

    int precision = std::max(0, DecimalFlexLimit - integral_digits);
    double rounding_factor = std::pow(10, precision);
    const auto rounded_value =
        std::round(static_cast<double>(value) * rounding_factor)
        / rounding_factor;

    os << std::fixed << std::setprecision(precision) << rounded_value;
    auto str = os.str();

    return str;
  }
}

template <std::floating_point T>
constexpr auto calculate_exponent(T value, unsigned int base) -> int {
  if (is_approximately_zero(value)) {
    return 0;
  }
  return static_cast<int>(
      std::floor(std::log(std::abs(static_cast<double>(value)))
                 / std::log(static_cast<double>(base))));
}

}  // namespace detail

template <unsigned int Base = 1024, typename T = double>
auto to_human(T value) -> std::string {
  using namespace std::string_view_literals;
  static auto units =
      std::array{"y"sv, "z"sv, "a"sv, "f"sv, "p"sv, "n"sv, "u"sv, "m"sv,
                 ""sv,  "K"sv, "M"sv, "G"sv, "T"sv, "P"sv, "E"sv};

  auto dbl_value = static_cast<double>(value);
  const int exponent = detail::calculate_exponent(dbl_value, Base);
  if (exponent < -8 || exponent > 6) {
    return detail::format_number<3>(value);
  }

  const double scaled_value = dbl_value / std::pow(Base, exponent);
  auto str = detail::format_number<3>(scaled_value);
  if (str.ends_with(".00")) {
    str.resize(str.size() - 3);
  }

  return str + std::string{units[static_cast<size_t>(exponent + 8)]};
}

template <typename T, unsigned int Base = 1024>
auto from_human(std::string_view sv) -> T {
  static const auto unit_map = std::unordered_map<char, std::pair<T, T>>{
      {'y', {1, power<Base, 8>()}}, {'z', {1, power<Base, 7>()}},
      {'a', {1, power<Base, 6>()}}, {'f', {1, power<Base, 5>()}},
      {'p', {1, power<Base, 4>()}}, {'n', {1, power<Base, 3>()}},
      {'u', {1, power<Base, 2>()}}, {'m', {1, power<Base, 1>()}},
      {'K', {power<Base, 1>(), 1}}, {'k', {power<Base, 1>(), 1}},
      {'M', {power<Base, 2>(), 1}}, {'G', {power<Base, 3>(), 1}},
      {'T', {power<Base, 4>(), 1}}, {'P', {power<Base, 5>(), 1}},
      {'E', {power<Base, 6>(), 1}},
  };

  T value;

  auto [ptr, ec] = std::from_chars(sv.data(), sv.end(), value);
  if (ec == std::errc::invalid_argument) {
    throw std::invalid_argument("Invalid number format");
  }

  if (ptr == sv.end()) {
    return value;
  }

  if (auto unit_iter = unit_map.find(*ptr); unit_iter != unit_map.end()) {
    return static_cast<T>(value * unit_iter->second.first
                          / unit_iter->second.second);
  }

  throw std::invalid_argument("Unknown unit suffix");
}

}  // namespace rdbench::v2::utils
