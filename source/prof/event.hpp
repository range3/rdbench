#pragma once
#include <chrono>
#include <string_view>

#include "stopwatch.hpp"

namespace prof {

template <std::size_t N>
struct fixed_string {
  // NOLINTNEXTLINE
  constexpr fixed_string(const char (&str)[N]) {
    std::copy_n(str, N, value.begin());
  }

  std::array<char, N> value{};

  friend constexpr auto operator<=>(const fixed_string&,
                                    const fixed_string&) = default;
  friend constexpr auto operator==(const fixed_string&,
                                   const fixed_string&) -> bool = default;

  [[nodiscard]]
  constexpr auto view() const -> std::string_view {
    return std::string_view(value.data(), value.size() - 1);
  }
};

// deduction guide
template <std::size_t N>
// NOLINTNEXTLINE
fixed_string(const char (&)[N]) -> fixed_string<N>;

class event_base {
 public:
  using stopwatch_type = basic_stopwatch<std::chrono::high_resolution_clock>;
  using duration = stopwatch_type::duration;

  event_base() = default;
  event_base(const event_base&) = delete;
  auto operator=(const event_base&) -> event_base& = delete;
  event_base(event_base&&) = delete;
  auto operator=(event_base&&) -> event_base& = delete;
  virtual ~event_base() = default;

  [[nodiscard]]
  virtual auto name() const -> std::string_view = 0;

  void add(duration d) noexcept {
    total_time_ += d;
    max_time_ = std::max(max_time_, d);
    ++count_;
  }

  void clear() noexcept {
    total_time_ = duration::zero();
    max_time_ = duration::zero();
    count_ = 0;
  }

  [[nodiscard]]
  constexpr auto total_time() const noexcept -> duration {
    return total_time_;
  }

  [[nodiscard]]
  constexpr auto max_time() const noexcept -> duration {
    return max_time_;
  }

  [[nodiscard]]
  constexpr auto count() const noexcept -> std::size_t {
    return count_;
  }

 private:
  duration total_time_{};
  duration max_time_{};
  std::size_t count_{};
};

template <fixed_string Name>
class event : public event_base {
 protected:
  event() = default;

 public:
  event(const event&) = delete;
  auto operator=(const event&) -> event& = delete;
  event(event&&) = delete;
  auto operator=(event&&) -> event& = delete;
  ~event() override = default;

  static auto instance() -> event& {
    static event instance;
    return instance;
  }

  [[nodiscard]]
  auto name() const -> std::string_view override {
    return Name.view();
  }
};

}  // namespace prof
