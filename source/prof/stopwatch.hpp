#pragma once
#include <chrono>
#include <type_traits>

namespace prof {

template <typename Clock = std::chrono::high_resolution_clock,
          typename Rep = typename Clock::rep,
          typename Period = typename Clock::period>
class basic_stopwatch {
 public:
  using clock = Clock;
  using rep = Rep;
  using period = Period;
  using duration = std::chrono::duration<rep, period>;
  using time_point = std::chrono::time_point<clock, duration>;

  basic_stopwatch() = default;

  void reset() { reset(clock::now()); }

  [[nodiscard]]
  auto get() const -> duration {
    if constexpr (std::is_same_v<duration, typename clock::duration>) {
      return clock::now() - start_time_;
    } else {
      return std::chrono::duration_cast<duration>(clock::now() - start_time_);
    }
  }

  [[nodiscard]]
  auto get_and_reset() -> duration {
    auto current_time = clock::now();
    auto elapsed_time = current_time - start_time_;
    reset(current_time);
    if constexpr (std::is_same_v<duration, typename clock::duration>) {
      return elapsed_time;
    } else {
      return std::chrono::duration_cast<duration>(elapsed_time);
    }
  }

 private:
  void reset(time_point new_start_time) noexcept {
    start_time_ = new_start_time;
  }

  typename clock::time_point start_time_{clock::now()};
};

using stopwatch =
    basic_stopwatch<std::chrono::high_resolution_clock, double, std::ratio<1>>;

}  // namespace prof
