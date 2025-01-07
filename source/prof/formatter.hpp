#pragma once

#include <functional>
#include <iomanip>
#include <sstream>

#include <prof/event.hpp>

namespace prof {

template <typename F>
concept has_static_format =
    requires(const std::vector<std::reference_wrapper<event_base>>& events,
             event_base::duration total_time) {
      {
        F::format(events, total_time)
      } -> std::same_as<typename F::result_type>;
    };

template <typename F>
concept has_instance_format =
    requires(F f,
             const std::vector<std::reference_wrapper<event_base>>& events,
             event_base::duration total_time) {
      { f.format(events, total_time) } -> std::same_as<typename F::result_type>;
    };

template <typename F>
concept formatter = requires {
  typename F::result_type;
  requires std::movable<typename F::result_type>;
} && (has_static_format<F> || has_instance_format<F>);

class text_formatter {
 public:
  using result_type = std::string;
  using duration = event_base::duration;

  [[nodiscard]]
  static auto format(
      const std::vector<std::reference_wrapper<event_base>>& events,
      duration total_time) -> result_type {
    using duration_sec = std::chrono::duration<double>;
    std::stringstream ss;
    for (const auto& e : events) {
      const auto& event = e.get();
      auto event_time_sec =
          std::chrono::duration_cast<duration_sec>(event.total_time()).count();
      auto total_time_sec =
          std::chrono::duration_cast<duration_sec>(total_time).count();
      auto count = event.count();
      auto event_max_sec =
          std::chrono::duration_cast<duration_sec>(event.max_time()).count();

      ss << std::fixed << std::setprecision(6) << "  " << std::setw(22)
         << std::left << event.name() << ": " << std::setw(10) << std::right
         << event_time_sec / total_time_sec * 100.0 << " % " << "( "
         << std::setw(15) << event_time_sec << " s / " << std::setw(15)
         << total_time_sec << " s ) " << "count: " << std::setw(10) << count
         << " ave: " << std::setw(8)
         << (count > 0 ? event_time_sec / static_cast<double>(count) : 0)
         << " s max: " << std::setw(8) << event_max_sec << " s\n";
    }
    return ss.str();
  }
};

}  // namespace prof
