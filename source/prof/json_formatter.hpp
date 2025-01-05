#pragma once

#include <nlohmann/json.hpp>

#include "prof/event.hpp"

namespace prof {

struct event_stats {
  std::string_view name;
  event_base::duration total_time;
  event_base::duration max_time;
  std::size_t count;

  explicit event_stats(const event_base& event)
      : name(event.name()),
        total_time(event.total_time()),
        max_time(event.max_time()),
        count(event.count()) {}
};

struct profile_results {
  event_base::duration total_time;
  std::vector<event_stats> events;
};

}  // namespace prof

namespace nlohmann {

template <>
struct adl_serializer<prof::profile_results> {
  static void to_json(ordered_json& j, const prof::profile_results& results) {
    using duration_sec = std::chrono::duration<double>;
    const auto total_time_sec =
        std::chrono::duration_cast<duration_sec>(results.total_time).count();

    j = ordered_json{{"totalTime", total_time_sec}};

    for (const auto& event : results.events) {
      const auto event_time_sec =
          std::chrono::duration_cast<duration_sec>(event.total_time).count();
      const auto max_time_sec =
          std::chrono::duration_cast<duration_sec>(event.max_time).count();
      const auto avg_time =
          event.count > 0 ? event_time_sec / static_cast<double>(event.count)
                          : 0.0;

      const auto name = std::string{event.name};
      j[name + "TotalTime"] = event_time_sec;
      j[name + "MaxTime"] = max_time_sec;
      j[name + "Count"] = event.count;
      j[name + "AvgTime"] = avg_time;
      j[name + "Percentage"] = (event_time_sec / total_time_sec * 100.0);
    }
  }
};

}  // namespace nlohmann

namespace prof {

class json_formatter {
 public:
  using result_type = nlohmann::ordered_json;
  using duration = event_base::duration;

  [[nodiscard]]
  static auto format(
      const std::vector<std::reference_wrapper<event_base>>& events,
      duration total_time) -> result_type {
    std::vector<event_stats> event_statistics;
    event_statistics.reserve(events.size());

    for (const auto& event_ref : events) {
      event_statistics.emplace_back(event_ref.get());
    }

    return profile_results{total_time, std::move(event_statistics)};
  }
};

}  // namespace prof
