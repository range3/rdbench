#pragma once

#include <prof/event.hpp>
#include <prof/recorder.hpp>

namespace rdbench::v2 {

class profiler {
 public:
  using io_phase = prof::event<"IO">;
  using compute_phase = prof::event<"Comp">;

  profiler() = default;

  static void init() {
    auto& recorder = prof::recorder::instance();
    recorder.register_event(io_phase::instance());
    recorder.register_event(compute_phase::instance());
  }
};

}  // namespace rdbench::v2
