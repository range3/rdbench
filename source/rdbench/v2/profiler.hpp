#pragma once

#include <prof/event.hpp>
#include <prof/recorder.hpp>

namespace rdbench::v2 {

class profiler {
 public:
  using io_phase = prof::event<"IO">;
  using compute_phase = prof::event<"Comp">;
  using compute_next_state = prof::event<"CompNextState">;
  using exchange_halos = prof::event<"ExchangeHalos">;

  profiler() = default;

  static void init() {
    auto& recorder = prof::recorder::instance();
    recorder.register_event(io_phase::instance());
    recorder.register_event(compute_phase::instance());
    recorder.register_event(compute_next_state::instance());
    recorder.register_event(exchange_halos::instance());
  }
};

}  // namespace rdbench::v2
