#pragma once

#include <functional>
#include <vector>

#include "event.hpp"
#include "formatter.hpp"

namespace prof {

class running_recorder;

class recorder {
  friend class prof::running_recorder;

 public:
  using stopwatch_type = event_base::stopwatch_type;
  using duration = event_base::duration;

 protected:
  recorder() = default;

 public:
  static auto instance() -> recorder& {
    static recorder instance;
    return instance;
  }

  void register_event(event_base& e) { events_.emplace_back(e); }

  [[nodiscard]]
  auto is_running() const noexcept -> bool {
    return running_;
  }

  [[nodiscard]] auto total_time() const -> duration { return total_time_; }

  auto start() -> running_recorder;

  void clear() {
    for (auto& e : events_) {
      e.get().clear();
    }
    total_time_ = duration::zero();
  }

  template <formatter F>
  auto format(F&& f) const -> typename F::result_type {
    return std::forward<F>(f).format(events_, total_time_);
  }

  template <formatter F>
  auto format() const -> typename F::result_type {
    return format(F{});
  }

 private:
  void stop() noexcept {
    auto d = sw_.get();
    if (current_phase_ != nullptr) {
      current_phase_->add(d);
      current_phase_ = nullptr;
    }
    total_time_ += d;
    running_ = false;
  }

  template <typename Phase>
  void switch_phase() {
    auto& next_phase = Phase::instance();
    if (current_phase_ != nullptr) {
      auto d = sw_.get_and_reset();
      current_phase_->add(d);
      total_time_ += d;
    } else {
      sw_.reset();
    }
    current_phase_ = &next_phase;
  }

  std::vector<std::reference_wrapper<event_base>> events_;
  stopwatch_type sw_;
  duration total_time_{};
  event_base* current_phase_{};
  bool running_{false};
};

class running_recorder {
 public:
  explicit running_recorder(recorder& r) : recorder_{r} {}
  ~running_recorder() noexcept { recorder_.get().stop(); }
  running_recorder(const running_recorder&) = delete;
  auto operator=(const running_recorder&) -> running_recorder& = delete;
  running_recorder(running_recorder&&) = delete;
  auto operator=(running_recorder&&) -> running_recorder& = delete;

  template <typename Phase>
  void switch_phase() {
    recorder_.get().switch_phase<Phase>();
  }

 private:
  std::reference_wrapper<recorder> recorder_;
};

inline auto recorder::start() -> running_recorder {
  if (is_running()) {
    throw std::logic_error{"recorder is already running"};
  }
  running_ = true;
  sw_.reset();
  return running_recorder{*this};
}

}  // namespace prof
