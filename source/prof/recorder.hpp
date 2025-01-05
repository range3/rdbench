#pragma once

#include <atomic>
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
    return running_count_ > 0;
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
  auto begin_recording() -> bool {
    if (running_count_.fetch_add(1, std::memory_order_acq_rel) == 0) {
      sw_.reset();
      return true;
    }
    return false;
  }

  auto end_recording() -> bool {
    if (running_count_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
      stop();
      return true;
    }
    return false;
  }

  void stop() noexcept {
    auto d = sw_.get();
    if (current_phase_ != nullptr) {
      current_phase_->add(d);
      current_phase_ = nullptr;
    }
    total_time_ += d;
  }

  template <typename Phase>
  void switch_phase() {
    auto d = sw_.get_and_reset();
    if (current_phase_ != nullptr) {
      current_phase_->add(d);
    }
    total_time_ += d;
    current_phase_ = &Phase::instance();
  }

  std::vector<std::reference_wrapper<event_base>> events_;
  stopwatch_type sw_;
  duration total_time_{};
  event_base* current_phase_{};
  std::atomic<std::size_t> running_count_{0};
};

class running_recorder {
 public:
  running_recorder() { recorder::instance().begin_recording(); }
  ~running_recorder() noexcept { recorder::instance().end_recording(); }
  running_recorder(const running_recorder&) = delete;
  auto operator=(const running_recorder&) -> running_recorder& = delete;
  running_recorder(running_recorder&&) = delete;
  auto operator=(running_recorder&&) -> running_recorder& = delete;

  template <typename Phase>
  void switch_phase() {
    recorder::instance().switch_phase<Phase>();
  }
};

// NOLINTNEXTLINE
inline auto recorder::start() -> running_recorder {
  return running_recorder{};
}

}  // namespace prof
