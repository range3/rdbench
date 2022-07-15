#pragma once

#include <chrono>
#include <vector>

class Stopwatch {
public:
  using clock = std::chrono::system_clock;
  using duration = clock::duration;
  using timepoint = clock::time_point;

  void reset() { start_ = clock::now(); }

  duration get_and_reset() {
    auto now = clock::now();
    auto elapsed = now - start_;
    start_ = std::move(now);
    return elapsed;
  }

  duration get() const { return clock::now() - start_; }

private:
  timepoint start_;
};
