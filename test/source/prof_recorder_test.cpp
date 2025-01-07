#include <chrono> // NOLINT
#include <functional>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "prof/event.hpp"
#include "prof/formatter.hpp"
#include "prof/recorder.hpp"

using namespace std::chrono_literals;

namespace {

class test_formatter {
 public:
  using result_type = std::string;

  static auto format(
      const std::vector<std::reference_wrapper<prof::event_base>>& /*events*/,
      prof::event_base::duration /*unused*/) -> result_type {
    return "test_format";
  }
};

class event_registry {
 public:
  using test_event_1 = prof::event<"test_event_1">;
  using test_event_2 = prof::event<"test_event_2">;

  event_registry() {
    auto& recorder = prof::recorder::instance();
    recorder.register_event(test_event_1::instance());
    recorder.register_event(test_event_2::instance());
  }
};

const event_registry registry{};  // NOLINT

}  // namespace

TEST_CASE("recorder basic functionality", "[recorder]") {
  auto& recorder = prof::recorder::instance();
  recorder.clear();

  SECTION("initial state") {
    CHECK_FALSE(recorder.is_running());
    CHECK(recorder.total_time() == prof::recorder::duration::zero());
  }

  SECTION("register events") {
    auto result = recorder.format<test_formatter>();
    CHECK(result == "test_format");
  }

  SECTION("start and stop recording") {
    {
      auto running = recorder.start();
      CHECK(recorder.is_running());
      std::this_thread::sleep_for(100ms); // NOLINT
    }  // running goes out of scope, stops recording

    CHECK_FALSE(recorder.is_running());
    CHECK(recorder.total_time() > prof::event_base::duration::zero());
  }

  SECTION("phase switching") {
    {
      auto running = recorder.start();
      std::this_thread::sleep_for(50ms); // NOLINT

      running.switch_phase<event_registry::test_event_1>();
      std::this_thread::sleep_for(50ms); // NOLINT

      running.switch_phase<event_registry::test_event_2>();
      std::this_thread::sleep_for(50ms); // NOLINT
    }

    CHECK(event_registry::test_event_1::instance().count() == 1);
    CHECK(event_registry::test_event_2::instance().count() == 1);
  }

  SECTION("clear functionality") {
    {
      auto running = recorder.start();
      std::this_thread::sleep_for(100ms); // NOLINT
    }

    CHECK(recorder.total_time() > prof::event_base::duration::zero());

    recorder.clear();
    CHECK(recorder.total_time() == prof::event_base::duration::zero());
  }
}

TEST_CASE("recorder with text formatter", "[recorder][formatter]") {
  auto& recorder = prof::recorder::instance();
  recorder.clear();

  {
    auto running = recorder.start();
    running.switch_phase<event_registry::test_event_1>();
    std::this_thread::sleep_for(100ms); // NOLINT

    running.switch_phase<event_registry::test_event_2>();
    std::this_thread::sleep_for(50ms); // NOLINT
  }
  {
    const prof::scoped_event<event_registry::test_event_1> event_scope;
    std::this_thread::sleep_for(100ms); // NOLINT
  }

  auto formatted = recorder.format<prof::text_formatter>();
  std::cout << formatted << '\n';
}
