#include <chrono>
#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>

#include <catch2/catch_test_macros.hpp>
#include <cxxabi.h>

#include "prof/event.hpp"

template <typename T>
auto type_name() -> std::string {
  int status = 0;
  const std::unique_ptr<char, void (*)(void*)> result(
      abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status),
      std::free);
  return status == 0 ? result.get() : typeid(T).name();
}

TEST_CASE("event basic operation", "[event]") {
  SECTION("event name") {
    auto& e1 = prof::event<"IO">::instance();
    auto& e2 = prof::event<"OI">::instance();
    CHECK(e1.name() == "IO");
    CHECK(e2.name() == "OI");
    CHECK(e1.name() != e2.name());

    // std::cout << "e1: " << type_name<decltype(e1)>() << '\n';
    // std::cout << "e2: " << type_name<decltype(e2)>() << '\n';

    CHECK(type_name<decltype(e1)>() != type_name<decltype(e2)>());
  }

  SECTION("type matching") {
    static_assert(!std::is_same_v<prof::event<"IO">, prof::event<"OI">>);
    static_assert(std::is_same_v<prof::event<"IO">, prof::event<"IO">>);
  }
}

TEST_CASE("Event basic functionality", "[event]") {
  using std::chrono_literals::operator""s;  // NOLINT

  // Define an event with compile-time name
  auto& event = prof::event<"test_event">::instance();
  using duration = prof::event<"test_event">::duration;

  SECTION("Check initial state") {
    CHECK(event.name() == "test_event");
    CHECK(event.count() == 0);
    CHECK(event.total_time() == duration::zero());
    CHECK(event.max_time() == duration::zero());
  }

  SECTION("Add duration and check statistics") {
    // Record 1 second processing time
    event.add(std::chrono::duration_cast<duration>(1.0s));

    CHECK(event.count() == 1);
    CHECK(event.total_time() == 1.0s);
    CHECK(event.max_time() == 1.0s);

    // Add another 0.5 second processing time
    event.add(std::chrono::duration_cast<duration>(0.5s));

    CHECK(event.count() == 2);
    CHECK(event.total_time() == 1.5s);
    CHECK(event.max_time() == 1.0s);  // max remains 1 second
  }

  SECTION("Check clear functionality") {
    event.add(std::chrono::duration_cast<duration>(1.0s));
    event.add(std::chrono::duration_cast<duration>(0.5s));

    event.clear();

    CHECK(event.count() == 0);
    CHECK(event.total_time() == duration::zero());
    CHECK(event.max_time() == duration::zero());
  }

  SECTION("Check max time updates") {
    event.add(std::chrono::duration_cast<duration>(1.0s));
    event.add(std::chrono::duration_cast<duration>(2.0s));  // Add larger value
    event.add(std::chrono::duration_cast<duration>(0.5s));

    CHECK(event.max_time() == 2.0s);  // max should be 2 seconds
    CHECK(event.total_time() == 3.5s);
    CHECK(event.count() == 3);
  }
}
