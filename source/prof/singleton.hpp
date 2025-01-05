#pragma once

#include <cassert>
#include <optional>
#include <utility>

namespace prof {

template <typename T>
class singleton {
  static_assert(!std::is_array_v<T>, "Array types are not supported");
  static_assert(std::is_destructible_v<T>, "T must be destructible");

 public:
  using instance_type = T;

  static auto instance() noexcept -> instance_type& {
    assert(initialized());
    return *instance_storage();
  }

  [[nodiscard]]
  static constexpr auto initialized() -> bool {
    return instance_storage().has_value();
  }

  template <typename... Args>
  static void init(Args&&... args) {
    assert(!initialized());
    instance_storage().emplace(std::forward<Args>(args)...);
  }

  static void fini() {
    assert(initialized());
    instance_storage().reset();
  }

 private:
  static auto instance_storage() -> std::optional<T>& {
    static std::optional<T> instance_storage;
    return instance_storage;
  }
};

}  // namespace prof
