#pragma once

#include <cstddef>
#include <memory>
#include <utility>

namespace cxxraii::inline v1 {

template <typename Handle>
struct handle_traits {
  static void free(Handle handle) noexcept {
    // default implementation does nothing
    (void)handle;
  }
};

template <typename NativeHandle, NativeHandle InvalidValue>
class basic_weak_handle {
 public:
  using native_type = NativeHandle;
  static constexpr auto invalid_value = InvalidValue;

  constexpr basic_weak_handle() noexcept = default;
  constexpr basic_weak_handle(std::nullptr_t) noexcept {}  // NOLINT
  constexpr explicit basic_weak_handle(native_type handle) noexcept
      : handle_(handle) {}

  [[nodiscard]]
  explicit operator bool() const noexcept {
    return handle_ != invalid_value;
  }

  // smart reference pattern
  [[nodiscard]]
  constexpr auto operator->() noexcept -> basic_weak_handle* {
    return this;
  }

  [[nodiscard]]
  constexpr auto operator->() const noexcept -> const basic_weak_handle* {
    return this;
  }

  [[nodiscard]]
  constexpr auto native() noexcept -> native_type& {
    return handle_;
  }

  [[nodiscard]]
  constexpr auto native() const noexcept -> native_type {
    return handle_;
  }

  [[nodiscard]]
  constexpr auto get() const noexcept -> basic_weak_handle {
    return *this;
  }

  [[maybe_unused]]
  constexpr auto release() noexcept -> native_type {
    return std::exchange(handle_, invalid_value);
  }

  [[nodiscard]]
  constexpr auto valid() const noexcept -> bool {
    return staic_cast<bool>(*this);
  }

  [[nodiscard]]
  constexpr friend auto operator==(const basic_weak_handle& lhs,
                                   const basic_weak_handle& rhs) noexcept
      -> bool {
    return lhs.handle_ == rhs.handle_;
  }

  [[nodiscard]]
  constexpr friend auto operator!=(const basic_weak_handle& lhs,
                                   const basic_weak_handle& rhs) noexcept
      -> bool {
    return !(lhs == rhs);
  }

 private:
  native_type handle_{invalid_value};
};

template <typename NativeHandle, NativeHandle InvalidValue>
struct basic_deleter {
  using pointer = basic_weak_handle<NativeHandle, InvalidValue>;

  void operator()(pointer handle) const noexcept {
    if (handle) {
      handle_traits<NativeHandle>::free(handle.release());
    }
  }
};

template <typename NativeHandle, NativeHandle InvalidValue>
using basic_managed_handle =
    std::unique_ptr<basic_weak_handle<NativeHandle, InvalidValue>,
                    basic_deleter<NativeHandle, InvalidValue>>;

// helper type bundle
template <typename NativeHandle, NativeHandle InvalidValue>
struct types {
  using native_type = NativeHandle;
  static constexpr auto invalid_value = InvalidValue;
  using weak_handle = basic_weak_handle<NativeHandle, InvalidValue>;
  using managed_handle = basic_managed_handle<NativeHandle, InvalidValue>;
  using deleter_type = basic_deleter<NativeHandle, InvalidValue>;
};

}  // namespace cxxraii::inline v1
