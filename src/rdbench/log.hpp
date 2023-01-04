// SPDX-License-Identifier: Apache-2.0
// Copyright 2022 range3 ( https://github.com/range3 )
#pragma once
#include <fmt/core.h>

#include <chrono>
#include <string>

#include "config.h"
#include "date.h"

#ifdef ENABLE_TRACE_LOG
#  define TRACE_LOG(rank, msg) fmt::print("{} {} {}\n", iso8601_gmt(), rank, msg)
#else
#  define TRACE_LOG(rank, msg)
#endif

inline auto iso8601_gmt() -> std::string {
  return date::format("%FT%TZ",
                      date::floor<std::chrono::milliseconds>(std::chrono::system_clock::now()));
}
