#pragma once

#include <chrono>
#include <string>

#include <cxxmpi/cart_comm.hpp>
#include <cxxmpi/comm.hpp>
#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include "rdbench/v2/options.hpp"
#include "rdbench/v2/utils/human_readable.hpp"

namespace rdbench::v2 {

namespace detail {

[[nodiscard]]
inline auto format_iso8601_with_timezone() -> std::string {
  auto now = std::chrono::system_clock::now();
  auto now_time_t = std::chrono::system_clock::to_time_t(now);
  auto now_tm = *std::localtime(&now_time_t);  // NOLINT not thread-safe

  // timezone offset
  std::time_t gmt = now_time_t;
  std::time_t local =
      std::mktime(std::localtime(&gmt));             // NOLINT not thread-safe
  std::time_t utc = std::mktime(std::gmtime(&gmt));  // NOLINT not thread-safe
  int offset_minutes = static_cast<int>(std::difftime(local, utc)) / 60;

  // ISO8601
  return fmt::format("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}{:+03d}:{:02}",
                     now_tm.tm_year + 1900, now_tm.tm_mon + 1, now_tm.tm_mday,
                     now_tm.tm_hour, now_tm.tm_min, now_tm.tm_sec,
                     offset_minutes / 60, std::abs(offset_minutes % 60));
}

}  // namespace detail

class metrics_handler {
 public:
  explicit metrics_handler(const options& opt,
                           const cxxmpi::weak_cart_comm comm)
      : opts_{opt},
        comm_{comm},
        start_time_{detail::format_iso8601_with_timezone()} {}

  [[nodiscard]]
  auto create_config_json() const -> nlohmann::ordered_json {
    auto j = nlohmann::ordered_json{};
    j["nprocs"] = comm_.size();
    j.merge_patch(nlohmann::ordered_json(opts_.get()));

    const auto dims = comm_.dims();
    j["nrTilesX"] = dims[1];
    j["nrTilesY"] = dims[0];

    const auto tile_size_bytes = this->tile_size_bytes();
    const auto file_size_bytes = this->file_size_bytes();
    const auto total_io_bytes = this->total_io_bytes();

    j["szTile"] = tile_size_bytes;
    j["szTileHuman"] = utils::to_human(tile_size_bytes) + "iB";
    j["szFile"] = file_size_bytes;
    j["szFileHuman"] = utils::to_human(file_size_bytes) + "iB";
    j["szTotalWrite"] = total_io_bytes;
    j["szTotalWriteHuman"] = utils::to_human(total_io_bytes) + "iB";

    return j;
  }

  [[nodiscard]]
  auto tile_size_bytes() const -> size_t {
    return opts_.get().sz_tile_x * opts_.get().sz_tile_y * sizeof(double);
  }

  [[nodiscard]]
  auto file_size_bytes() const -> size_t {
    const auto dims = comm_.dims();
    return tile_size_bytes() * static_cast<size_t>(dims[0])
         * static_cast<size_t>(dims[1]);
  }

  [[nodiscard]]
  auto total_io_bytes() const -> size_t {
    return file_size_bytes() * opts_.get().nr_files();
  }

  [[nodiscard]]
  auto create_timestamp_json() const -> nlohmann::ordered_json {
    auto j = nlohmann::ordered_json{};
    j["startTime"] = start_time_;
    j["endTime"] = detail::format_iso8601_with_timezone();
    return j;
  }

  [[nodiscard]]
  auto create_output_json(const nlohmann::ordered_json& prof_result) const
      -> nlohmann::ordered_json {
    auto result_json = create_timestamp_json();
    const auto config_json = create_config_json();

    result_json.merge_patch(config_json);
    result_json.merge_patch(prof_result);

    auto io_time = prof_result["IOTotalTime"].get<double>();
    auto io_bytes = this->total_io_bytes();
    auto bw = static_cast<double>(io_bytes) / io_time;
    result_json["writeBandwidth"] = bw;
    result_json["writeBandwidthHuman"] = utils::to_human(bw) + "iB/s";

    return result_json;
  }

 private:
  std::reference_wrapper<const options> opts_;
  cxxmpi::weak_cart_comm comm_;
  std::string start_time_;
};

}  // namespace rdbench::v2
