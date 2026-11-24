#include "bench.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <utility>

namespace ssk::pipelines {

Timing time_runs(const simplify::Simplifier<2>& ours, const simplify::Curve<2>& curve) {
  Timing out;
  double total = 0.0;
  while (out.seconds.empty() || (total < kMinSeconds && out.seconds.size() < kMaxRepeats)) {
    const auto started = std::chrono::steady_clock::now();
    simplify::Curve<2> simplified = ours.simplify(curve);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    out.simplified = std::move(simplified);
    out.seconds.push_back(seconds);
    total += seconds;
  }
  return out;
}

std::string digest(const simplify::Curve<2>& curve) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const simplify::Point<2>& q : curve) {
    for (const double c : q) {
      std::array<unsigned char, sizeof(double)> bytes{};
      std::memcpy(bytes.data(), &c, sizeof(double));
      for (const unsigned char b : bytes) {
        hash ^= b;
        hash *= 1099511628211ULL;
      }
    }
  }
  std::array<char, 17> text{};
  std::snprintf(text.data(), text.size(), "%016llx", static_cast<unsigned long long>(hash));
  return text.data();
}

json::Value timing_to_json(const std::string& ours, double eps, double delta,
                           std::size_t input_size, const Timing& timing) {
  json::Array seconds;
  for (const double s : timing.seconds) {
    seconds.emplace_back(s);
  }
  return json::Object{
      {"ours", json::Value(ours)},
      {"eps", json::Value(eps)},
      {"delta", json::Value(delta)},
      {"input_size", json::Value(input_size)},
      {"output_size", json::Value(timing.simplified.size())},
      {"digest", json::Value(digest(timing.simplified))},
      {"seconds", json::Value(std::move(seconds))},
      {"best", json::Value(*std::min_element(timing.seconds.begin(), timing.seconds.end()))},
  };
}

}  // namespace ssk::pipelines
