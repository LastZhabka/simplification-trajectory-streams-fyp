// One implementation of ours on one trajectory at one delta and eps, timed. scripts/benchmark.py
// runs each such cell as its own process, so that a run past the time limit can be killed: a
// thread cannot be, and a deadline inside the stabbers would put the benchmark in the algorithm.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "io/json.hpp"
#include "simplify/simplifier.hpp"

namespace ssk::pipelines {

// A run of milliseconds is repeated, so it is not timed once against the clock's resolution.
inline constexpr double kMinSeconds = 1.0;
inline constexpr std::size_t kMaxRepeats = 5;

struct Timing {
  simplify::Curve<2> simplified;
  std::vector<double> seconds;  // one per repeat
};

// Repeats until the repeats add up to kMinSeconds, or there are kMaxRepeats of them.
[[nodiscard]] Timing time_runs(const simplify::Simplifier<2>& ours,
                               const simplify::Curve<2>& curve);

// FNV-1a over the coordinates' bytes, so two implementations meant to agree to the bit can be
// compared without writing out their points.
[[nodiscard]] std::string digest(const simplify::Curve<2>& curve);

// {"ours", "eps", "delta", "input_size", "output_size", "digest", "seconds": [...], "best"}
[[nodiscard]] json::Value timing_to_json(const std::string& ours, double eps, double delta,
                                         std::size_t input_size, const Timing& timing);

}  // namespace ssk::pipelines
