// Timing ours for docs/benchmark.md, and the line ssk_bench prints.
#include "bench.hpp"
#include "simplify/greedy_stabber/greedy_stabber_decomp.hpp"
#include "test_main.hpp"

#include <cmath>
#include <string>

namespace {

using ssk::pipelines::digest;
using ssk::pipelines::time_runs;
using ssk::pipelines::timing_to_json;
using ssk::simplify::Curve;
using ssk::simplify::GreedyStabbingShared;

Curve<2> wave(std::size_t n) {
  Curve<2> c;
  for (std::size_t i = 0; i < n; ++i) {
    const auto t = static_cast<double>(i);
    c.push_back({t, std::sin(t / 10.0)});
  }
  return c;
}

TEST("bench/a fast run is repeated, and the output is the algorithm's") {
  const GreedyStabbingShared ours({{"delta", 0.1}, {"eps", 0.5}});
  const Curve<2> curve = wave(200);
  const auto timing = time_runs(ours, curve);
  CHECK_MSG(timing.seconds.size() == ssk::pipelines::kMaxRepeats,
            "a run of microseconds must be repeated up to the cap");
  CHECK(timing.simplified == ours.simplify(curve));
}

TEST("bench/the digest tells curves apart to the bit") {
  const Curve<2> a = wave(50);
  Curve<2> b = a;
  CHECK(digest(a) == digest(b));
  b[17][1] = std::nextafter(b[17][1], 2.0);
  CHECK(digest(a) != digest(b));
  CHECK(digest(a).size() == 16);
}

TEST("bench/the line carries the cell and the best repeat") {
  const GreedyStabbingShared ours({{"delta", 0.1}, {"eps", 0.5}});
  const Curve<2> curve = wave(200);
  const auto timing = time_runs(ours, curve);
  const auto doc = timing_to_json("greedy-stabbing-shared", 0.5, 0.1, curve.size(), timing);
  CHECK(doc.find("ours")->str() == "greedy-stabbing-shared");
  CHECK(doc.find("input_size")->num() == 200.0);
  CHECK(doc.find("output_size")->num() == static_cast<double>(timing.simplified.size()));
  CHECK(doc.find("seconds")->arr().size() == timing.seconds.size());
  double best = timing.seconds.front();
  for (const double s : timing.seconds) {
    best = std::min(best, s);
  }
  CHECK(doc.find("best")->num() == best);
}

}  // namespace
