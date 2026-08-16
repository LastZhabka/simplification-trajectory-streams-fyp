// The sweep drives each baseline to ceil(N / 2^m): exactly, for the two that take a budget,
// and by a search on its threshold for DOTS.
#include "simplify/dots.hpp"
#include "sweep.hpp"
#include "test_main.hpp"

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace {

using ssk::pipelines::Run;
using ssk::pipelines::Sweep;
using ssk::simplify::Context;
using ssk::simplify::Curve;

Curve<2> wander(std::size_t n) {
  Curve<2> c;
  std::uint32_t s = 7654321;
  double x = 0.0;
  double y = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({x, y});
    s = s * 1664525u + 1013904223u;
    x += static_cast<double>((s >> 16) % 32u) / 16.0;
    s = s * 1664525u + 1013904223u;
    y += static_cast<double>((s >> 16) % 32u) / 16.0 - 1.0;
  }
  return c;
}

Context uniform_time(std::size_t n) {
  Context in;
  for (std::size_t i = 0; i < n; ++i) {
    in.t.push_back(static_cast<double>(i));
  }
  return in;
}

void check_runs(const Sweep& sweep, const Curve<2>& c, const std::string& who) {
  for (const Run& run : sweep.runs) {
    const auto& idx = run.indices;
    ::ssk::test::check(idx.size() >= 2, who + ": must keep both endpoints");
    ::ssk::test::check(idx.front() == 0 && idx.back() == c.size() - 1,
                       who + ": endpoints must survive");
    for (std::size_t i = 1; i < idx.size(); ++i) {
      ::ssk::test::check(idx[i - 1] < idx[i], who + ": indices must be increasing");
    }
    ::ssk::test::check(idx.back() < c.size(), who + ": index out of range");
    ::ssk::test::check(!run.params.empty(), who + ": a run must record its parameters");
  }
}

TEST("sweep/targets halve until they hit the floor") {
  using Sizes = std::vector<std::size_t>;
  CHECK(ssk::pipelines::budget_targets(100, 6, 2) == Sizes({50, 25, 13, 7, 4, 2}));
  CHECK(ssk::pipelines::budget_targets(100, 6, 5) == Sizes({50, 25, 13, 7}));
  CHECK(ssk::pipelines::budget_targets(100, 3, 2) == Sizes({50, 25, 13}));
  CHECK_MSG(ssk::pipelines::budget_targets(3, 6, 2) == Sizes({2}),
            "three points still support one rate: the two endpoints");
  CHECK_MSG(ssk::pipelines::budget_targets(2, 6, 2).empty(),
            "a curve too short for any rate contributes no runs");
  CHECK_MSG(ssk::pipelines::budget_targets(9, 6, 5) == Sizes({5}),
            "SQUISH's floor of 5 is reachable from 9 points, but only once");
  CHECK_MSG(ssk::pipelines::budget_targets(8, 6, 5).empty(),
            "halving 8 lands on 4, which is below SQUISH's floor");
}

TEST("sweep/dpn hits every budget exactly") {
  const Curve<2> c = wander(200);
  const Sweep sweep = ssk::pipelines::sweep_dpn(c, 6);
  check_runs(sweep, c, "dpn");
  CHECK(sweep.algorithm == "douglas-peucker-n");
  CHECK(sweep.runs.size() == 6);
  for (const Run& run : sweep.runs) {
    CHECK_MSG(run.indices.size() == run.target, "DPn delivers the budget exactly");
    CHECK(run.params.front().first == "count");
    CHECK_CLOSE(run.params.front().second, static_cast<double>(run.target), 1e-12);
  }
  CHECK(sweep.runs.front().m == 1);
  CHECK(sweep.runs.front().target == 100);
  CHECK(sweep.runs.back().target == 4);
}

TEST("sweep/squish stays within every budget and stops at its floor") {
  const Curve<2> c = wander(200);
  const Context in = uniform_time(c.size());
  const Sweep sweep = ssk::pipelines::sweep_squish(c, in, 6);
  check_runs(sweep, c, "squish");
  CHECK(sweep.algorithm == "squish");
  for (const Run& run : sweep.runs) {
    CHECK_MSG(run.indices.size() <= run.target, "output must stay within the buffer");
    CHECK(run.params.front().first == "buffer_size");
  }
  CHECK_MSG(sweep.runs.back().target >= 5, "no run may ask for a budget of 4 or less");
}

TEST("sweep/dots searches its threshold to the budget") {
  const Curve<2> c = wander(200);
  const Context in = uniform_time(c.size());
  const Sweep sweep = ssk::pipelines::sweep_dots(c, in, 6);
  check_runs(sweep, c, "dots");
  CHECK(sweep.algorithm == "dots");
  CHECK(sweep.runs.size() == 6);

  double previous = 0.0;
  for (const Run& run : sweep.runs) {
    CHECK_MSG(run.indices.size() <= run.target, "the search must not overshoot the budget");
    CHECK(run.params.front().first == "lssd_threshold");
    const double threshold = run.params.front().second;
    CHECK_MSG(threshold > previous, "a smaller budget needs a larger threshold");
    previous = threshold;
    CHECK_MSG(run.algorithm_runs > 1, "a searched budget takes more than one run");
  }
}

TEST("sweep/dots reaches the same budgets as an independent search") {
  const Curve<2> c = wander(150);
  const Context in = uniform_time(c.size());
  const Sweep sweep = ssk::pipelines::sweep_dots(c, in, 4);

  for (const Run& run : sweep.runs) {
    // The same budget, found by a plain bisection that shares nothing between rates.
    double lo = 1e-22;
    double hi = 1e4;
    std::size_t flat = 2;
    for (int it = 0; it < 80; ++it) {
      const double mid = std::exp(0.5 * (std::log(lo) + std::log(hi)));
      const std::size_t k =
          ssk::simplify::Dots(ssk::simplify::Params{{"lssd_threshold", mid}}).indices(c, in).size();
      if (k <= run.target) {
        flat = k;
        hi = mid;
      } else {
        lo = mid;
      }
    }
    CHECK_MSG(run.indices.size() == flat, "the swept search must find what a flat one finds");
  }
}

TEST("sweep/the document carries its parameters and rehydrates") {
  const auto doc = ssk::io::read_trajectory(ssk::json::parse(
      R"({"dim": 2, "name": "toy/1", "t": [0,1,2,3,4,5,6,7,8,9],
          "points": [[0,0],[1,4],[2,0],[3,5],[4,0],[5,6],[6,0],[7,3],[8,0],[9,1]]})"));
  const auto c = ssk::simplify::curve_of<2>(doc);
  const Sweep sweep = ssk::pipelines::sweep_dpn(c, 2);
  const auto written = ssk::pipelines::to_json(doc, "toy/toy-1.json", sweep);

  // Through a round trip, so the test covers what actually lands on disk.
  const auto read = ssk::json::parse(written.dump());
  CHECK(read.find("algorithm")->str() == "douglas-peucker-n");
  CHECK(read.find("source")->str() == "toy/toy-1.json");
  CHECK(read.find("name")->str() == "toy/1");
  CHECK(read.find("input_points")->num() == 10.0);
  CHECK(read.find("dim")->num() == 2.0);

  const auto& runs = read.find("runs")->arr();
  CHECK(runs.size() == sweep.runs.size());
  const auto& first = runs.front();
  CHECK(first.find("m")->num() == 1.0);
  CHECK(first.find("rate")->num() == 2.0);
  CHECK(first.find("target")->num() == 5.0);
  CHECK(first.find("kept")->num() == 5.0);
  CHECK_MSG(first.find("params")->find("count")->num() == 5.0,
            "the hyper-parameter used must travel with the result");

  const auto& indices = first.find("indices")->arr();
  CHECK(indices.size() == 5);
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const auto at = static_cast<std::size_t>(indices[i].num());
    ::ssk::test::check(at < doc.points.size(), "index must address the source document");
    ::ssk::test::check_close(doc.points[at][0], c[sweep.runs.front().indices[i]][0], 1e-12,
                             "an index must rehydrate to the point it named");
  }
}

TEST("sweep/a curve too short to halve produces no runs") {
  const Curve<2> tiny = wander(2);
  CHECK(ssk::pipelines::sweep_dpn(tiny, 6).runs.empty());
  CHECK(ssk::pipelines::sweep_squish(tiny, uniform_time(2), 6).runs.empty());
  CHECK(ssk::pipelines::sweep_dots(tiny, uniform_time(2), 6).runs.empty());

  // A short curve drops out of the deeper rates rather than being clamped into them, and
  // out of SQUISH entirely below its floor of 5.
  const Curve<2> ten = wander(10);
  CHECK(ssk::pipelines::sweep_dpn(ten, 6).runs.size() == 3);
  CHECK(ssk::pipelines::sweep_squish(ten, uniform_time(10), 6).runs.size() == 1);
}

}  // namespace
