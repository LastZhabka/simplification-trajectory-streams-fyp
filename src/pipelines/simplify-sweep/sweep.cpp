#include "sweep.hpp"

#include <cmath>
#include <map>

#include "simplify/dots.hpp"
#include "simplify/douglas_peucker.hpp"
#include "simplify/squish.hpp"

namespace ssk::pipelines {
namespace {

using simplify::Params;

json::Array times_at(const std::vector<double>& t, const std::vector<std::size_t>& indices) {
  json::Array out;
  out.reserve(indices.size());
  for (const std::size_t i : indices) {
    out.emplace_back(t[i]);
  }
  return out;
}

}  // namespace

std::vector<std::size_t> budget_targets(std::size_t n, int max_m, std::size_t floor) {
  std::vector<std::size_t> out;
  for (int m = 1; m <= max_m; ++m) {
    const auto target = static_cast<std::size_t>(
        std::ceil(static_cast<double>(n) / std::pow(2.0, m)));
    if (target < floor || target >= n) {
      break;
    }
    out.push_back(target);
  }
  return out;
}

Sweep sweep_dpn(const Curve<2>& curve, int max_m) {
  Sweep sweep{"douglas-peucker-n", {}};
  int m = 0;
  for (const std::size_t target : budget_targets(curve.size(), max_m, 2)) {
    ++m;
    const auto count = static_cast<double>(target);
    sweep.runs.push_back({m,
                          target,
                          {{"count", count}},
                          simplify::DouglasPeucker<2>(Params{{"count", count}}).indices(curve),
                          1});
  }
  return sweep;
}

Sweep sweep_squish(const Curve<2>& curve, const Context& in, int max_m) {
  Sweep sweep{"squish", {}};
  int m = 0;
  for (const std::size_t target : budget_targets(curve.size(), max_m, 5)) {
    ++m;
    const auto size = static_cast<double>(target);
    sweep.runs.push_back(
        {m,
         target,
         {{"buffer_size", size}},
         simplify::Squish(Params{{"buffer_size", size}}).indices(curve, in),
         1});
  }
  return sweep;
}

// DOTS has no budget knob, so every rate is a search on lssd_threshold. Three things keep it
// cheap: the count is monotone in the threshold, so the search bisects; all rates are solved
// as one problem, smallest threshold first, each solved bound flooring the next; and every
// evaluation is memoised, with the bracket seeded by interpolating the counts already seen.
// Measured at 6.1 DOTS runs per budget, against 40 for an independent flat bisection.
//
// Not every budget is reachable. The output size is a step function of the threshold, so a
// target can fall in a gap between two attained sizes -- on geolife-013552 the m=1 budget is
// 32 242 while DOTS never emits more than 31 197, because it always takes some shortcut. Once
// the bracket's two ends stop moving, the search is only refining the location of a single
// step it has already found on both sides, and every further evaluation returns a size it has
// already seen. kStable stops that; without it the search ran the bracket down to kTight,
// which is meaningless precision for an integer-valued step function and cost hours on the
// longest trajectories.
Sweep sweep_dots(const Curve<2>& curve, const Context& in, int max_m) {
  constexpr double kFirstHi = 1e4;
  constexpr double kGrow = 1e6;      // widen the upper bound until it collapses the curve
  constexpr double kTight = 1e-6;    // bracket width in log space at which to give up
  constexpr int kStable = 3;         // halvings with neither end moving before giving up

  Sweep sweep{"dots", {}};
  int runs = 0;
  std::map<double, std::size_t> seen;  // log(threshold) -> output size
  const auto eval = [&](double log_th) {
    const auto it = seen.find(log_th);
    if (it != seen.end()) {
      return it->second;
    }
    ++runs;
    const std::size_t k =
        simplify::Dots(Params{{"lssd_threshold", std::exp(log_th)}}).indices(curve, in).size();
    seen.emplace(log_th, k);
    return k;
  };

  double lo = std::log(1e-22);
  int m = 0;
  for (const std::size_t target : budget_targets(curve.size(), max_m, 2)) {
    ++m;
    runs = 0;

    double hi = std::log(kFirstHi);
    while (eval(hi) > target && hi < std::log(1e300)) {
      hi += std::log(kGrow);
    }

    double a = lo;
    double b = hi;
    double best_log = hi;
    for (const auto& [log_th, k] : seen) {
      if (k > target && log_th > a) {
        a = log_th;
      }
      if (k <= target && log_th < b) {
        b = log_th;
        best_log = log_th;
      }
    }

    int stable = 0;
    while (b - a > kTight && stable < kStable) {
      double mid = 0.5 * (a + b);
      const auto ia = seen.find(a);
      const auto ib = seen.find(b);
      if (ia != seen.end() && ib != seen.end() && ia->second > ib->second) {
        const double f = static_cast<double>(ia->second - target) /
                         static_cast<double>(ia->second - ib->second);
        const double guess = a + f * (b - a);
        if (guess > a && guess < b) {
          mid = 0.5 * (mid + guess);
        }
      }
      const std::size_t k = eval(mid);
      const std::size_t was = (ib != seen.end()) ? ib->second : 0;
      if (k <= target) {
        best_log = mid;
        b = mid;
        if (k == target) {
          break;
        }
        stable = (k == was) ? stable + 1 : 0;
      } else {
        a = mid;
        stable = (ia != seen.end() && k == ia->second) ? stable + 1 : 0;
      }
    }

    const double threshold = std::exp(best_log);
    sweep.runs.push_back(
        {m,
         target,
         {{"lssd_threshold", threshold}},
         simplify::Dots(Params{{"lssd_threshold", threshold}}).indices(curve, in),
         runs + 1});
    lo = a;
  }
  return sweep;
}

bool ordered(const Run& run) {
  for (std::size_t i = 1; i < run.indices.size(); ++i) {
    if (run.indices[i - 1] >= run.indices[i]) {
      return false;
    }
  }
  return true;
}

json::Value run_to_json(const io::Document& doc, const std::string& source,
                        const std::string& algorithm, const Run& run) {
  io::Trajectory points;
  points.reserve(run.indices.size());
  for (const std::size_t i : run.indices) {
    points.push_back(doc.points[i]);
  }

  json::Object params;
  for (const auto& [key, value] : run.params) {
    params[key] = json::Value(value);
  }

  const auto kept = static_cast<double>(run.indices.size());
  json::Object stats{
      {"input_size", json::Value(doc.points.size())},
      {"output_size", json::Value(run.indices.size())},
      {"compression", json::Value(kept == 0.0 ? 0.0 : static_cast<double>(doc.points.size()) / kept)},
      {"m", json::Value(run.m)},
      {"rate", json::Value(std::pow(2.0, run.m))},
      {"target", json::Value(run.target)},
      {"algorithm_runs", json::Value(run.algorithm_runs)},
  };

  json::Object out{
      {"dim", json::Value(doc.dim)},
      {"points", io::points_to_json(points)},
      {"algorithm", json::Value(algorithm)},
      {"mode", json::Value("budget")},
      {"params", json::Value(std::move(params))},
      {"stats", json::Value(std::move(stats))},
      {"source", json::Value(source)},
  };
  if (!doc.name.empty()) {
    out["name"] = json::Value(doc.name);
  }
  if (!doc.t.empty()) {
    out["t"] = json::Value(times_at(doc.t, run.indices));
    out["t_unit"] = json::Value(doc.t_unit);
  }
  return json::Value(std::move(out));
}

}  // namespace ssk::pipelines
