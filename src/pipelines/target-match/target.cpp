#include "target.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "match.hpp"
#include "simplify/dots.hpp"
#include "simplify/douglas_peucker.hpp"
#include "simplify/squish.hpp"
#include "simplify/squish_e.hpp"

namespace ssk::pipelines {
namespace {

using simplify::Params;
using Indices = std::vector<std::size_t>;

std::vector<std::size_t> all(std::size_t n) {
  Indices out(n);
  std::iota(out.begin(), out.end(), std::size_t{0});
  return out;
}

simplify::Curve<2> at(const simplify::Curve<2>& curve, const Indices& indices) {
  simplify::Curve<2> out;
  out.reserve(indices.size());
  for (const std::size_t i : indices) {
    out.push_back(curve[i]);
  }
  return out;
}

// A budget k in [floor, n], n being the input itself. The midpoint is geometric: the targets
// sit at 4x to 64x compression, and a decision costs n * k, so the large budgets are the ones
// to visit least.
std::vector<Best> search_budget(const std::string& name, std::size_t floor,
                                const std::function<Indices(std::size_t)>& run,
                                const frechet::Frechet<2>& computer,
                                const simplify::Curve<2>& curve,
                                const std::vector<Target>& targets) {
  const std::size_t n = curve.size();
  std::map<std::size_t, Indices> seen;
  int runs = 0;
  const auto indices = [&](std::size_t k) -> const Indices& {
    auto it = seen.find(k);
    if (it == seen.end()) {
      ++runs;
      it = seen.emplace(k, k >= n ? all(n) : run(k)).first;
    }
    return it->second;
  };

  std::vector<Best> out;
  std::size_t hi = n;
  for (const Target& target : targets) {
    runs = 0;
    const auto within = [&](std::size_t k) {
      return k >= n || computer.within(curve, at(curve, indices(k)), target.distance);
    };
    // Every budget at or below lo is taken as not within: under the floor it cannot be run.
    std::size_t lo = std::min(floor - 1, hi - 1);
    while (hi - lo > 1) {
      const auto mid = std::clamp(
          static_cast<std::size_t>(std::llround(
              std::sqrt(static_cast<double>(std::max<std::size_t>(lo, 1)) *
                        static_cast<double>(hi)))),
          lo + 1, hi - 1);
      if (within(mid)) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    out.push_back({target, {{name, static_cast<double>(hi)}}, indices(hi), runs, hi >= n});
  }
  return out;
}

// DOTS bounds error, not size, so its search is on log(lssd_threshold): the output shrinks as
// the threshold grows. It stops once the two ends give the same size, after kStable halvings
// that each landed on a size already at that end, or at kTight, as in sweep_dots.
std::vector<Best> search_dots(const frechet::Frechet<2>& computer, const simplify::Curve<2>& curve,
                              const simplify::Context& in, const std::vector<Target>& targets) {
  constexpr double kTight = 1e-6;
  constexpr int kStable = 3;
  const double floor_log = std::log(1e-22);
  const double cap_log = std::log(1e300);

  // The low end stands for the input itself, within every target, and is never run: DOTS at a
  // threshold near 0 did not finish in 10 minutes on geolife-001523's 92 645 points.
  std::map<double, Indices> seen{{floor_log, all(curve.size())}};
  int runs = 0;
  const auto indices = [&](double log_th) -> const Indices& {
    auto it = seen.find(log_th);
    if (it == seen.end()) {
      ++runs;
      it = seen.emplace(log_th, simplify::Dots(Params{{"lssd_threshold", std::exp(log_th)}})
                                    .indices(curve, in))
               .first;
    }
    return it->second;
  };

  std::vector<Best> out;
  double a = floor_log;  // within every target so far
  for (const Target& target : targets) {
    runs = 0;
    const auto within = [&](double log_th) {
      return log_th == floor_log ||
             computer.within(curve, at(curve, indices(log_th)), target.distance);
    };
    // A shortcut over L points that strays up to T sums about L * T^2, so n * T^2 is near the
    // answer. Starting high instead, at 1e4 as sweep_dots does, cost 900 s in one run on
    // geolife-001523: DOTS is slowest where it collapses a long curve to a few points.
    const double guess = std::log(static_cast<double>(curve.size()) * target.distance *
                                  target.distance);
    double b = std::max(guess, a + std::log(4.0));
    while (b < cap_log && indices(b).size() > 2 && within(b)) {
      a = b;
      b += std::log(100.0);
    }
    if (b < cap_log && within(b)) {
      a = b;
    } else {
      int stable = 0;
      while (b - a > kTight && stable < kStable && indices(a).size() != indices(b).size()) {
        const double mid = 0.5 * (a + b);
        const std::size_t size = indices(mid).size();
        if (within(mid)) {
          stable = (size == indices(a).size()) ? stable + 1 : 0;
          a = mid;
        } else {
          stable = (size == indices(b).size()) ? stable + 1 : 0;
          b = mid;
        }
      }
    }
    const bool input = a == floor_log;
    out.push_back({target, {{"lssd_threshold", input ? 0.0 : std::exp(a)}}, indices(a), runs, input});
  }
  return out;
}

std::string level_name(double fraction) {
  std::ostringstream out;
  out << "f" << fraction;
  return out.str();
}

json::Value target_json(const Target& t) {
  return json::Object{{"level", json::Value(t.level)},
                      {"fraction", json::Value(t.fraction)},
                      {"diagonal", json::Value(t.diagonal)},
                      {"distance", json::Value(t.distance)}};
}

}  // namespace

std::vector<Target> targets(double diagonal) {
  std::vector<Target> out;
  if (diagonal <= 0.0) {
    return out;
  }
  for (const double f : kFractions) {
    out.push_back({level_name(f), f, diagonal, f * diagonal});
  }
  return out;
}

std::vector<Best> search(const std::string& algorithm, const frechet::Frechet<2>& computer,
                         const simplify::Curve<2>& curve, const simplify::Context& in,
                         const std::vector<Target>& targets) {
  const std::size_t n = curve.size();
  if (algorithm == "douglas-peucker-n") {
    return search_budget(
        "count", 2,
        [&](std::size_t k) {
          return simplify::DouglasPeucker<2>(Params{{"count", static_cast<double>(k)}})
              .indices(curve);
        },
        computer, curve, targets);
  }
  if (algorithm == "squish") {
    return search_budget(
        "buffer_size", 5,
        [&](std::size_t k) {
          return simplify::Squish(Params{{"buffer_size", static_cast<double>(k)}})
              .indices(curve, in);
        },
        computer, curve, targets);
  }
  if (algorithm == "squish-e") {
    // It keeps about (n - 1) / cmp_ratio points, so a budget k asks for the ratio that keeps k.
    return search_budget(
        "budget", 3,
        [&](std::size_t k) {
          const double ratio = std::ceil(static_cast<double>(n - 1) / static_cast<double>(k));
          return simplify::SquishE(Params{{"cmp_ratio", ratio}, {"sed_error", 0.0}})
              .indices(curve, in);
        },
        computer, curve, targets);
  }
  if (algorithm == "dots") {
    return search_dots(computer, curve, in, targets);
  }
  throw std::runtime_error("target-search: no baseline named '" + algorithm + "'");
}

json::Value best_to_json(const io::Document& doc, const std::string& dataset,
                         const std::string& source, const std::string& algorithm,
                         const Best& best, const Measurement& measured) {
  io::Trajectory points;
  points.reserve(best.indices.size());
  for (const std::size_t i : best.indices) {
    points.push_back(doc.points[i]);
  }
  json::Object params;
  for (const auto& [key, value] : best.params) {
    params[key] = json::Value(value);
  }
  const std::size_t n = doc.points.size();
  const std::size_t k = best.indices.size();
  json::Object out{
      {"dim", json::Value(doc.dim)},
      {"points", io::points_to_json(points)},
      {"algorithm", json::Value(algorithm)},
      {"mode", json::Value("error-target")},
      {"params", json::Value(std::move(params))},
      {"stats", json::Object{{"input_size", json::Value(n)},
                             {"output_size", json::Value(k)},
                             {"compression", json::Value(static_cast<double>(n) /
                                                         static_cast<double>(k))},
                             {"algorithm_runs", json::Value(best.algorithm_runs)},
                             {"input", json::Value(best.input)}}},
      {"target", target_json(best.target)},
      {"frechet", json::Object{{"distance", json::Value(measured.distance)},
                               {"tolerance", json::Value(measured.tolerance)},
                               {"computer", json::Value("dv-gis-cup-2017")}}},
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
  };
  if (!doc.name.empty()) {
    out["name"] = json::Value(doc.name);
  }
  if (!doc.t.empty()) {
    json::Array t;
    for (const std::size_t i : best.indices) {
      t.emplace_back(doc.t[i]);
    }
    out["t"] = json::Value(std::move(t));
    out["t_unit"] = json::Value(doc.t_unit);
  }
  return json::Value(std::move(out));
}

json::Value search_summary_to_json(const io::Document& doc, const std::string& dataset,
                                   const std::string& source, const std::vector<Target>& targets,
                                   const std::vector<Found>& found) {
  json::Array levels;
  for (const Target& t : targets) {
    levels.emplace_back(target_json(t));
  }
  json::Array results;
  for (const Found& f : found) {
    results.emplace_back(json::Object{
        {"algorithm", json::Value(f.algorithm)},
        {"level", json::Value(f.best.target.level)},
        {"output_size", json::Value(f.best.indices.size())},
        {"distance", json::Value(f.measured.distance)},
        {"tolerance", json::Value(f.measured.tolerance)},
        {"algorithm_runs", json::Value(f.best.algorithm_runs)},
        {"input", json::Value(f.best.input)},
    });
  }
  return json::Object{
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
      {"name", json::Value(doc.name)},
      {"input_size", json::Value(doc.points.size())},
      {"targets", json::Value(std::move(levels))},
      {"results", json::Value(std::move(results))},
  };
}

std::vector<Target> targets_of(const json::Value& search_summary) {
  std::vector<Target> out;
  for (const json::Value& t : search_summary.find("targets")->arr()) {
    out.push_back({t.find("level")->str(), t.find("fraction")->num(), t.find("diagonal")->num(),
                   t.find("distance")->num()});
  }
  return out;
}

OursRun run_ours(const frechet::Frechet<2>& computer, const std::string& ours,
                 const simplify::Curve<2>& curve, const Target& target, double eps,
                 bool unscaled) {
  OursRun out;
  out.target = target;
  out.delta = unscaled ? target.distance : target.distance / (1.0 + eps);
  out.simplified = make_ours(ours, {{"delta", out.delta}, {"eps", eps}})->simplify(curve);
  out.measured = measure(computer, curve, out.simplified, 1e-9);
  out.within = computer.within(curve, out.simplified, (1.0 + eps) * out.delta);
  out.within_target = computer.within(curve, out.simplified, target.distance);
  return out;
}

json::Value ours_to_json(const io::Document& doc, const std::string& dataset,
                         const std::string& source, const std::string& ours, double eps,
                         bool unscaled, const OursRun& run) {
  io::Trajectory points;
  points.reserve(run.simplified.size());
  for (const auto& q : run.simplified) {
    points.push_back({q[0], q[1]});
  }
  const std::size_t n = doc.points.size();
  const std::size_t k = points.size();
  json::Object out{
      {"dim", json::Value(doc.dim)},
      {"points", io::points_to_json(points)},
      {"algorithm", json::Value(ours)},
      {"mode", json::Value("error-target")},
      {"params", json::Object{{"delta", json::Value(run.delta)}, {"eps", json::Value(eps)}}},
      {"delta_rule", json::Value(unscaled ? "T" : "T/(1+eps)")},
      {"stats", json::Object{{"input_size", json::Value(n)},
                             {"output_size", json::Value(k)},
                             {"compression", json::Value(static_cast<double>(n) /
                                                         static_cast<double>(k))}}},
      {"target", target_json(run.target)},
      {"frechet", json::Object{{"distance", json::Value(run.measured.distance)},
                               {"tolerance", json::Value(run.measured.tolerance)},
                               {"bound", json::Value((1.0 + eps) * run.delta)},
                               {"within_bound", json::Value(run.within)},
                               {"within_target", json::Value(run.within_target)},
                               {"computer", json::Value("dv-gis-cup-2017")}}},
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
  };
  if (!doc.name.empty()) {
    out["name"] = json::Value(doc.name);
  }
  return json::Value(std::move(out));
}

json::Value match_summary_to_json(const io::Document& doc, const std::string& dataset,
                                  const std::string& source, const std::string& ours, double eps,
                                  bool unscaled, const json::Value& search_summary,
                                  const std::vector<OursRun>& runs) {
  json::Array pairs;
  for (const OursRun& run : runs) {
    for (const json::Value& r : search_summary.find("results")->arr()) {
      if (r.find("level")->str() != run.target.level) {
        continue;
      }
      pairs.emplace_back(json::Object{
          {"algorithm", *r.find("algorithm")},
          {"level", json::Value(run.target.level)},
          {"fraction", json::Value(run.target.fraction)},
          {"target", json::Value(run.target.distance)},
          {"baseline", json::Object{{"output_size", *r.find("output_size")},
                                    {"distance", *r.find("distance")},
                                    {"tolerance", *r.find("tolerance")},
                                    {"algorithm_runs", *r.find("algorithm_runs")},
                                    {"input", *r.find("input")}}},
          {"ours", json::Object{{"delta", json::Value(run.delta)},
                                {"output_size", json::Value(run.simplified.size())},
                                {"distance", json::Value(run.measured.distance)},
                                {"tolerance", json::Value(run.measured.tolerance)},
                                {"within_bound", json::Value(run.within)},
                                {"within_target", json::Value(run.within_target)}}},
      });
    }
  }
  return json::Object{
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
      {"name", json::Value(doc.name)},
      {"ours", json::Value(ours)},
      {"eps", json::Value(eps)},
      {"delta_rule", json::Value(unscaled ? "T" : "T/(1+eps)")},
      {"input_size", json::Value(doc.points.size())},
      {"pairs", json::Value(std::move(pairs))},
  };
}

}  // namespace ssk::pipelines
