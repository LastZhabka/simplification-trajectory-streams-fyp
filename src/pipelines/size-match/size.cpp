#include "size.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>

#include "match.hpp"

namespace ssk::pipelines {
namespace {

// ssk_frechet's, so ours' delta is found as finely as both errors are measured.
constexpr double kRelativeTol = 1e-9;

const json::Value& field(const json::Value& doc, const char* key) {
  const json::Value* value = doc.find(key);
  if (value == nullptr) {
    throw std::runtime_error(std::string("size-match: baseline document has no '") + key + "'");
  }
  return *value;
}

}  // namespace

std::size_t budget_of(const json::Value& baseline) {
  return static_cast<std::size_t>(field(field(baseline, "stats"), "output_size").num());
}

std::vector<Fit> fit(const frechet::Frechet<2>& computer, const std::string& ours,
                     const simplify::Curve<2>& curve, double eps,
                     const std::vector<std::size_t>& budgets) {
  const double scale = diagonal(curve);
  const double tol = (scale > 0.0) ? scale * kRelativeTol : kRelativeTol;
  const auto run = [&](double delta) {
    return make_ours(ours, {{"delta", delta}, {"eps", eps}})->simplify(curve);
  };
  // Sizes only: on the longest trajectory a run under the budget is tens of thousands of points.
  std::map<double, std::size_t> seen;
  int runs = 0;
  const auto size_at = [&](double delta) {
    auto [it, fresh] = seen.try_emplace(delta, 0);
    if (fresh) {
      ++runs;
      it->second = run(delta).size();
    }
    return it->second;
  };

  std::map<std::size_t, Fit> done;
  std::vector<Fit> out;
  for (const std::size_t budget : budgets) {
    auto found = done.find(budget);
    if (found == done.end()) {
      if (budget < 2) {
        throw std::runtime_error("size-match: ours cannot fit in fewer than 2 vertices");
      }
      runs = 0;
      // Every delta at or below lo is taken as over the budget: at 0 ours cannot be run. At the
      // diagonal of the curve's box every segment in the box is within delta, so kappa is 2 and
      // ours has at most 2 vertices; the doubling is for the margin run() takes off delta.
      double lo = 0.0;
      double hi = std::max(scale, tol);
      while (size_at(hi) > budget) {
        lo = hi;
        hi *= 2.0;
      }
      while (hi - lo > tol) {
        const double mid = 0.5 * (lo + hi);
        if (size_at(mid) <= budget) {
          hi = mid;
        } else {
          lo = mid;
        }
      }
      Fit f;
      f.budget = budget;
      f.delta = hi;
      f.below = lo;
      f.simplified = run(hi);
      f.measured = measure(computer, curve, f.simplified, kRelativeTol);
      f.within_bound = computer.within(curve, f.simplified, (1.0 + eps) * hi);
      f.runs = runs;
      found = done.emplace(budget, std::move(f)).first;
    }
    out.push_back(found->second);
  }
  return out;
}

json::Value fit_to_json(const io::Document& input, const std::string& dataset,
                        const std::string& source, const std::string& ours, double eps,
                        const Fit& fit, const json::Value& baseline) {
  io::Trajectory points;
  points.reserve(fit.simplified.size());
  for (const auto& q : fit.simplified) {
    points.push_back({q[0], q[1]});
  }
  json::Object matched;
  for (const char* key : {"algorithm", "mode", "params", "stats", "frechet", "simplified"}) {
    if (const json::Value* value = baseline.find(key)) {
      matched[key] = *value;
    }
  }
  const std::size_t n = input.points.size();
  const std::size_t k = points.size();
  json::Object out{
      {"dim", json::Value(input.dim)},
      {"points", io::points_to_json(points)},
      {"algorithm", json::Value(ours)},
      {"mode", json::Value("size-matched")},
      {"params", json::Object{{"delta", json::Value(fit.delta)}, {"eps", json::Value(eps)}}},
      {"delta_rule", json::Value("size")},
      {"stats", json::Object{{"input_size", json::Value(n)},
                             {"output_size", json::Value(k)},
                             {"budget", json::Value(fit.budget)},
                             {"compression", json::Value(static_cast<double>(n) /
                                                         static_cast<double>(k))},
                             {"algorithm_runs", json::Value(fit.runs)}}},
      {"frechet", json::Object{{"distance", json::Value(fit.measured.distance)},
                               {"tolerance", json::Value(fit.measured.tolerance)},
                               {"bound", json::Value((1.0 + eps) * fit.delta)},
                               {"within_bound", json::Value(fit.within_bound)},
                               {"computer", json::Value("dv-gis-cup-2017")}}},
      {"matched", json::Value(std::move(matched))},
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
  };
  if (!input.name.empty()) {
    out["name"] = json::Value(input.name);
  }
  return json::Value(std::move(out));
}

json::Value fit_summary_to_json(const io::Document& input, const std::string& dataset,
                                const std::string& source, const std::string& ours, double eps,
                                const std::vector<Fitted>& fitted) {
  json::Array pairs;
  for (const Fitted& f : fitted) {
    const json::Value& stats = field(f.baseline, "stats");
    const json::Value& frechet = field(f.baseline, "frechet");
    pairs.emplace_back(json::Object{
        {"algorithm", field(f.baseline, "algorithm")},
        {"m", field(stats, "m")},
        {"baseline", json::Object{{"output_size", field(stats, "output_size")},
                                  {"target", field(stats, "target")},
                                  {"algorithm_runs", field(stats, "algorithm_runs")},
                                  {"distance", field(frechet, "distance")},
                                  {"tolerance", field(frechet, "tolerance")}}},
        {"ours", json::Object{{"delta", json::Value(f.fit.delta)},
                              {"output_size", json::Value(f.fit.simplified.size())},
                              {"distance", json::Value(f.fit.measured.distance)},
                              {"tolerance", json::Value(f.fit.measured.tolerance)},
                              {"within_bound", json::Value(f.fit.within_bound)},
                              {"algorithm_runs", json::Value(f.fit.runs)}}},
    });
  }
  return json::Object{
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
      {"name", json::Value(input.name)},
      {"ours", json::Value(ours)},
      {"eps", json::Value(eps)},
      {"delta_rule", json::Value("size")},
      {"input_size", json::Value(input.points.size())},
      {"pairs", json::Value(std::move(pairs))},
  };
}

}  // namespace ssk::pipelines
