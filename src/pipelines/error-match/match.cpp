#include "match.hpp"

#include <stdexcept>
#include <utility>

#include "simplify/convex_shadow.hpp"

namespace ssk::pipelines {
namespace {

const json::Value& field(const json::Value& doc, const char* key) {
  const json::Value* value = doc.find(key);
  if (value == nullptr) {
    throw std::runtime_error(std::string("error-match: baseline document has no '") + key + "'");
  }
  return *value;
}

// `bound` and `within_bound` are what viz/plot.py reports.
json::Value frechet_to_json(const Match& match, double eps) {
  return json::Object{
      {"distance", json::Value(match.measured.distance)},
      {"tolerance", json::Value(match.measured.tolerance)},
      {"bound", json::Value((1.0 + eps) * match.delta)},
      {"within_bound", json::Value(match.within_bound)},
      {"computer", json::Value("dv-gis-cup-2017")},
  };
}

}  // namespace

Match match(const frechet::Frechet<2>& computer, const simplify::Curve<2>& input,
            json::Value baseline, double eps, double relative_tol) {
  Match out;
  out.delta = field(field(baseline, "frechet"), "distance").num() / (1.0 + eps);
  out.simplified = simplify::ConvexShadow({{"delta", out.delta}, {"eps", eps}}).simplify(input);
  out.measured = measure(computer, input, out.simplified, relative_tol);
  out.within_bound = computer.within(input, out.simplified, (1.0 + eps) * out.delta);
  out.baseline = std::move(baseline);
  return out;
}

json::Value match_to_json(const io::Document& input, const std::string& dataset,
                          const std::string& source, double eps, const Match& match) {
  io::Trajectory points;
  points.reserve(match.simplified.size());
  for (const auto& q : match.simplified) {
    points.push_back({q[0], q[1]});
  }

  json::Object matched;
  for (const char* key : {"algorithm", "mode", "params", "stats", "frechet", "simplified"}) {
    if (const json::Value* value = match.baseline.find(key)) {
      matched[key] = *value;
    }
  }

  const std::size_t n = input.points.size();
  const std::size_t k = points.size();
  json::Object out{
      {"dim", json::Value(input.dim)},
      {"points", io::points_to_json(points)},
      {"algorithm", json::Value("convex-shadow")},
      {"mode", json::Value("error-matched")},
      {"params", json::Object{{"delta", json::Value(match.delta)}, {"eps", json::Value(eps)}}},
      {"stats", json::Object{{"input_size", json::Value(n)},
                             {"output_size", json::Value(k)},
                             {"compression", json::Value(static_cast<double>(n) /
                                                         static_cast<double>(k))}}},
      {"frechet", frechet_to_json(match, eps)},
      {"matched", json::Value(std::move(matched))},
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
  };
  if (!input.name.empty()) {
    out["name"] = json::Value(input.name);
  }
  return json::Value(std::move(out));
}

json::Value summary_to_json(const io::Document& input, const std::string& dataset,
                            const std::string& source, double eps,
                            const std::vector<Match>& matches) {
  json::Array pairs;
  for (const Match& match : matches) {
    const json::Value& stats = field(match.baseline, "stats");
    const json::Value& frechet = field(match.baseline, "frechet");
    pairs.emplace_back(json::Object{
        {"algorithm", field(match.baseline, "algorithm")},
        {"m", field(stats, "m")},
        {"baseline", json::Object{{"output_size", field(stats, "output_size")},
                                  {"target", field(stats, "target")},
                                  {"algorithm_runs", field(stats, "algorithm_runs")},
                                  {"distance", field(frechet, "distance")},
                                  {"tolerance", field(frechet, "tolerance")}}},
        {"ours", json::Object{{"delta", json::Value(match.delta)},
                              {"output_size", json::Value(match.simplified.size())},
                              {"distance", json::Value(match.measured.distance)},
                              {"tolerance", json::Value(match.measured.tolerance)},
                              {"within_bound", json::Value(match.within_bound)}}},
    });
  }
  return json::Object{
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
      {"name", json::Value(input.name)},
      {"eps", json::Value(eps)},
      {"input_size", json::Value(input.points.size())},
      {"pairs", json::Value(std::move(pairs))},
  };
}

}  // namespace ssk::pipelines
