#include "measure.hpp"

#include <algorithm>
#include <cmath>

namespace ssk::pipelines {

double diagonal(const frechet::Curve<2>& curve) {
  if (curve.empty()) {
    return 0.0;
  }
  double min_x = curve[0][0];
  double max_x = curve[0][0];
  double min_y = curve[0][1];
  double max_y = curve[0][1];
  for (const auto& p : curve) {
    min_x = std::min(min_x, p[0]);
    max_x = std::max(max_x, p[0]);
    min_y = std::min(min_y, p[1]);
    max_y = std::max(max_y, p[1]);
  }
  return std::hypot(max_x - min_x, max_y - min_y);
}

Measurement measure(const frechet::Frechet<2>& computer, const frechet::Curve<2>& input,
                    const frechet::Curve<2>& simplified, double relative_tol) {
  Measurement m;
  m.input_points = input.size();
  m.output_points = simplified.size();
  const double scale = diagonal(input);
  // A curve with no extent still needs a positive tolerance for the bisection to terminate.
  m.tolerance = (scale > 0.0) ? scale * relative_tol : relative_tol;
  m.distance = computer.distance(input, simplified, m.tolerance);
  return m;
}

json::Value to_json(const io::Document& input, const json::Value& simplified,
                    const std::string& dataset, const std::string& source,
                    const std::string& simplified_path, const Measurement& measurement) {
  json::Object frechet{
      {"distance", json::Value(measurement.distance)},
      {"tolerance", json::Value(measurement.tolerance)},
      {"computer", json::Value("dv-gis-cup-2017")},
  };

  json::Object out{
      {"dim", json::Value(input.dim)},
      {"dataset", json::Value(dataset)},
      {"source", json::Value(source)},
      {"simplified", json::Value(simplified_path)},
      {"frechet", json::Value(std::move(frechet))},
  };
  if (!input.name.empty()) {
    out["name"] = json::Value(input.name);
  }
  // Carried through from the simplification, so a results table needs only this document.
  for (const char* key : {"algorithm", "mode", "params", "stats"}) {
    if (const json::Value* v = simplified.find(key)) {
      out[key] = *v;
    }
  }
  return json::Value(std::move(out));
}

}  // namespace ssk::pipelines
