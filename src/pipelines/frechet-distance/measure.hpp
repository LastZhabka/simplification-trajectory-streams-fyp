// Measuring how far a simplification is from the trajectory it came from.
//
// The tolerance is relative, not absolute, because the coordinate unit differs by six orders
// of magnitude across the corpus -- Mopsi is degrees, NGSIM is State Plane feet, MOT is
// pixels -- so one absolute tolerance is either meaningless or ruinous depending on the
// dataset. It is taken against the input's bounding-box diagonal, which is in whatever unit
// the input happens to use.
#pragma once

#include <cstddef>
#include <string>

#include "frechet/frechet.hpp"
#include "io/json.hpp"
#include "io/trajectory.hpp"

namespace ssk::pipelines {

struct Measurement {
  double distance = 0.0;
  double tolerance = 0.0;
  std::size_t input_points = 0;
  std::size_t output_points = 0;
};

// Largest extent of the curve, as a scale for the tolerance.
[[nodiscard]] double diagonal(const frechet::Curve<2>& curve);

[[nodiscard]] Measurement measure(const frechet::Frechet<2>& computer,
                                  const frechet::Curve<2>& input,
                                  const frechet::Curve<2>& simplified, double relative_tol);

// The output document: the measurement, plus the identity of what was measured and the
// hyper-parameters that produced it, carried over from the simplification's own document so
// a results table needs nothing else.
[[nodiscard]] json::Value to_json(const io::Document& input, const json::Value& simplified,
                                  const std::string& dataset, const std::string& source,
                                  const std::string& simplified_path,
                                  const Measurement& measurement);

}  // namespace ssk::pipelines
