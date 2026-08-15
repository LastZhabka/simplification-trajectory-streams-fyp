// Reading and writing curves in the project's trajectory format:
//
//     {"dim": 2, "name": "geolife/000/20081023025304", "points": [[x, y], ...]}
//
// Dimension-agnostic on purpose: this is the one layer every algorithm track shares, and
// the same document `src/trajio` exports and `src/viz` draws.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "io/json.hpp"

namespace ssk::io {

using Point = std::vector<double>;
using Trajectory = std::vector<Point>;

struct Document {
  std::size_t dim = 0;
  std::string name;  // provenance; empty when the file carries none
  Trajectory points;
};

// `dim` is taken from the field when present and inferred from the first point otherwise;
// either way every point must match it. Throws std::runtime_error on a document without a
// "points" array, on a point that is not an array of numbers, and on a dimension mismatch.
Document read_trajectory(const json::Value& doc);

Document read_trajectory_file(const std::string& path);

// Points alone, as an array of arrays -- for building a result document around them.
json::Value points_to_json(const Trajectory& points);

json::Value to_json(const Document& doc);

void write_trajectory_file(const std::string& path, const Document& doc);

}  // namespace ssk::io
