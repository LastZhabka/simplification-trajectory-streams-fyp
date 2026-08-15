#include "io/trajectory.hpp"

#include <fstream>
#include <stdexcept>
#include <utility>

namespace ssk::io {

Document read_trajectory(const json::Value& doc) {
  const json::Value* pts = doc.find("points");
  if (pts == nullptr || !pts->is_array()) {
    throw std::runtime_error("trajectory: document must have an array field \"points\"");
  }

  Document out;
  const json::Value* d = doc.find("dim");
  if (d != nullptr && d->is_number()) {
    out.dim = static_cast<std::size_t>(d->num());
  }
  const json::Value* name = doc.find("name");
  if (name != nullptr && name->is_string()) {
    out.name = name->str();
  }

  out.points.reserve(pts->arr().size());
  std::size_t index = 0;
  for (const json::Value& p : pts->arr()) {
    if (!p.is_array()) {
      throw std::runtime_error("trajectory: point " + std::to_string(index) +
                               " is not an array of coordinates");
    }
    Point q;
    q.reserve(p.arr().size());
    for (const json::Value& c : p.arr()) {
      if (!c.is_number()) {
        throw std::runtime_error("trajectory: point " + std::to_string(index) +
                                 " has a coordinate that is not a number");
      }
      q.push_back(c.num());
    }
    if (out.dim == 0) {
      out.dim = q.size();  // no "dim" field: infer from the first point
    }
    if (q.size() != out.dim) {
      throw std::runtime_error("trajectory: point " + std::to_string(index) + ": expected " +
                               std::to_string(out.dim) + " coordinates, found " +
                               std::to_string(q.size()));
    }
    out.points.push_back(std::move(q));
    ++index;
  }

  const json::Value* unit = doc.find("t_unit");
  if (unit != nullptr && unit->is_string()) {
    out.t_unit = unit->str();
  }
  const json::Value* times = doc.find("t");
  if (times != nullptr && !times->is_null()) {
    if (!times->is_array()) {
      throw std::runtime_error("trajectory: \"t\" must be an array");
    }
    if (times->arr().size() != out.points.size()) {
      throw std::runtime_error("trajectory: \"t\" has " +
                               std::to_string(times->arr().size()) + " entries but there are " +
                               std::to_string(out.points.size()) + " points");
    }
    out.t.reserve(times->arr().size());
    for (const json::Value& v : times->arr()) {
      if (!v.is_number()) {
        throw std::runtime_error("trajectory: \"t\" entry " +
                                 std::to_string(out.t.size()) + " is not a number");
      }
      out.t.push_back(v.num());
    }
  }
  return out;
}

Document read_trajectory_file(const std::string& path) {
  return read_trajectory(json::parse_file(path));
}

json::Value points_to_json(const Trajectory& points) {
  json::Array a;
  a.reserve(points.size());
  for (const Point& p : points) {
    json::Array c;
    c.reserve(p.size());
    for (const double v : p) {
      c.push_back(json::Value(v));
    }
    a.push_back(json::Value(std::move(c)));
  }
  return json::Value(std::move(a));
}

json::Value to_json(const Document& doc) {
  json::Object out;
  out.emplace("dim", json::Value(static_cast<double>(doc.dim)));
  if (!doc.name.empty()) {
    out.emplace("name", json::Value(doc.name));
  }
  if (!doc.t_unit.empty()) {
    out.emplace("t_unit", json::Value(doc.t_unit));
  }
  if (!doc.t.empty()) {
    json::Array times;
    times.reserve(doc.t.size());
    for (const double v : doc.t) {
      times.push_back(json::Value(v));
    }
    out.emplace("t", json::Value(std::move(times)));
  }
  out.emplace("points", points_to_json(doc.points));
  return json::Value(std::move(out));
}

void write_trajectory_file(const std::string& path, const Document& doc) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("trajectory: cannot write '" + path + "'");
  }
  file << to_json(doc).dump(2) << "\n";
}

}  // namespace ssk::io
