#include "geometry/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <vector>

namespace ssk::geometry {
namespace {

// Half a turn or more from a to b, counter-clockwise. Exactly half rounds either way, so the
// threshold sits just above it.
bool half_turn_apart(Vec a, Vec b) { return cross(a, b) <= 1e-9; }

bool before(Vec a, Vec b) { return a.x < b.x || (a.x == b.x && a.y < b.y); }

}  // namespace

ConvexPolygon regular_polygon(std::size_t m, double circumradius) {
  ConvexPolygon out;
  for (std::size_t k = 0; k < m; ++k) {
    const double angle = 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(m);
    out.normals.push_back({std::cos(angle), std::sin(angle)});
    out.offsets.push_back(circumradius * std::cos(std::numbers::pi / static_cast<double>(m)));
  }
  return out;
}

ConvexPolygon translated(const ConvexPolygon& poly, Vec v) {
  ConvexPolygon out = poly;
  for (std::size_t k = 0; k < out.offsets.size(); ++k) {
    out.offsets[k] = poly.offsets[k] + dot(poly.normals[k], v);
  }
  return out;
}

ConvexPolygon from_corners(const Polygon& corners) {
  ConvexPolygon out;
  for (std::size_t i = 0; i < corners.size(); ++i) {
    const Vec a = corners[i];
    const Vec edge = corners[(i + 1) % corners.size()] - a;
    const Vec n = Vec{edge.y, -edge.x} * (1.0 / norm(edge));
    out.normals.push_back(n);
    out.offsets.push_back(dot(n, a));
  }
  return out;
}

Vec meet(Vec n1, double b1, Vec n2, double b2) {
  const double d = cross(n1, n2);
  return {(b1 * n2.y - b2 * n1.y) / d, (n1.x * b2 - n2.x * b1) / d};
}

std::vector<std::size_t> edges(const ConvexPolygon& poly) {
  const std::vector<Vec>& n = poly.normals;
  const std::vector<double>& b = poly.offsets;
  const auto cuts = [&](std::size_t i, std::size_t j, std::size_t k) {
    return dot(n[k], meet(n[i], b[i], n[j], b[j])) > b[k];
  };
  std::deque<std::size_t> q;
  for (std::size_t k = 0; k < n.size(); ++k) {
    while (q.size() >= 2 && cuts(q[q.size() - 2], q.back(), k)) {
      q.pop_back();
    }
    while (q.size() >= 2 && cuts(q[0], q[1], k)) {
      q.pop_front();
    }
    if (!q.empty() && half_turn_apart(n[q.back()], n[k])) {
      return {};
    }
    q.push_back(k);
  }
  while (q.size() >= 3 && cuts(q[q.size() - 2], q.back(), q.front())) {
    q.pop_back();
  }
  while (q.size() >= 3 && cuts(q[0], q[1], q.back())) {
    q.pop_front();
  }
  if (q.size() < 3 || half_turn_apart(n[q.back()], n[q.front()])) {
    return {};
  }
  return {q.begin(), q.end()};
}

Base base(const ConvexPolygon& poly, Vec p) {
  const std::vector<Vec>& n = poly.normals;
  const std::vector<double>& b = poly.offsets;
  const std::vector<std::size_t> lines = edges(poly);
  const auto faces = [&](std::size_t k) { return dot(n[k], p) > b[k]; };
  Base out;
  std::size_t facing = 0;
  for (std::size_t j = 0; j < lines.size(); ++j) {
    const std::size_t here = lines[j];
    const std::size_t next = lines[(j + 1) % lines.size()];
    out.corners.push_back(meet(n[here], b[here], n[next], b[next]));
    facing += faces(here) ? 1 : 0;
    if (faces(here) && !faces(next)) {
      out.right = out.corners.back();
    }
    if (!faces(here) && faces(next)) {
      out.left = out.corners.back();
    }
  }
  double area = 0.0;
  for (std::size_t j = 1; j + 1 < out.corners.size(); ++j) {
    area += cross(out.corners[j] - out.corners[0], out.corners[j + 1] - out.corners[0]);
  }
  if (area <= 0.0 || facing == lines.size()) {
    out.corners.clear();
  }
  out.contains = facing == 0;
  return out;
}

Wedge cone(const Base& b, Vec p) {
  if (b.contains) {
    return Wedge{};
  }
  return {false, b.right - p, b.left - p};
}

bool within(const Wedge& w, Vec d) { return cross(w.right, d) >= 0.0 && cross(d, w.left) >= 0.0; }

std::optional<Wedge> intersect(const Wedge& a, const Wedge& b) {
  if (b.full) {
    return a;
  }
  Wedge out = b;
  if (!a.full) {
    out.right = within(a, b.right) ? b.right : a.right;
    out.left = within(a, b.left) ? b.left : a.left;
    if (!within(b, out.right) || !within(b, out.left)) {
      return std::nullopt;
    }
  }
  if (cross(out.right, out.left) <= 0.0) {
    return std::nullopt;
  }
  return out;
}

Polygon clip_left_of(const Polygon& poly, Vec p, Vec d) {
  const Vec n{d.y, -d.x};
  const double c = dot(n, p);
  Polygon out;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec cur = poly[i];
    const Vec nxt = poly[(i + 1) % poly.size()];
    const double sc = c - dot(n, cur);
    const double sn = c - dot(n, nxt);
    if (sc >= 0.0) {
      out.push_back(cur);
    }
    if ((sc > 0.0 && sn < 0.0) || (sc < 0.0 && sn > 0.0)) {
      out.push_back(cur + (nxt - cur) * (sc / (sc - sn)));
    }
  }
  return out;
}

Polygon cut(Polygon poly, Vec p, const Wedge& w) {
  if (!w.full) {
    poly = clip_left_of(poly, p, w.right);
    poly = clip_left_of(poly, p, w.left * -1.0);
  }
  return poly;
}

Vec closest(const Polygon& poly, Vec q) {
  bool inside = true;
  Vec best = poly[0];
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec ab = poly[(i + 1) % poly.size()] - a;
    const double length = norm(ab);
    const double reach = norm(q - a) + norm(q - a - ab);
    inside = inside && (length < 1e-9 * reach || cross(ab, q - a) >= 1e-9 * length * reach);
    const double t = length > 0.0 ? std::clamp(dot(q - a, ab) / dot(ab, ab), 0.0, 1.0) : 0.0;
    const Vec on = a + ab * t;
    if (dot(q - on, q - on) < dot(q - best, q - best)) {
      best = on;
    }
  }
  return inside ? q : best;
}

Polygon full_convex_hull(Polygon points, bool sorted) {
  if (!sorted) {
    std::sort(points.begin(), points.end(), before);
  }
  points.erase(std::unique(points.begin(), points.end()), points.end());
  if (points.size() < 3) {
    return points;
  }
  Polygon hull(2 * points.size());
  std::size_t k = 0;
  for (const Vec q : points) {
    while (k >= 2 && cross(hull[k - 1] - hull[k - 2], q - hull[k - 2]) <= 0.0) {
      --k;
    }
    hull[k++] = q;
  }
  const std::size_t lower = k + 1;
  for (std::size_t i = points.size() - 1; i-- > 0;) {
    while (k >= lower && cross(hull[k - 1] - hull[k - 2], points[i] - hull[k - 2]) <= 0.0) {
      --k;
    }
    hull[k++] = points[i];
  }
  hull.resize(k - 1);
  if (hull.size() >= 3 && !is_convex(hull)) {
    throw std::logic_error("full_convex_hull: the hull built is not convex");
  }
  return hull;
}

bool is_convex(const Polygon& poly) {
  if (poly.size() < 3) {
    return false;
  }
  double turning = 0.0;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec b = poly[(i + 1) % poly.size()];
    const Vec c = poly[(i + 2) % poly.size()];
    const double left = cross(b - a, c - a);
    if (left <= 0.0) {
      return false;
    }
    turning += std::atan2(left, dot(b - a, c - b));
  }
  return turning < 3.0 * std::numbers::pi;
}

}  // namespace ssk::geometry
