// Plane geometry the greedy stabbers share: vectors, convex polygons as half-planes, cones from a
// point, clipping and hulls.
#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace ssk::geometry {

struct Vec {
  double x = 0.0;
  double y = 0.0;

  bool operator==(const Vec&) const = default;
};

using Polygon = std::vector<Vec>;

inline Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.y + b.y}; }
inline Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y}; }
inline Vec operator*(Vec a, double s) { return {a.x * s, a.y * s}; }
inline double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y; }
inline double cross(Vec a, Vec b) { return a.x * b.y - a.y * b.x; }
inline double norm(Vec a) { return std::sqrt(dot(a, a)); }

// {x : normals[k] . x <= offsets[k] for every k}, the normals unit and in counter-clockwise order.
struct ConvexPolygon {
  std::vector<Vec> normals;
  std::vector<double> offsets;
};

// The cone from a point counter-clockwise from `right` to `left`, or the whole plane.
struct Wedge {
  bool full = true;
  Vec right;
  Vec left;
};

// A convex polygon's corners and, seen from p, the corners where the sides facing p stop (`right`)
// and start (`left`). They depend on p only through which sides face it. No corners when it is
// empty or has no area.
struct Base {
  Polygon corners;
  bool contains = false;
  Vec right;
  Vec left;
};

[[nodiscard]] ConvexPolygon regular_polygon(std::size_t m, double circumradius);
[[nodiscard]] ConvexPolygon translated(const ConvexPolygon& poly, Vec v);
// The half-planes of a counter-clockwise convex polygon's edges, given its corners.
[[nodiscard]] ConvexPolygon from_corners(const Polygon& corners);

[[nodiscard]] Vec meet(Vec n1, double b1, Vec n2, double b2);
// The sides that are edges, counter-clockwise; none when the polygon is empty or unbounded.
[[nodiscard]] std::vector<std::size_t> edges(const ConvexPolygon& poly);
[[nodiscard]] Base base(const ConvexPolygon& poly, Vec p);
// The cone from p that the polygon spans: the whole plane when p is in it.
[[nodiscard]] Wedge cone(const Base& b, Vec p);

[[nodiscard]] bool within(const Wedge& w, Vec d);
// A cone of zero width is a tangency: the witness of the size bound always leaves area.
[[nodiscard]] std::optional<Wedge> intersect(const Wedge& a, const Wedge& b);

[[nodiscard]] Polygon clip_left_of(const Polygon& poly, Vec p, Vec d);
[[nodiscard]] Polygon cut(Polygon poly, Vec p, const Wedge& w);
// Only clearly inside counts as inside, as the point found becomes an output vertex.
[[nodiscard]] Vec closest(const Polygon& poly, Vec q);

// Upper and lower hull together, counter-clockwise from the leftmost point, corners only. Andrew's
// monotone chain, the form of Graham's scan that sorts by x then y: O(n) when `sorted`, which
// the caller vouches for, O(n log n) otherwise. Throws if what it built is not convex.
[[nodiscard]] Polygon full_convex_hull(Polygon points, bool sorted);
// Every corner a strict left turn, once around: a star turns left at every corner too.
[[nodiscard]] bool is_convex(const Polygon& poly);

}  // namespace ssk::geometry
