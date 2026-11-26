// The plane geometry the greedy stabbers stand on, against what geometry.hpp promises. Where a
// routine can be checked independently it is: brute force over the half-planes, exact
// arithmetic on lattice points, angles for the cones, and the obtuse-angle property of a
// projection.
#include "geometry/geometry.hpp"
#include "test_main.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace geo = ssk::geometry;
using geo::Base;
using geo::ConvexPolygon;
using geo::Polygon;
using geo::Vec;
using geo::Wedge;

constexpr double kPi = std::numbers::pi;
constexpr double kDegree = kPi / 180.0;

std::string str(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%g", x);
  return buf;
}

std::string str(Vec v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "(%.17g, %.17g)", v.x, v.y);
  return buf;
}

Vec polar(double angle, double r) { return {r * std::cos(angle), r * std::sin(angle)}; }

double distance(Vec a, Vec b) { return geo::norm(a - b); }

double area(const Polygon& poly) {
  double twice = 0.0;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    twice += geo::cross(poly[i], poly[(i + 1) % poly.size()]);
  }
  return twice / 2.0;
}

// How far x is beyond the side it is furthest beyond; at most 0 in the polygon.
double beyond(const ConvexPolygon& poly, Vec x) {
  double worst = -std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < poly.normals.size(); ++k) {
    worst = std::max(worst, geo::dot(poly.normals[k], x) - poly.offsets[k]);
  }
  return worst;
}

// [x0, x1] x [y0, y1], its normals at 0, 90, 180 and 270 degrees.
ConvexPolygon box(double x0, double y0, double x1, double y1) {
  return {{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}, {x1, y1, -x0, -y0}};
}

// n corners on a circle, counter-clockwise, consecutive ones 0.2 to 1.8 of an even spacing apart.
// Points on a circle in angular order are a convex polygon, so this needs no hull.
Polygon on_circle(std::mt19937& rng, std::size_t n, Vec center, double r) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  const double spacing = 2.0 * kPi / static_cast<double>(n);
  const double turn = 2.0 * kPi * unit(rng);
  Polygon out;
  for (std::size_t k = 0; k < n; ++k) {
    out.push_back(center + polar(turn + spacing * (static_cast<double>(k) + 0.8 * unit(rng)), r));
  }
  return out;
}

Wedge arc(double from, double to) { return {false, polar(from, 1.0), polar(to, 1.0)}; }

bool same(const Wedge& a, const Wedge& b) {
  return a.full == b.full && a.right == b.right && a.left == b.left;
}

bool lexicographic(Vec a, Vec b) { return a.x < b.x || (a.x == b.x && a.y < b.y); }

TEST("geometry/vectors: sum, difference, scaling, dot, cross and norm") {
  const Vec a{3, 4};
  const Vec b{-1, 2};
  CHECK((a + b == Vec{2, 6}));
  CHECK((a - b == Vec{4, 2}));
  CHECK((a * -0.5 == Vec{-1.5, -2}));
  CHECK(geo::dot(a, b) == 5.0);
  CHECK(geo::cross(a, b) == 10.0);
  CHECK(geo::cross(b, a) == -10.0);
  CHECK(geo::cross(a, a * 2.0) == 0.0);
  CHECK(geo::norm(a) == 5.0);
  CHECK(geo::cross(Vec{1, 0}, Vec{0, 1}) == 1.0);
}

TEST("geometry/a regular polygon: normals from (1, 0), counter-clockwise, corners on a circle") {
  for (const std::size_t m : {3u, 4u, 6u, 7u, 32u}) {
    const double r = 2.5;
    const double turn = 2.0 * kPi / static_cast<double>(m);
    const ConvexPolygon p = geo::regular_polygon(m, r);
    CHECK(p.normals.size() == m && p.offsets.size() == m);
    CHECK((p.normals[0] == Vec{1, 0}));
    for (std::size_t k = 0; k < m; ++k) {
      CHECK_CLOSE(geo::norm(p.normals[k]), 1.0, 1e-15);
      CHECK_CLOSE(p.offsets[k], r * std::cos(kPi / static_cast<double>(m)), 1e-15);
      CHECK_CLOSE(geo::cross(p.normals[k], p.normals[(k + 1) % m]), std::sin(turn), 1e-15);
      CHECK_CLOSE(geo::dot(p.normals[k], p.normals[(k + 1) % m]), std::cos(turn), 1e-15);
    }
    const Polygon corners = geo::base(p, {0, 0}).corners;
    CHECK(corners.size() == m);
    for (std::size_t k = 0; k < corners.size(); ++k) {
      const Vec expected = polar(kPi * static_cast<double>(2 * k + 1) / static_cast<double>(m), r);
      CHECK_MSG(distance(corners[k], expected) <= 1e-14 * r,
                "m " + std::to_string(m) + " corner " + std::to_string(k) + " " +
                    str(corners[k]) + " against " + str(expected));
    }
  }
}

// What the m-gon region is for: it holds the ball Lemma 2's witness needs, and reaches the error
// bound, at its corners, without passing it.
TEST("geometry/a regular polygon holds the disk of its apothem, and lies in its circumcircle") {
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> angle(0.0, 2.0 * kPi);
  for (const std::size_t m : {3u, 5u, 6u, 13u, 40u}) {
    const double r = 0.75;
    const ConvexPolygon p = geo::regular_polygon(m, r);
    const double apothem = r * std::cos(kPi / static_cast<double>(m));
    for (int i = 0; i < 2000; ++i) {
      const double a = angle(rng);
      CHECK(beyond(p, polar(a, apothem * (1.0 - 1e-12))) <= 0.0);
      CHECK(beyond(p, polar(a, r * (1.0 + 1e-12))) > 0.0);
    }
    for (std::size_t k = 0; k < m; ++k) {
      const double corner = kPi * static_cast<double>(2 * k + 1) / static_cast<double>(m);
      CHECK(beyond(p, polar(corner, r * (1.0 - 1e-12))) <= 0.0);
    }
  }
}

TEST("geometry/translating by v keeps the normals and moves every corner by v") {
  const ConvexPolygon p = geo::regular_polygon(7, 1.3);
  const Polygon at_origin = geo::base(p, {0, 0}).corners;
  for (const Vec v : {Vec{0, 0}, Vec{3.5, -7.25}, Vec{-1e6, 2e5}}) {
    const ConvexPolygon t = geo::translated(p, v);
    CHECK(t.normals == p.normals);
    const Polygon moved = geo::base(t, v).corners;
    CHECK(moved.size() == at_origin.size());
    for (std::size_t k = 0; k < moved.size(); ++k) {
      CHECK_MSG(distance(moved[k], at_origin[k] + v) <= 1e-14 * (geo::norm(v) + 1.3),
                "by " + str(v) + ": " + str(moved[k]) + " against " + str(at_origin[k] + v));
    }
  }
  CHECK(geo::translated(p, {0, 0}).offsets == p.offsets);
  const Vec u{0.25, -4};
  const Vec w{-3, 1.5};
  const ConvexPolygon twice = geo::translated(geo::translated(p, u), w);
  const ConvexPolygon once = geo::translated(p, u + w);
  for (std::size_t k = 0; k < p.offsets.size(); ++k) {
    CHECK_CLOSE(twice.offsets[k], once.offsets[k], 1e-14);
  }
}

TEST("geometry/from_corners makes each counter-clockwise edge an outward unit half-plane") {
  const ConvexPolygon square = geo::from_corners({{0, 0}, {2, 0}, {2, 2}, {0, 2}});
  CHECK((square.normals == std::vector<Vec>{{0, -1}, {1, 0}, {0, 1}, {-1, 0}}));
  CHECK((square.offsets == std::vector<double>{0, 2, 2, 0}));

  std::mt19937 rng(3);
  std::uniform_real_distribution<double> place(-500.0, 500.0);
  for (int trial = 0; trial < 300; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 10);
    const Vec center{place(rng), place(rng)};
    const double r = 0.1 + trial % 7;
    const Polygon corners = on_circle(rng, n, center, r);
    const ConvexPolygon p = geo::from_corners(corners);
    const double tol = 1e-13 * (geo::norm(center) + r);
    CHECK(p.normals.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
      CHECK_CLOSE(geo::norm(p.normals[i]), 1.0, 1e-15);
      CHECK(std::abs(geo::dot(p.normals[i], corners[i]) - p.offsets[i]) <= tol);
      CHECK(std::abs(geo::dot(p.normals[i], corners[(i + 1) % n]) - p.offsets[i]) <= tol);
      CHECK(beyond(p, corners[i]) <= tol);
    }
    if (n >= 4) {
      CHECK(beyond(p, center) < -0.1 * r);
    }
  }
}

TEST("geometry/from_corners then base gives the corners back") {
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> place(-3000.0, 3000.0);
  for (int trial = 0; trial < 300; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 12);
    const Vec center{place(rng), place(rng)};
    const double r = 1.0 + trial % 5;
    const Polygon corners = on_circle(rng, n, center, r);
    Vec mean{0, 0};
    for (const Vec c : corners) {
      mean = mean + c * (1.0 / static_cast<double>(n));
    }
    const Polygon back = geo::base(geo::from_corners(corners), mean).corners;
    CHECK(back.size() == n);
    for (std::size_t i = 0; i < back.size(); ++i) {
      CHECK(distance(back[i], corners[(i + 1) % n]) <= 1e-12 * (geo::norm(center) + r));
    }
  }
}

TEST("geometry/meet is the point on both lines, the same whichever comes first") {
  CHECK((geo::meet({1, 0}, 3, {0, 1}, -2) == Vec{3, -2}));
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> angle(0.0, 2.0 * kPi);
  std::uniform_real_distribution<double> offset(-100.0, 100.0);
  std::uniform_real_distribution<double> apart(0.05, kPi - 0.05);
  for (int i = 0; i < 10000; ++i) {
    const double a = angle(rng);
    const Vec n1 = polar(a, 1.0);
    const Vec n2 = polar(a + apart(rng), 1.0);
    const double b1 = offset(rng);
    const double b2 = offset(rng);
    const Vec x = geo::meet(n1, b1, n2, b2);
    const double tol = 1e-12 * (geo::norm(x) + 100.0);
    CHECK(std::abs(geo::dot(n1, x) - b1) <= tol);
    CHECK(std::abs(geo::dot(n2, x) - b2) <= tol);
    CHECK((geo::meet(n2, b2, n1, b1) == x));
  }
}

TEST("geometry/edges of a polygon each of whose sides bounds it are all of them, in order") {
  for (std::size_t m = 3; m <= 40; ++m) {
    std::vector<std::size_t> all(m);
    for (std::size_t k = 0; k < m; ++k) {
      all[k] = k;
    }
    CHECK_MSG(geo::edges(geo::regular_polygon(m, 1.0)) == all, "m " + std::to_string(m));
  }
  // A half-plane tangent to the unit circle touches it where every other one holds, so each is
  // an edge, however unevenly they are spread, while no gap reaches half a turn.
  std::mt19937 rng(9);
  for (int trial = 0; trial < 500; ++trial) {
    const std::size_t n = 4 + static_cast<std::size_t>(trial % 30);
    ConvexPolygon p;
    for (const Vec c : on_circle(rng, n, {0, 0}, 1.0)) {
      p.normals.push_back(c);
      p.offsets.push_back(1.0);
    }
    std::vector<std::size_t> all(n);
    for (std::size_t k = 0; k < n; ++k) {
      all[k] = k;
    }
    CHECK_MSG(geo::edges(p) == all, "trial " + std::to_string(trial));
  }
}

TEST("geometry/edges leaves out a side that bounds nothing and keeps one that cuts a corner") {
  const Vec diagonal{std::sqrt(0.5), std::sqrt(0.5)};
  const ConvexPolygon far{{{1, 0}, diagonal, {0, 1}, {-1, 0}, {0, -1}}, {1, 10, 1, 1, 1}};
  CHECK((geo::edges(far) == std::vector<std::size_t>{0, 2, 3, 4}));
  CHECK_CLOSE(area(geo::base(far, {0, 0}).corners), 4.0, 1e-14);

  ConvexPolygon cutting = far;
  cutting.offsets[1] = 1.5 * std::sqrt(0.5);
  CHECK((geo::edges(cutting) == std::vector<std::size_t>{0, 1, 2, 3, 4}));
  CHECK_CLOSE(area(geo::base(cutting, {0, 0}).corners), 4.0 - 0.125, 1e-14);
}

TEST("geometry/edges are none when the half-planes leave the polygon unbounded or empty") {
  CHECK(geo::edges({}).empty());
  CHECK(geo::edges({{{1, 0}, {0, 1}}, {1, 1}}).empty());
  CHECK(geo::edges({{{1, 0}, {0, 1}, polar(135 * kDegree, 1.0)}, {1, 1, 1}}).empty());
  CHECK(geo::edges(box(1, -1, -1, 1)).empty());
  CHECK(geo::edges(box(-1, 3, 1, 2)).empty());
}

// A strip open at one end has two sides half a turn apart: exactly, or as cos and sin leave it,
// a cross product of +1.2e-16 rather than 0.
TEST("geometry/edges are none for a strip open at one end, however its half turn rounds") {
  CHECK(geo::edges({{{1, 0}, {0, 1}, {0, -1}}, {1, 1, 1}}).empty());
  const ConvexPolygon rounded{{polar(0, 1.0), polar(90 * kDegree, 1.0), polar(270 * kDegree, 1.0)},
                              {1, 1, 1}};
  CHECK(geo::cross(rounded.normals[1], rounded.normals[2]) > 0.0);
  CHECK(geo::edges(rounded).empty());
}

// Side k bounds the polygon when the part of its line inside every other side is longer than
// rounding: that needs no order and no deque. Negative when none of the line is inside.
double chord(const ConvexPolygon& poly, std::size_t k) {
  const Vec n = poly.normals[k];
  const Vec a = n * poly.offsets[k];
  const Vec d{-n.y, n.x};
  double lo = -std::numeric_limits<double>::infinity();
  double hi = std::numeric_limits<double>::infinity();
  for (std::size_t j = 0; j < poly.normals.size(); ++j) {
    if (j == k) {
      continue;
    }
    const double toward = geo::dot(poly.normals[j], d);
    const double room = poly.offsets[j] - geo::dot(poly.normals[j], a);
    if (toward > 0.0) {
      hi = std::min(hi, room / toward);
    } else if (toward < 0.0) {
      lo = std::max(lo, room / toward);
    } else if (room < 0.0) {
      return -1.0;
    }
  }
  return hi - lo;
}

TEST("geometry/edges agrees with brute force on random half-planes, empty ones too") {
  std::mt19937 rng(13);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int bounded = 0;
  int empty = 0;
  int unbounded = 0;
  for (int trial = 0; trial < 20000; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 14);
    std::vector<double> angles(n);
    for (double& a : angles) {
      a = 2.0 * kPi * unit(rng);
    }
    std::sort(angles.begin(), angles.end());
    double gap = angles.front() + 2.0 * kPi - angles.back();
    for (std::size_t k = 1; k < n; ++k) {
      gap = std::max(gap, angles[k] - angles[k - 1]);
    }
    ConvexPolygon p;
    for (const double a : angles) {
      p.normals.push_back(polar(a, 1.0));
      p.offsets.push_back(-0.5 + 2.0 * unit(rng));
    }
    const std::vector<std::size_t> got = geo::edges(p);
    if (gap > kPi + 1e-6) {
      CHECK_MSG(got.empty(), "trial " + std::to_string(trial) + ": unbounded");
      ++unbounded;
      continue;
    }
    if (gap > kPi - 1e-6) {
      continue;
    }
    std::vector<std::size_t> expected;
    bool clear = true;
    for (std::size_t k = 0; k < n; ++k) {
      const double length = chord(p, k);
      clear = clear && std::abs(length) > 1e-9;
      if (length > 0.0) {
        expected.push_back(k);
      }
    }
    if (!clear) {
      continue;
    }
    CHECK_MSG(got == expected, "trial " + std::to_string(trial) + ": " +
                                   std::to_string(got.size()) + " edges, brute force " +
                                   std::to_string(expected.size()));
    ++(expected.empty() ? empty : bounded);
  }
  CHECK_MSG(bounded > 3000 && empty > 1000 && unbounded > 1000,
            std::to_string(bounded) + " bounded, " + std::to_string(empty) + " empty, " +
                std::to_string(unbounded) + " unbounded");
}

TEST("geometry/base from inside: every corner, counter-clockwise, and contains") {
  const Base b = geo::base(box(-1, -1, 1, 1), {0.3, -0.2});
  CHECK(b.contains);
  CHECK((b.corners == Polygon{{1, 1}, {-1, 1}, {-1, -1}, {1, -1}}));
}

// Seen from p the polygon spans the cone from `right` to `left`, counter-clockwise. From (3, 0),
// looking along -x, the right-hand side is +y.
TEST("geometry/base from outside: the corners where the sides facing p stop and start") {
  struct Case {
    Vec p;
    Vec right;
    Vec left;
  };
  for (const auto& [p, right, left] :
       {Case{{3, 0}, {1, 1}, {1, -1}}, Case{{3, 3}, {-1, 1}, {1, -1}},
        Case{{0, -5}, {1, -1}, {-1, -1}}, Case{{-2, 4}, {-1, -1}, {1, 1}}}) {
    const Base b = geo::base(box(-1, -1, 1, 1), p);
    CHECK_MSG(!b.contains, "from " + str(p));
    CHECK_MSG(b.corners.size() == 4, "from " + str(p));
    CHECK_MSG(b.right == right, "from " + str(p) + ": right " + str(b.right));
    CHECK_MSG(b.left == left, "from " + str(p) + ": left " + str(b.left));
  }
}

TEST("geometry/base takes a point on a side or at a corner as in the polygon") {
  for (const Vec p : {Vec{1, 0}, Vec{1, 1}, Vec{-1, -1}, Vec{0, -1}, Vec{-1, 0.5}}) {
    CHECK_MSG(geo::base(box(-1, -1, 1, 1), p).contains, str(p));
  }
}

// What the shared stabber relies on to build Q[p] once for all the starts between two marks.
TEST("geometry/base depends on p only through which sides face it") {
  const Base a = geo::base(box(-1, -1, 1, 1), {3, 0});
  const Base b = geo::base(box(-1, -1, 1, 1), {7.5, 0.9});
  CHECK(a.corners == b.corners && a.right == b.right && a.left == b.left);

  const ConvexPolygon p = geo::translated(geo::regular_polygon(9, 1.0), {4, -2});
  std::mt19937 rng(19);
  std::uniform_real_distribution<double> place(-10.0, 10.0);
  std::map<std::vector<bool>, Base> first;
  int repeats = 0;
  for (int i = 0; i < 5000; ++i) {
    const Vec q{4.0 + place(rng), -2.0 + place(rng)};
    std::vector<bool> facing;
    for (std::size_t k = 0; k < p.normals.size(); ++k) {
      facing.push_back(geo::dot(p.normals[k], q) > p.offsets[k]);
    }
    const Base got = geo::base(p, q);
    const auto [it, fresh] = first.try_emplace(facing, got);
    if (!fresh) {
      ++repeats;
      CHECK_MSG(it->second.corners == got.corners && it->second.right == got.right &&
                    it->second.left == got.left && it->second.contains == got.contains,
                "from " + str(q));
    }
  }
  CHECK(repeats > 4000);
}

TEST("geometry/from outside, the cone through base's two corners holds the whole polygon") {
  std::mt19937 rng(23);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int trial = 0; trial < 3000; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 12);
    const Vec center{200.0 * unit(rng) - 100.0, 200.0 * unit(rng) - 100.0};
    const double r = 0.5 + 5.0 * unit(rng);
    const Polygon corners = on_circle(rng, n, center, r);
    const Vec p = center + polar(2.0 * kPi * unit(rng), r * (1.05 + 20.0 * unit(rng)));
    const Base b = geo::base(geo::from_corners(corners), p);
    CHECK(!b.contains);
    CHECK(b.corners.size() == n);
    CHECK(std::find(b.corners.begin(), b.corners.end(), b.right) != b.corners.end());
    CHECK(std::find(b.corners.begin(), b.corners.end(), b.left) != b.corners.end());
    const Vec right = b.right - p;
    const Vec left = b.left - p;
    CHECK(geo::cross(right, left) > 0.0);
    for (const Vec c : b.corners) {
      const double tol = 1e-12 * geo::norm(c - p) * (geo::norm(p) + r);
      CHECK_MSG(geo::cross(right, c - p) >= -tol && geo::cross(c - p, left) >= -tol,
                "trial " + std::to_string(trial) + ": corner " + str(c) + " outside the cone");
    }
  }
}

TEST("geometry/base of an empty polygon, or of one without area, has no corners") {
  CHECK(geo::base(box(1, -1, -1, 1), {0, 0}).corners.empty());
  CHECK(geo::base(box(0, 0, 0, 0), {0, 0}).corners.empty());
  CHECK(geo::base(box(-1, 0, 1, 0), {0, 0}).corners.empty());
  CHECK(geo::base(box(-1, 0, 1, 0), {0, 5}).corners.empty());
}

TEST("geometry/translating the polygon and the point together moves base and nothing else") {
  const ConvexPolygon p = geo::regular_polygon(11, 2.0);
  for (const Vec q : {Vec{0.5, 0.3}, Vec{5, 1}, Vec{-3, -4}}) {
    const Base a = geo::base(p, q);
    for (const Vec v : {Vec{10, -20}, Vec{-1e4, 3e3}}) {
      const Base b = geo::base(geo::translated(p, v), q + v);
      const double tol = 1e-14 * (geo::norm(v) + 2.0) * 4.0;
      CHECK(a.contains == b.contains);
      CHECK(a.corners.size() == b.corners.size());
      for (std::size_t k = 0; k < a.corners.size(); ++k) {
        CHECK(distance(b.corners[k], a.corners[k] + v) <= tol);
      }
      if (!a.contains) {
        CHECK(distance(b.right, a.right + v) <= tol);
        CHECK(distance(b.left, a.left + v) <= tol);
      }
    }
  }
}

TEST("geometry/the cone from inside is the whole plane, from outside the rays to base's corners") {
  const ConvexPolygon square = box(-1, -1, 1, 1);
  const Vec inside{0.5, 0.5};
  CHECK(geo::cone(geo::base(square, inside), inside).full);
  const Vec p{3, 0};
  const Wedge w = geo::cone(geo::base(square, p), p);
  CHECK(!w.full);
  CHECK((w.right == Vec{-2, 1}));
  CHECK((w.left == Vec{-2, -1}));
}

TEST("geometry/within a cone is between its rays, the rays included, whatever their length") {
  for (const Wedge& quarter : {Wedge{false, {1, 0}, {0, 1}}, Wedge{false, {7, 0}, {0, 0.25}}}) {
    for (const Vec d : {Vec{1, 1}, Vec{1, 0}, Vec{0, 1}, Vec{5, 0.1}, Vec{1e-9, 1}, Vec{3, 3}}) {
      CHECK_MSG(geo::within(quarter, d), str(d) + " is in");
    }
    for (const Vec d : {Vec{-1, 0}, Vec{0, -1}, Vec{-1, -1}, Vec{1, -1e-9}, Vec{-1e-9, 1},
                        Vec{-1, 1}, Vec{1, -1}}) {
      CHECK_MSG(!geo::within(quarter, d), str(d) + " is out");
    }
  }
  for (const Vec d : {Vec{1, 0}, Vec{-1, -1}, Vec{0, -3}, Vec{-2, 0.5}}) {
    CHECK_MSG(geo::within(Wedge{}, d), str(d) + " is in the plane");
  }
}

TEST("geometry/intersect: what two cones share, and none when that is a ray or nothing") {
  const Wedge plane{};
  const Wedge a = arc(0, 60 * kDegree);
  const Wedge b = arc(30 * kDegree, 90 * kDegree);
  const Wedge inner = arc(10 * kDegree, 20 * kDegree);

  const std::optional<Wedge> ab = geo::intersect(a, b);
  CHECK(ab && same(*ab, {false, b.right, a.left}));
  const std::optional<Wedge> ba = geo::intersect(b, a);
  CHECK(ba && same(*ba, {false, b.right, a.left}));
  CHECK(geo::intersect(a, inner) && same(*geo::intersect(a, inner), inner));
  CHECK(geo::intersect(inner, a) && same(*geo::intersect(inner, a), inner));
  CHECK(geo::intersect(a, a) && same(*geo::intersect(a, a), a));

  CHECK(geo::intersect(plane, plane) && geo::intersect(plane, plane)->full);
  CHECK(geo::intersect(plane, a) && same(*geo::intersect(plane, a), a));
  CHECK(geo::intersect(a, plane) && same(*geo::intersect(a, plane), a));

  CHECK(!geo::intersect(arc(0, 30 * kDegree), arc(60 * kDegree, 90 * kDegree)));
  CHECK(!geo::intersect(arc(60 * kDegree, 90 * kDegree), arc(0, 30 * kDegree)));
  CHECK(!geo::intersect(arc(0, 60 * kDegree), arc(180 * kDegree, 240 * kDegree)));

  const Vec shared = polar(30 * kDegree, 1.0);
  CHECK(!geo::intersect({false, {1, 0}, shared}, {false, shared, polar(60 * kDegree, 1.0)}));
  CHECK(!geo::intersect({false, shared, polar(60 * kDegree, 1.0)}, {false, {1, 0}, shared}));
}

// Cones narrower than half a turn share one arc or none, so the answer is also an interval of
// angles, worked out without cross products.
TEST("geometry/intersect agrees with the arc the two cones share, worked out in angles") {
  std::mt19937 rng(29);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int shared = 0;
  int apart = 0;
  for (int trial = 0; trial < 20000; ++trial) {
    const double s1 = 2.0 * kPi * unit(rng);
    const double w1 = 0.01 + (kPi - 0.02) * unit(rng);
    const double s2 = 2.0 * kPi * unit(rng);
    const double w2 = 0.01 + (kPi - 0.02) * unit(rng);
    const Wedge a{false, polar(s1, 0.1 + 10.0 * unit(rng)), polar(s1 + w1, 0.1 + 10.0 * unit(rng))};
    const Wedge b{false, polar(s2, 0.1 + 10.0 * unit(rng)), polar(s2 + w2, 0.1 + 10.0 * unit(rng))};

    const double ahead = std::fmod(s2 - s1 + 4.0 * kPi, 2.0 * kPi);
    std::optional<std::pair<double, double>> arc_of;
    if (ahead < w1) {
      arc_of = std::pair{s2, std::min(w1 - ahead, w2)};
    } else if (ahead > 2.0 * kPi - w2) {
      arc_of = std::pair{s1, std::min(w2 - (2.0 * kPi - ahead), w1)};
    }
    if (std::abs(ahead - w1) < 1e-9 || std::abs(ahead - (2.0 * kPi - w2)) < 1e-9 ||
        (arc_of && arc_of->second < 1e-9)) {
      continue;
    }
    for (const auto& got : {geo::intersect(a, b), geo::intersect(b, a)}) {
      if (!arc_of) {
        CHECK_MSG(!got, "trial " + std::to_string(trial) + ": the cones share nothing");
        continue;
      }
      CHECK_MSG(got && !got->full, "trial " + std::to_string(trial) + ": the cones share an arc");
      if (got && !got->full) {
        const Vec right = got->right * (1.0 / geo::norm(got->right));
        const Vec left = got->left * (1.0 / geo::norm(got->left));
        CHECK(distance(right, polar(arc_of->first, 1.0)) <= 1e-12);
        CHECK(distance(left, polar(arc_of->first + arc_of->second, 1.0)) <= 1e-12);
      }
    }
    ++(arc_of ? shared : apart);
  }
  CHECK(shared > 5000 && apart > 5000);
}

TEST("geometry/clip_left_of keeps what is left of the line through p along d, the line included") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  CHECK((geo::clip_left_of(square, {1, 0}, {0, 1}) == Polygon{{0, 0}, {1, 0}, {1, 2}, {0, 2}}));
  CHECK((geo::clip_left_of(square, {1, 0}, {0, -1}) == Polygon{{1, 0}, {2, 0}, {2, 2}, {1, 2}}));
  CHECK(geo::clip_left_of(square, {5, 0}, {0, 1}) == square);
  CHECK(geo::clip_left_of(square, {-5, 0}, {0, 1}).empty());
  CHECK(geo::clip_left_of({}, {0, 0}, {1, 0}).empty());
}

TEST("geometry/clip_left_of along an edge keeps the polygon or that edge, no corner twice") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  CHECK(geo::clip_left_of(square, {0, 0}, {1, 0}) == square);
  CHECK((geo::clip_left_of(square, {0, 0}, {-1, 0}) == Polygon{{0, 0}, {2, 0}}));
  CHECK(geo::clip_left_of(square, {2, 2}, {0, 1}) == square);
}

TEST("geometry/the two sides of a line split a polygon's area between them") {
  std::mt19937 rng(31);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int trial = 0; trial < 3000; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 12);
    const Vec center{200.0 * unit(rng) - 100.0, 200.0 * unit(rng) - 100.0};
    const double r = 0.5 + 10.0 * unit(rng);
    const Polygon poly = on_circle(rng, n, center, r);
    const Vec p = center + polar(2.0 * kPi * unit(rng), 1.2 * r * unit(rng));
    const Vec d = polar(2.0 * kPi * unit(rng), 0.1 + 10.0 * unit(rng));
    const Polygon left = geo::clip_left_of(poly, p, d);
    const Polygon right = geo::clip_left_of(poly, p, d * -1.0);
    const double scale = geo::norm(center) + r;
    const double tol = 1e-13 * scale * scale * static_cast<double>(n);
    CHECK_MSG(std::abs(area(left) + area(right) - area(poly)) <= tol,
              "trial " + std::to_string(trial));
    CHECK(area(left) >= -tol && area(right) >= -tol);
    for (const Vec q : left) {
      CHECK(geo::cross(d, q - p) >= -1e-12 * geo::norm(d) * scale);
    }
    for (const Vec q : right) {
      CHECK(geo::cross(d, q - p) <= 1e-12 * geo::norm(d) * scale);
    }
  }
}

TEST("geometry/cut by the whole plane is the polygon as it is") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  CHECK(geo::cut(square, {7, 7}, Wedge{}) == square);
}

TEST("geometry/cut keeps the part of the polygon in the cone from p, nothing if it misses") {
  const Polygon rect{{0, -1}, {2, -1}, {2, 1}, {0, 1}};
  const Polygon got = geo::cut(rect, {0, 0}, Wedge{false, {1, -1}, {1, 1}});
  CHECK((got == Polygon{{1, -1}, {2, -1}, {2, 1}, {1, 1}, {0, 0}}));
  CHECK(area(got) == 3.0);
  const Polygon square{{1, 1}, {2, 1}, {2, 2}, {1, 2}};
  CHECK(geo::cut(square, {0, 0}, Wedge{false, {-1, 0.5}, {-1, -0.5}}).empty());
}

TEST("geometry/cut by the polygon's own cone keeps all of it, and its two halves add up to it") {
  std::mt19937 rng(37);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int trial = 0; trial < 2000; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 12);
    const Vec center{200.0 * unit(rng) - 100.0, 200.0 * unit(rng) - 100.0};
    const double r = 0.5 + 5.0 * unit(rng);
    const Polygon poly = on_circle(rng, n, center, r);
    const Vec p = center + polar(2.0 * kPi * unit(rng), r * (1.5 + 10.0 * unit(rng)));
    const Wedge w = geo::cone(geo::base(geo::from_corners(poly), p), p);
    const double whole = area(poly);
    const double tol = 1e-12 * (geo::norm(center) + r) * (geo::norm(center) + r);
    CHECK_MSG(std::abs(area(geo::cut(poly, p, w)) - whole) <= tol,
              "trial " + std::to_string(trial));
    const Vec middle = w.right * (1.0 / geo::norm(w.right)) + w.left * (1.0 / geo::norm(w.left));
    const double one = area(geo::cut(poly, p, {false, w.right, middle}));
    const double other = area(geo::cut(poly, p, {false, middle, w.left}));
    CHECK(one > 0.0 && other > 0.0);
    CHECK_MSG(std::abs(one + other - whole) <= tol, "trial " + std::to_string(trial));
  }
}

TEST("geometry/closest returns a point clearly inside as it is") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  for (const Vec q : {Vec{1, 1}, Vec{0.5, 1.9}, Vec{1.999, 0.001}, Vec{1, 1e-6}}) {
    CHECK_MSG(geo::closest(square, q) == q, str(q));
  }
  std::mt19937 rng(41);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int trial = 0; trial < 2000; ++trial) {
    const Vec center{200.0 * unit(rng) - 100.0, 200.0 * unit(rng) - 100.0};
    const double r = 0.5 + 5.0 * unit(rng);
    const Polygon poly = on_circle(rng, 6 + static_cast<std::size_t>(trial % 10), center, r);
    const Vec q = center + polar(2.0 * kPi * unit(rng), 0.5 * r * unit(rng));
    CHECK_MSG(geo::closest(poly, q) == q, "trial " + std::to_string(trial));
  }
}

TEST("geometry/closest from outside is the nearest point of the boundary") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  struct Case {
    Vec q;
    Vec nearest;
  };
  for (const auto& [q, nearest] : {Case{{3, 1}, {2, 1}}, Case{{3, 3}, {2, 2}},
                                   Case{{-1, -2}, {0, 0}}, Case{{1, -5}, {1, 0}},
                                   Case{{1, 2.5}, {1, 2}}, Case{{-0.5, 1.25}, {0, 1.25}}}) {
    const Vec got = geo::closest(square, q);
    CHECK_MSG(got == nearest, "from " + str(q) + ": " + str(got));
  }
}

// r is the point of a convex set nearest q exactly when r is in the set and q - r makes an
// obtuse angle with x - r for every x in it; for a polygon its corners are enough.
TEST("geometry/closest from outside is the projection onto the polygon") {
  std::mt19937 rng(43);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int trial = 0; trial < 3000; ++trial) {
    const std::size_t n = 3 + static_cast<std::size_t>(trial % 14);
    const Vec center{200.0 * unit(rng) - 100.0, 200.0 * unit(rng) - 100.0};
    const double r = 0.5 + 5.0 * unit(rng);
    const Polygon poly = on_circle(rng, n, center, r);
    const Vec q = center + polar(2.0 * kPi * unit(rng), r * (1.01 + 4.0 * unit(rng)));
    const Vec x = geo::closest(poly, q);
    const double scale = geo::norm(center) + r;
    CHECK_MSG(beyond(geo::from_corners(poly), x) <= 1e-12 * scale,
              "trial " + std::to_string(trial) + ": " + str(x) + " is not in the polygon");
    for (const Vec c : poly) {
      CHECK_MSG(geo::dot(q - x, c - x) <= 1e-12 * scale * scale,
                "trial " + std::to_string(trial) + ": corner " + str(c) + " is nearer");
      CHECK(distance(q, x) <= distance(q, c) + 1e-12 * scale);
    }
  }
}

// Only clearly inside counts as inside, as the point found becomes an output vertex: a query on
// a side, or inside it by less than a relative 1e-9, is put on it.
TEST("geometry/closest puts a point on the boundary, or just inside it, on the boundary") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  CHECK((geo::closest(square, {1, 1e-12}) == Vec{1, 0}));
  CHECK((geo::closest(square, {1, 0}) == Vec{1, 0}));
  CHECK((geo::closest(square, {2, 2}) == Vec{2, 2}));
  CHECK((geo::closest(square, {2 - 1e-12, 1}) == Vec{2, 1}));
}

// The shadow a few ulps wide where the regions around A and B of a spike meet corner to corner:
// skipping every edge shorter than 1e-9 of its distance to the query took a query a delta away
// as inside it (UPDATES, 2026-11-24).
TEST("geometry/closest to a tiny polygon a unit away is a point of it, not the query") {
  for (const double size : {1e-8, 1e-10, 1e-12, 1e-13}) {
    const Polygon tiny{{1, 0}, {1 + size, 0}, {1, size}};
    const Vec got = geo::closest(tiny, {0, 0});
    CHECK_MSG(distance(got, {1, 0}) <= 2.0 * size, "size " + str(size) + ": " + str(got));
  }
}

// Far from the origin the same polygon has every edge a few ulps of its coordinates long, so
// every edge is skipped, and q, however far, was taken as inside it.
TEST("geometry/closest to a polygon a few ulps wide is a point of it, however far it lies") {
  for (const double at : {1.0, 1e3, 1e6, 1e9}) {
    for (const double size : {1e-14, 2e-15}) {
      const Polygon tiny{{at, 0}, {at + size * at, 0}, {at, size * at}};
      for (const Vec q : {Vec{0, 0}, Vec{at - 1, 0}, Vec{at, -0.5 * at}}) {
        const Vec got = geo::closest(tiny, q);
        CHECK_MSG(distance(got, {at, 0}) <= 2.0 * size * at,
                  "at " + str(at) + ", size " + str(size) + ", from " + str(q) + ": " + str(got));
      }
    }
  }
  CHECK((geo::closest({{5, 5}}, {0, 0}) == Vec{5, 5}));
  CHECK((geo::closest({{5, 5}}, {5, 5}) == Vec{5, 5}));
}

TEST("geometry/the hull is the corners only, counter-clockwise from the lowest leftmost point") {
  const Polygon points{{1, 1}, {2, 2}, {0, 0}, {1, 0},     {2, 0}, {0, 2},
                       {2, 1}, {0, 1}, {1, 2}, {0.5, 1.5}, {2, 2}};
  CHECK((geo::full_convex_hull(points, false) == Polygon{{0, 0}, {2, 0}, {2, 2}, {0, 2}}));
}

TEST("geometry/the hull drops repeats, and fewer than three points come back sorted") {
  CHECK(geo::full_convex_hull({}, false).empty());
  CHECK((geo::full_convex_hull({{1, 1}, {1, 1}, {1, 1}}, false) == Polygon{{1, 1}}));
  CHECK((geo::full_convex_hull({{2, 0}, {1, 5}, {2, 0}}, false) == Polygon{{1, 5}, {2, 0}}));
  CHECK((geo::full_convex_hull({{0, 0}, {1, 0}, {0, 1}, {0, 0}, {1, 0}}, false) ==
         Polygon{{0, 0}, {1, 0}, {0, 1}}));
}

TEST("geometry/collinear points hull to their two ends") {
  CHECK((geo::full_convex_hull({{0, 0}, {3, 3}, {1, 1}, {2, 2}}, false) ==
         Polygon{{0, 0}, {3, 3}}));
  CHECK((geo::full_convex_hull({{1, 2}, {1, 0}, {1, 1}}, false) == Polygon{{1, 0}, {1, 2}}));
}

TEST("geometry/the hull of points already sorted is the same either way") {
  std::mt19937 rng(47);
  std::uniform_real_distribution<double> place(-50.0, 50.0);
  for (int trial = 0; trial < 500; ++trial) {
    Polygon points(3 + static_cast<std::size_t>(trial % 40));
    for (Vec& q : points) {
      q = {place(rng), place(rng)};
    }
    Polygon sorted = points;
    std::sort(sorted.begin(), sorted.end(), lexicographic);
    CHECK(geo::full_convex_hull(sorted, true) == geo::full_convex_hull(points, false));
  }
}

bool on_segment(Vec a, Vec b, Vec p) {
  return geo::cross(b - a, p - a) == 0.0 && std::min(a.x, b.x) <= p.x &&
         p.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
}

bool in_triangle(Vec a, Vec b, Vec c, Vec p) {
  const double turn = geo::cross(b - a, c - a);
  if (turn == 0.0) {
    return false;
  }
  const double d1 = geo::cross(b - a, p - a);
  const double d2 = geo::cross(c - b, p - b);
  const double d3 = geo::cross(a - c, p - c);
  return turn > 0.0 ? d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0 : d1 <= 0.0 && d2 <= 0.0 && d3 <= 0.0;
}

// p is in the hull of `others` when one of their triangles or segments holds it.
bool held(const Polygon& others, Vec p) {
  for (std::size_t i = 0; i < others.size(); ++i) {
    for (std::size_t j = i; j < others.size(); ++j) {
      if (on_segment(others[i], others[j], p)) {
        return true;
      }
      for (std::size_t k = j + 1; k < others.size(); ++k) {
        if (in_triangle(others[i], others[j], others[k], p)) {
          return true;
        }
      }
    }
  }
  return false;
}

// On small integers every cross product is exact, so points in line are exactly in line, and a
// point is a corner exactly when no triangle or segment of the other points holds it.
TEST("geometry/the hull of lattice points has exactly the corners brute force finds") {
  std::mt19937 rng(53);
  std::uniform_int_distribution<int> coordinate(-4, 4);
  for (int trial = 0; trial < 1500; ++trial) {
    Polygon points(1 + static_cast<std::size_t>(trial % 24));
    for (Vec& q : points) {
      q = {static_cast<double>(coordinate(rng)), static_cast<double>(coordinate(rng))};
    }
    Polygon distinct = points;
    std::sort(distinct.begin(), distinct.end(), lexicographic);
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    Polygon corners;
    for (const Vec p : distinct) {
      Polygon others;
      std::copy_if(distinct.begin(), distinct.end(), std::back_inserter(others),
                   [p](Vec q) { return !(q == p); });
      if (!held(others, p)) {
        corners.push_back(p);
      }
    }

    const Polygon hull = geo::full_convex_hull(points, false);
    Polygon got = hull;
    std::sort(got.begin(), got.end(), lexicographic);
    CHECK_MSG(got == corners, "trial " + std::to_string(trial) + ": " +
                                  std::to_string(hull.size()) + " corners, brute force " +
                                  std::to_string(corners.size()));
    CHECK(hull.empty() || hull[0] == distinct[0]);
    if (hull.size() >= 3) {
      CHECK(geo::is_convex(hull));
      for (const Vec q : points) {
        for (std::size_t i = 0; i < hull.size(); ++i) {
          CHECK(geo::cross(hull[(i + 1) % hull.size()] - hull[i], q - hull[i]) >= 0.0);
        }
      }
    }
  }
}

// The hull is checked once built, so a caller wrong about the order finds out.
TEST("geometry/the hull throws when points said to be sorted are not") {
  const Polygon unsorted{{-1, 5}, {2, 5}, {-5, -4}, {-2, 5}};
  bool threw = false;
  try {
    (void)geo::full_convex_hull(unsorted, true);
  } catch (const std::logic_error&) {
    threw = true;
  }
  CHECK(threw);
  CHECK((geo::full_convex_hull(unsorted, false) == Polygon{{-5, -4}, {2, 5}, {-2, 5}}));
}

TEST("geometry/is_convex: counter-clockwise with a strict left turn at every corner") {
  const Polygon square{{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  CHECK(geo::is_convex(square));
  CHECK(geo::is_convex({{2, 2}, {0, 2}, {0, 0}, {2, 0}}));
  CHECK(geo::is_convex({{0, 0}, {1, 0}, {0, 1}}));
  CHECK(!geo::is_convex({{0, 2}, {2, 2}, {2, 0}, {0, 0}}));
  CHECK(!geo::is_convex({}));
  CHECK(!geo::is_convex({{0, 0}}));
  CHECK(!geo::is_convex({{0, 0}, {1, 0}}));
  CHECK(!geo::is_convex({{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}}));
  CHECK(!geo::is_convex({{0, 0}, {2, 0}, {2, 0}, {2, 2}, {0, 2}}));
  CHECK(!geo::is_convex({{0, 0}, {2, 0}, {1, 1}, {2, 2}, {0, 2}}));
  Polygon many;
  for (int k = 0; k < 1000; ++k) {
    many.push_back(polar(2.0 * kPi * k / 1000.0, 1.0));
  }
  CHECK(geo::is_convex(many));
}

TEST("geometry/a star turns left at every corner and still is not convex") {
  Polygon star;
  for (const int k : {0, 2, 4, 1, 3}) {
    star.push_back(polar(0.1 + 2.0 * kPi * k / 5.0, 1.0));
  }
  for (std::size_t i = 0; i < star.size(); ++i) {
    const Vec a = star[i];
    const Vec b = star[(i + 1) % star.size()];
    const Vec c = star[(i + 2) % star.size()];
    CHECK(geo::cross(b - a, c - a) > 0.0);
  }
  CHECK(!geo::is_convex(star));
}

}  // namespace
