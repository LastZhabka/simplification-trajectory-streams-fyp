// The decision procedure is the primitive; the value is bisection on top. Cases here are
// ones whose Frechet distance can be worked out by hand, plus the properties the measure has
// to satisfy for a results table built on it to mean anything.
#include "frechet/frechet.hpp"
#include "test_main.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

using ssk::frechet::Curve;
using ssk::frechet::Frechet;
using ssk::frechet::Point;

Curve<2> line(std::size_t n, double dy = 0.0) {
  Curve<2> c;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({static_cast<double>(i), dy});
  }
  return c;
}

TEST("frechet/identical curves are at distance zero") {
  const Frechet<2> f;
  const Curve<2> c = line(10);
  CHECK(f.within(c, c, 0.0));
  CHECK_CLOSE(f.distance(c, c, 1e-9), 0.0, 1e-6);
}

TEST("frechet/parallel lines are their separation apart") {
  const Frechet<2> f;
  CHECK_CLOSE(f.distance(line(10), line(10, 1.0), 1e-9), 1.0, 1e-6);
  CHECK_CLOSE(f.distance(line(10), line(10, 0.25), 1e-9), 0.25, 1e-6);
}

TEST("frechet/a straight line matches its own endpoints exactly") {
  const Frechet<2> f;
  const Curve<2> ends{{0.0, 0.0}, {9.0, 0.0}};
  CHECK_CLOSE(f.distance(line(10), ends, 1e-9), 0.0, 1e-6);
}

TEST("frechet/a spike costs its own height against its chord") {
  const Frechet<2> f;
  const Curve<2> spike{{0.0, 0.0}, {1.0, 2.0}, {2.0, 0.0}};
  const Curve<2> chord{{0.0, 0.0}, {2.0, 0.0}};
  CHECK_CLOSE(f.distance(spike, chord, 1e-9), 2.0, 1e-6);
}

TEST("frechet/the matching must stay monotone") {
  const Frechet<2> f;
  const Curve<2> forward{{0.0, 0.0}, {1.0, 0.0}};
  const Curve<2> backward{{1.0, 0.0}, {0.0, 0.0}};
  CHECK_MSG(!f.within(forward, backward, 0.5),
            "a reversed curve cannot be matched by walking both forwards");
  CHECK_CLOSE(f.distance(forward, backward, 1e-9), 1.0, 1e-6);
}

TEST("frechet/decision and value agree either side of the answer") {
  const Frechet<2> f;
  const Curve<2> spike{{0.0, 0.0}, {1.0, 2.0}, {2.0, 0.0}};
  const Curve<2> chord{{0.0, 0.0}, {2.0, 0.0}};
  const double d = f.distance(spike, chord, 1e-9);
  CHECK_MSG(f.within(spike, chord, d * 1.001), "must be within a radius above the distance");
  CHECK_MSG(!f.within(spike, chord, d * 0.999), "must not be within a radius below it");
}

TEST("frechet/works in three dimensions") {
  const Frechet<3> f;
  Curve<3> a, b;
  for (std::size_t i = 0; i < 20; ++i) {
    const auto t = static_cast<double>(i);
    a.push_back({t, 0.0, 0.0});
    b.push_back({t, 0.0, 3.0});
  }
  CHECK_CLOSE(f.distance(a, b, 1e-9), 3.0, 1e-6);
}

TEST("frechet/degenerate curves do not break it") {
  const Frechet<2> f;
  const Curve<2> one{{0.0, 0.0}};
  const Curve<2> same_point{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  CHECK_CLOSE(f.distance(one, one, 1e-9), 0.0, 1e-6);
  CHECK_CLOSE(f.distance(one, same_point, 1e-9), 0.0, 1e-6);
  CHECK_CLOSE(f.distance(same_point, line(4), 1e-9), 3.0, 1e-6);

  const Curve<2> pair{{0.0, 0.0}, {0.0, 4.0}};
  CHECK_MSG(!f.within(one, pair, 3.0), "a single point must cover the whole other curve");
  CHECK(f.within(one, pair, 4.0));
}

TEST("frechet/a simplification is never further than its worst kept point") {
  // Dropping vertices can only move the curve by so much: the distance to a subsequence is
  // bounded by the largest gap it introduces. A sanity bound, not a tight one.
  const Frechet<2> f;
  std::uint32_t s = 999;
  Curve<2> c;
  double x = 0.0;
  double y = 0.0;
  for (std::size_t i = 0; i < 60; ++i) {
    c.push_back({x, y});
    s = s * 1664525u + 1013904223u;
    x += static_cast<double>((s >> 16) % 16u) / 8.0;
    s = s * 1664525u + 1013904223u;
    y += static_cast<double>((s >> 16) % 16u) / 8.0 - 1.0;
  }
  Curve<2> every_other;
  for (std::size_t i = 0; i < c.size(); i += 2) {
    every_other.push_back(c[i]);
  }
  if (every_other.back() != c.back()) {
    every_other.push_back(c.back());
  }

  double worst = 0.0;
  for (std::size_t i = 0; i + 1 < c.size(); ++i) {
    worst = std::max(worst, std::hypot(c[i + 1][0] - c[i][0], c[i + 1][1] - c[i][1]));
  }
  const double d = f.distance(c, every_other, 1e-9);
  CHECK_MSG(d <= worst + 1e-9, "dropping every other vertex cannot cost more than a step");
  CHECK_MSG(d > 0.0, "and it does cost something");
}

TEST("frechet/a point 1e-9 of a long segment off it is decided at its offset") {
  // The scale of an exact baseline beside a long NGSIM step, where b^2 - 4ac cancelled to
  // noise and about half of these came out wrong, both ways.
  const Frechet<2> f;
  std::uint32_t s = 2026;
  const auto next = [&s] {
    s = s * 1664525u + 1013904223u;
    return static_cast<double>(s >> 8) / 16777216.0;
  };
  for (int i = 0; i < 200; ++i) {
    const Point<2> p1{2.23e6 + 1e4 * next(), 1.37e6 + 1e4 * next()};
    const double angle = 6.283185307179586 * next();
    const Point<2> p2{p1[0] + 500.0 * std::cos(angle), p1[1] + 500.0 * std::sin(angle)};
    const double along = 0.05 + 0.9 * next();
    const double off = 5e-7 * (0.5 + next());
    const Point<2> q{p1[0] + along * (p2[0] - p1[0]) - off * std::sin(angle),
                     p1[1] + along * (p2[1] - p1[1]) + off * std::cos(angle)};

    const double ux = p2[0] - p1[0];
    const double uy = p2[1] - p1[1];
    const double d = std::fabs(ux * (q[1] - p1[1]) - uy * (q[0] - p1[0])) / std::hypot(ux, uy);
    const Curve<2> segment{p1, p2};
    const Curve<2> bent{p1, q, p2};
    CHECK_MSG(f.within(segment, bent, 1.01 * d), "within just above the offset");
    CHECK_MSG(!f.within(segment, bent, 0.99 * d), "and not just below it");
    CHECK_CLOSE(f.distance(segment, bent, 1e-12), d, 1e-3 * d);
  }
}

// Each is the smallest found on which one constraint of the diagram, taken out, says within when
// the distance is more: a point's whole curve, a matching that may not run back up a cell, the
// shortcut through a free top-right corner, and a bottom row that a far vertex cuts. The range
// is the discrete distance of the curves cut 400 times per edge, and the piece's length below.
TEST("frechet/curves on which a constraint of the free-space diagram decides") {
  struct Case {
    const char* what;
    Curve<2> a;
    Curve<2> b;
    double low;
    double high;
  };
  const Frechet<2> f;
  for (const Case& c : {
           Case{"the point's curve", {{4, 2}}, {{3, 3}, {0, 4}, {3, 3}}, 4.4642, 4.4721},
           Case{"monotone in a cell", {{3, 4}, {3, 2}, {0, 2}}, {{3, 4}, {1, 4}, {2, 4}, {0, 2}},
                2.0541, 2.0616},
           Case{"the free corner", {{1, 3}, {3, 2}, {2, 3}},
                {{1, 3}, {4, 1}, {1, 1}, {0, 3}, {2, 3}}, 2.1123, 2.1213},
           Case{"the bottom row", {{1, 2}, {1, 3}, {1, 2}, {2, 4}}, {{4, 2}, {4, 4}, {4, 2}},
                3.0358, 3.0414},
       }) {
    for (const auto& [a, b] : {std::pair{c.a, c.b}, std::pair{c.b, c.a}}) {
      CHECK_MSG(!f.within(a, b, c.low - 0.01), std::string(c.what) + ": within below it");
      CHECK_MSG(f.within(a, b, c.high + 0.01), std::string(c.what) + ": not within above it");
    }
  }
}

Curve<2> subdivide(const Curve<2>& c, int k) {
  Curve<2> out;
  for (std::size_t i = 0; i + 1 < c.size(); ++i) {
    for (int j = 0; j < k; ++j) {
      const double t = static_cast<double>(j) / k;
      out.push_back({c[i][0] + t * (c[i + 1][0] - c[i][0]), c[i][1] + t * (c[i + 1][1] - c[i][1])});
    }
  }
  out.push_back(c.back());
  return out;
}

double longest_edge(const Curve<2>& c) {
  double longest = 0.0;
  for (std::size_t i = 0; i + 1 < c.size(); ++i) {
    longest = std::max(longest, std::hypot(c[i + 1][0] - c[i][0], c[i + 1][1] - c[i][1]));
  }
  return longest;
}

double discrete_frechet(const Curve<2>& a, const Curve<2>& b) {
  std::vector<double> prev(b.size());
  std::vector<double> cur(b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    for (std::size_t j = 0; j < b.size(); ++j) {
      const double d = std::hypot(a[i][0] - b[j][0], a[i][1] - b[j][1]);
      if (i == 0) {
        cur[j] = j == 0 ? d : std::max(cur[j - 1], d);
      } else {
        cur[j] = std::max(j == 0 ? prev[0] : std::min({prev[j], prev[j - 1], cur[j - 1]}), d);
      }
    }
    std::swap(prev, cur);
  }
  return prev.back();
}

// The discrete distance between the curves cut into k pieces per edge is at least the
// continuous one, and at most a piece longer: two bounds from a dynamic programme that shares
// nothing with the free-space diagram. Curves that walk back along a line, repeat lattice
// points or zig-zag are where a matching most easily goes wrong.
TEST("frechet/the decision lies between the discrete distances of the curves cut finer") {
  const Frechet<2> f;
  std::mt19937 rng(61);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int trial = 0; trial < 1200; ++trial) {
    const int family = trial % 5;
    const double off = trial % 7 == 0 ? 1e6 : 0.0;
    const auto make = [&](std::size_t n) {
      Curve<2> c;
      for (std::size_t i = 0; i < n; ++i) {
        const double x = static_cast<double>(i);
        switch (family) {
          case 0: c.push_back({off + u(rng), off + u(rng)}); break;
          case 1: c.push_back({off + u(rng), off}); break;
          case 2: c.push_back({off + std::round(3 * u(rng)), off + std::round(3 * u(rng))}); break;
          case 3: c.push_back({off + x + 0.3 * u(rng), off + 0.3 * u(rng)}); break;
          default: c.push_back({off + (i % 2) + 0.01 * u(rng), off + 0.01 * u(rng)}); break;
        }
      }
      return c;
    };
    const Curve<2> a = make(2 + trial % 6);
    Curve<2> b = make(2 + (trial / 6) % 6);
    if (family >= 3) {
      b.front() = a.front();
      b.back() = a.back();
    }
    const Curve<2> fine_a = subdivide(a, 16);
    const Curve<2> fine_b = subdivide(b, 16);
    const double upper = discrete_frechet(fine_a, fine_b);
    const double lower = upper - std::max(longest_edge(fine_a), longest_edge(fine_b));
    const double slack = 1e-9 * std::max(longest_edge(a), longest_edge(b)) + 1e-14 * off;
    const std::string where = "trial " + std::to_string(trial);
    CHECK_MSG(f.within(a, b, upper + slack), where + ": not within the upper bound");
    CHECK_MSG(lower - slack <= 0.0 || !f.within(a, b, lower - slack),
              where + ": within under the lower bound");
    const double tol = 1e-12;
    const double d = f.distance(a, b, tol);
    for (const double r : {0.9 * d, 0.99 * d, 1.01 * d, 1.1 * d}) {
      if (r < d && r >= d - tol) {
        continue;
      }
      const Curve<2> back_a(a.rbegin(), a.rend());
      const Curve<2> back_b(b.rbegin(), b.rend());
      CHECK_MSG(f.within(a, b, r) == f.within(b, a, r), where + ": not symmetric");
      CHECK_MSG(f.within(a, b, r) == f.within(back_a, back_b, r), where + ": not reversible");
      CHECK_MSG(f.within(a, b, r) == (r > d), where + ": the decision and the value disagree");
    }
  }
}

}  // namespace
