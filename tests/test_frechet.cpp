// The decision procedure is the primitive; the value is bisection on top. Cases here are
// ones whose Frechet distance can be worked out by hand, plus the properties the measure has
// to satisfy for a results table built on it to mean anything.
#include "frechet/frechet.hpp"
#include "test_main.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

namespace {

using ssk::frechet::Curve;
using ssk::frechet::Frechet;

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

}  // namespace
