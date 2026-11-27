// Every greedy stabber keeps the guarantee, d_F <= (1 + eps) * delta, decided exactly, and
// Algorithm 2 only shares Algorithm 1's computations, so its output is Algorithm 1's to the bit.
#include "frechet/frechet.hpp"
#include "simplify/greedy_stabber/greedy_stabber_decomp.hpp"
#include "simplify/greedy_stabber/greedy_stabber_no_decomp.hpp"
#include "test_main.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using ssk::simplify::Curve;
using ssk::simplify::Params;
using ssk::simplify::Simplifier;

std::vector<Curve<2>> curves() {
  std::vector<Curve<2>> out(3);
  for (int i = 0; i < 300; ++i) {
    const double t = 0.05 * i;
    out[0].push_back({t * std::cos(t), t * std::sin(t)});
    out[1].push_back({0.5 * i, (i % 2 == 0) ? 0.0 : 1.2});
    out[2].push_back({10.0 * std::sin(0.11 * i) + 0.01 * i, 7.0 * std::sin(0.07 * i)});
  }
  return out;
}

std::vector<std::unique_ptr<Simplifier<2>>> stabbers(const Params& params) {
  std::vector<std::unique_ptr<Simplifier<2>>> out;
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingNaive>(params));
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingNaivePaper>(params));
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingNaivePaperBoundary>(params));
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingShared>(params));
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingNoDecomp>(params));
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingNoDecompPaper>(params));
  out.push_back(std::make_unique<ssk::simplify::GreedyStabbingNoDecompPaperBoundary>(params));
  return out;
}

TEST("greedy-stabbing/every stabber stays within (1 + eps) * delta") {
  const ssk::frechet::Frechet<2> frechet;
  for (const Curve<2>& curve : curves()) {
    for (const double eps : {0.5, 0.1}) {
      for (const double delta : {0.05, 0.4, 2.0}) {
        for (const auto& s : stabbers({{"delta", delta}, {"eps", eps}})) {
          CHECK_MSG(frechet.within(curve, s->simplify(curve), (1.0 + eps) * delta),
                    s->name() + " at eps " + std::to_string(eps) + ", delta " +
                        std::to_string(delta));
        }
      }
    }
  }
}

// A, B, A with B where the m-gons around A and B all but touch corner to corner: the shadow is a
// sliver near the one point within (1 + eps) * delta of both, far from A. Found at eps 0.5 in
// Mopsi and MOT20, whose grid steps run along the hexagon's corners.
TEST("greedy-stabbing/a spike to where two corners meet stays within (1 + eps) * delta") {
  const ssk::frechet::Frechet<2> frechet;
  const double pi = std::acos(-1.0);
  for (const auto& [eps, corner] : {std::pair{0.5, -pi / 2}, std::pair{0.2, pi / 8}}) {
    for (const double gap : {1e-9, 1e-11, 1e-13}) {
      const double r = 2.0 * (1.0 + eps) * (1.0 - gap);
      const Curve<2> spike{{0, 0}, {r * std::cos(corner), r * std::sin(corner)}, {0, 0}};
      for (const auto& s : stabbers({{"delta", 1.0}, {"eps", eps}})) {
        CHECK_MSG(frechet.within(spike, s->simplify(spike), (1.0 + eps)),
                  s->name() + " at eps " + std::to_string(eps) + ", gap " + std::to_string(gap));
      }
    }
  }
  const Curve<2> mopsi{{29.7773, 62.6202}, {29.7773, 62.6201}, {29.7773, 62.6202}};
  const double delta = 3.333333356009115e-05;
  for (const auto& s : stabbers({{"delta", delta}, {"eps", 0.5}})) {
    CHECK_MSG(frechet.within(mopsi, s->simplify(mopsi), 1.5 * delta), s->name() + " on mopsi");
  }
}

// The same spike far from the first vertex, its gap just over run()'s precision margin. With the
// margin taken from the whole curve the shadow was a few ulps of its coordinates wide, closest()
// skipped every edge of it, and A was emitted, a delta from B: 1.72 delta against 1.5. The margin
// is now each vertex's, ulps of |o| + 2 delta + 64 |v - o| for the run's first vertex o, here the
// vertex before the spike.
TEST("greedy-stabbing/the spike far from the first vertex, just over the margin, keeps the bound") {
  const ssk::frechet::Frechet<2> frechet;
  const double pi = std::acos(-1.0);
  const double ulp = std::numeric_limits<double>::epsilon();
  for (const auto& [eps, corner] : {std::pair{0.5, -pi / 2}, std::pair{0.2, pi / 8}}) {
    for (const double far : {1e2, 1e4, 1e6}) {
      for (const double over : {1.01, 1.1, 1.3, 1.5}) {
        const double gap = over * ulp * (far + 2.0 + 64.0 * (1.0 + eps));
        const double r = 2.0 * (1.0 + eps) * (1.0 - gap);
        const double x = far - 10.0;
        const Curve<2> spike{{0, 0},  {far - 20.0, 0},
                             {x, 0},  {x, 10},
                             {x + r * std::cos(corner), 10 + r * std::sin(corner)},
                             {x, 10}};
        for (const auto& s : stabbers({{"delta", 1.0}, {"eps", eps}})) {
          CHECK_MSG(frechet.within(spike, s->simplify(spike), 1.0 + eps),
                    s->name() + " at eps " + std::to_string(eps) + ", " + std::to_string(far) +
                        " from the first vertex, " + std::to_string(over) + " x the margin");
        }
      }
    }
  }
}

// MOT20's mot-001866 at eps 0.25: a staircase of half-pixel steps, after two points that fix the
// origin and the precision margin and end the first run. A clip left a start's shadow a sliver
// between its two tangents, without the side that closed it, so no side faced p; the no-decomp
// step took p back into its shadow, gave the start the whole region, and went over the bound.
TEST("greedy-stabbing/a sliver between the tangents does not put p back in its shadow") {
  const ssk::frechet::Frechet<2> frechet;
  const Curve<2> stairs{{599.5, 1070.5}, {599.5, 721.5}, {840, 840},   {840, 840.5}, {839, 840.5},
                        {839, 841},      {838, 841},     {838, 841.5}, {838, 842}};
  const double delta = 0.26764416432842186;
  for (const auto& s : stabbers({{"delta", delta}, {"eps", 0.25}})) {
    CHECK_MSG(frechet.within(stairs, s->simplify(stairs), 1.25 * delta), s->name());
  }
}

TEST("greedy-stabbing/sharing the computations changes nothing") {
  for (const Curve<2>& curve : curves()) {
    for (const double eps : {0.5, 0.3, 0.1, 0.05}) {
      for (const double delta : {0.05, 0.4, 2.0}) {
        const Params params{{"delta", delta}, {"eps", eps}};
        CHECK(ssk::simplify::GreedyStabbingShared(params).simplify(curve) ==
              ssk::simplify::GreedyStabbingNaive(params).simplify(curve));
      }
    }
  }
}

// The shared stabber builds Q[p] only at the starts it marks, so a missed mark is an output that
// differs from the naive one. Lattice walks, also at coordinates in the millions, where the angle
// that says where to look and the dot product that decides round furthest apart.
TEST("greedy-stabbing/sharing changes nothing on lattice walks, at any offset") {
  std::mt19937 rng(7);
  for (const double offset : {0.0, 6451925.0}) {
    for (int walk = 0; walk < 6; ++walk) {
      Curve<2> curve;
      double x = offset, y = offset / 3.0;
      for (int i = 0; i < 200; ++i) {
        curve.push_back({x, y});
        x += static_cast<double>(static_cast<int>(rng() % 5) - 1);
        y += static_cast<double>(static_cast<int>(rng() % 5) - 2);
      }
      for (const double eps : {0.5, 0.2, 0.1, 0.05}) {
        for (const double delta : {0.5, 1.0, 2.0, 3.0}) {
          const Params params{{"delta", delta}, {"eps", eps}};
          CHECK_MSG(ssk::simplify::GreedyStabbingShared(params).simplify(curve) ==
                        ssk::simplify::GreedyStabbingNaive(params).simplify(curve),
                    "walk " + std::to_string(walk) + " at offset " + std::to_string(offset) +
                        ", eps " + std::to_string(eps) + ", delta " + std::to_string(delta));
        }
      }
    }
  }
}

// The paper's own step keeps S[p] as a polygon, Algorithm 1 as offsets and a cone, and the two
// are the same set, so on the same region and starts only rounding tells them apart.
TEST("greedy-stabbing/the paper's own step agrees with Algorithm 1 on each region and starts") {
  for (const Curve<2>& curve : curves()) {
    for (const double eps : {0.5, 0.2}) {
      for (const double delta : {0.05, 0.4, 2.0}) {
        const Params params{{"delta", delta}, {"eps", eps}};
        std::vector<std::unique_ptr<Simplifier<2>>> naive;
        naive.push_back(std::make_unique<ssk::simplify::GreedyStabbingNaive>(params));
        naive.push_back(std::make_unique<ssk::simplify::GreedyStabbingNaivePaper>(params));
        naive.push_back(std::make_unique<ssk::simplify::GreedyStabbingNaivePaperBoundary>(params));
        std::vector<std::unique_ptr<Simplifier<2>>> no_decomp;
        no_decomp.push_back(std::make_unique<ssk::simplify::GreedyStabbingNoDecomp>(params));
        no_decomp.push_back(std::make_unique<ssk::simplify::GreedyStabbingNoDecompPaper>(params));
        no_decomp.push_back(
            std::make_unique<ssk::simplify::GreedyStabbingNoDecompPaperBoundary>(params));
        for (std::size_t k = 0; k < naive.size(); ++k) {
          const Curve<2> a = naive[k]->simplify(curve);
          const Curve<2> b = no_decomp[k]->simplify(curve);
          bool close = a.size() == b.size();
          for (std::size_t i = 0; close && i < a.size(); ++i) {
            close = std::hypot(a[i][0] - b[i][0], a[i][1] - b[i][1]) <= 1e-9 * delta;
          }
          CHECK_MSG(close, no_decomp[k]->name() + " at eps " + std::to_string(eps) +
                               ", delta " + std::to_string(delta) + ": " +
                               std::to_string(a.size()) + " and " + std::to_string(b.size()) +
                               " vertices");
        }
      }
    }
  }
}

// MOT20's mot-001866 at eps 0.25, after two points that fix the origin and the precision margin
// and end the first run: a start of the grid on a side of the region, whose tangent is that side
// to the ulp. Kept as two lines, they left the paper's step a corner far off or no point inside,
// 1.6 delta from Algorithm 1; as one line, it is Algorithm 1 to the bit.
TEST("greedy-stabbing/the paper's step takes a tangent along a side of the region as that side") {
  const Curve<2> steps{
      {599.5, 1070.5}, {599.5, 721.5}, {643, 1039}, {643.5, 1038}, {643.5, 1037.5}};
  const Params params{{"delta", 0.26764416432842186}, {"eps", 0.25}};
  const Curve<2> a = ssk::simplify::GreedyStabbingNaivePaperBoundary(params).simplify(steps);
  const Curve<2> b = ssk::simplify::GreedyStabbingNoDecompPaperBoundary(params).simplify(steps);
  bool close = a.size() == b.size();
  for (std::size_t i = 0; close && i < a.size(); ++i) {
    close = std::hypot(a[i][0] - b[i][0], a[i][1] - b[i][1]) <= 1e-9 * 0.26764416432842186;
  }
  CHECK(close);
}

TEST("greedy-stabbing/a straight line is one segment, a point repeated at most one") {
  Curve<2> line;
  Curve<2> still;
  for (int i = 0; i < 100; ++i) {
    line.push_back({0.3 * i, -0.1 * i});
    still.push_back({4.0, 2.0});
  }
  for (const auto& s : stabbers({{"delta", 0.2}, {"eps", 0.3}})) {
    CHECK_MSG(s->simplify(line).size() == 2, s->name() + " on the line");
    CHECK_MSG(s->simplify(still).size() <= 2, s->name() + " on the repeated point");
  }
}

// At State Plane feet a coordinate's ulp is 9e-10, so a delta of 1e-10 cannot place a point;
// step 7 asked for such deltas wherever a baseline's error was 0 within measurement.
TEST("greedy-stabbing/below what its coordinates resolve a stabber drops exact repeats only") {
  const ssk::frechet::Frechet<2> frechet;
  const double x0 = 6451203.729;
  const double y0 = 1873252.549;
  const Curve<2> curve{{x0, y0},       {x0, y0},       {x0, y0},       {x0 + 3.1, y0},
                       {x0 + 3.1, y0}, {x0 + 3.1, y0 + 2.0}, {x0 + 3.1, y0 + 2.0},
                       {x0 + 6.0, y0 + 2.5}};
  const Curve<2> no_repeats{{x0, y0}, {x0 + 3.1, y0}, {x0 + 3.1, y0 + 2.0}, {x0 + 6.0, y0 + 2.5}};
  const Curve<2> still(7, {x0, y0});
  for (const auto& s : stabbers({{"delta", 1e-10}, {"eps", 0.5}})) {
    CHECK_MSG(s->simplify(curve) == no_repeats, s->name());
    CHECK_MSG(frechet.within(curve, s->simplify(curve), 0.0), s->name() + " is exact");
    CHECK_MSG(s->simplify(still).size() == 2, s->name() + " on a point repeated");
  }
}

// NGSIM's ngsim-002153: runs of fixes along one axis, collinear to the bit, which Douglas-Peucker
// drops at error 0, so step 9 asked for delta = its measurement tolerance, 3.6 ulps of the
// coordinates. Ours fell back to dropping repeats, 170 vertices against Theorem 3's 2 * 61 - 2.
TEST("greedy-stabbing/a few ulps above what its coordinates resolve, a stabber still stabs") {
  const ssk::frechet::Frechet<2> frechet;
  const double x0 = 6451925.53;
  const double y0 = 1872806.873;
  Curve<2> corner;
  for (int i = 0; i < 20; ++i) {
    corner.push_back({x0 - 0.016 * (20 - i), y0});
  }
  for (int i = 0; i < 20; ++i) {
    corner.push_back({x0, y0 + 0.008 * i});
  }
  const double delta = 3.389505126375869e-09;
  for (const double eps : {0.5, 0.05}) {
    for (const auto& s : stabbers({{"delta", delta}, {"eps", eps}})) {
      const Curve<2> out = s->simplify(corner);
      CHECK_MSG(out.size() <= 2 * 3 - 2, s->name() + " at eps " + std::to_string(eps) + ": " +
                                             std::to_string(out.size()) + " vertices");
      CHECK_MSG(frechet.within(corner, out, (1.0 + eps) * delta),
                s->name() + " at eps " + std::to_string(eps));
    }
  }
}

TEST("greedy-stabbing/at coordinates in the millions, just above that, the bound holds") {
  const ssk::frechet::Frechet<2> frechet;
  Curve<2> curve;
  for (int i = 0; i < 200; ++i) {
    curve.push_back({6451203.729 + 2e-8 * i, 1873252.549 + 5e-8 * std::sin(i / 3.0)});
  }
  for (const double delta : {1.6e-9, 2e-9, 3.4e-9, 1e-8, 2.5e-8, 1e-7}) {
    for (const double eps : {0.5, 0.1}) {
      for (const auto& s : stabbers({{"delta", delta}, {"eps", eps}})) {
        CHECK_MSG(frechet.within(curve, s->simplify(curve), (1.0 + eps) * delta),
                  s->name() + " at delta " + std::to_string(delta) + ", eps " +
                      std::to_string(eps));
      }
    }
  }
}

TEST("greedy-stabbing/delta must be positive and eps in (0, 1)") {
  for (const Params& bad : {Params{{"delta", 0.0}}, Params{{"delta", 1.0}, {"eps", 0.0}},
                            Params{{"delta", 1.0}, {"eps", 1.0}}}) {
    bool threw = false;
    try {
      (void)ssk::simplify::GreedyStabbingShared(bad);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
  }
}

TEST("greedy-stabbing/each stabber names itself") {
  const auto all = stabbers({{"delta", 1.0}});
  CHECK(all[0]->name() == "greedy-stabbing-naive");
  CHECK(all[1]->name() == "greedy-stabbing-naive-paper");
  CHECK(all[2]->name() == "greedy-stabbing-naive-paper-boundary");
  CHECK(all[3]->name() == "greedy-stabbing-shared");
  CHECK(all[4]->name() == "greedy-stabbing-no-decomp");
  CHECK(all[5]->name() == "greedy-stabbing-no-decomp-paper");
  CHECK(all[6]->name() == "greedy-stabbing-no-decomp-paper-boundary");
}

Curve<2> lattice_walk(std::mt19937& rng, std::size_t n) {
  Curve<2> c;
  double x = 0.0;
  double y = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({x, y});
    x += static_cast<double>(static_cast<int>(rng() % 5) - 1);
    y += static_cast<double>(static_cast<int>(rng() % 5) - 2);
  }
  return c;
}

// The pipeline builds one simplifier per match, but nothing says it must: a stabber keeps
// nothing of one curve into the next.
TEST("greedy-stabbing/a simplifier run on one curve and another gives the first the same again") {
  std::mt19937 rng(109);
  const Curve<2> a = lattice_walk(rng, 80);
  const Curve<2> b = lattice_walk(rng, 50);
  for (const auto& s : stabbers({{"delta", 1.3}, {"eps", 0.3}})) {
    const Curve<2> first = s->simplify(a);
    (void)s->simplify(b);
    CHECK_MSG(s->simplify(a) == first, s->name());
  }
}

// run() works relative to each run's first vertex, which for integer coordinates is exact, so moving
// the curve moves the output, but for the precision margin, which grows with the coordinates.
TEST("greedy-stabbing/moving an integer curve by whole units moves the output with it") {
  std::mt19937 rng(113);
  const Curve<2> here = lattice_walk(rng, 70);
  const double dx = 1048576.0;
  const double dy = -524288.0;
  Curve<2> there;
  for (const auto& p : here) {
    there.push_back({p[0] + dx, p[1] + dy});
  }
  const double delta = 0.73;
  for (const auto& s : stabbers({{"delta", delta}, {"eps", 0.3}})) {
    const Curve<2> a = s->simplify(here);
    const Curve<2> b = s->simplify(there);
    bool moved = a.size() == b.size();
    for (std::size_t i = 0; moved && i < a.size(); ++i) {
      moved = std::hypot(b[i][0] - dx - a[i][0], b[i][1] - dy - a[i][1]) <= 1e-6 * delta;
    }
    CHECK_MSG(moved, s->name() + ": " + std::to_string(a.size()) + " and " +
                         std::to_string(b.size()) + " vertices");
  }
}

// A stream does not know what comes next, so nothing emitted may depend on it. The margin was
// taken from the whole curve's largest coordinate and extent, and a suffix far away moved every
// vertex before it. Only the last run, which the prefix ends early, may differ.
TEST("greedy-stabbing/what is emitted for a prefix does not depend on what follows it") {
  std::mt19937 rng(127);
  const Curve<2> prefix = lattice_walk(rng, 60);
  Curve<2> whole = prefix;
  for (const auto& p : lattice_walk(rng, 40)) {
    whole.push_back({p[0] + 3.0e6, p[1] - 1.0e6});
  }
  for (const auto& s : stabbers({{"delta", 0.61}, {"eps", 0.3}})) {
    const Curve<2> a = s->simplify(prefix);
    const Curve<2> b = s->simplify(whole);
    const std::size_t settled = a.size() - 2;
    CHECK_MSG(b.size() >= settled && std::equal(a.begin(), a.begin() + settled, b.begin()),
              s->name() + ": " + std::to_string(a.size()) + " and " + std::to_string(b.size()) +
                  " vertices");
  }
}

// A run is worked out relative to its first vertex, so the margin grows along it with how far the
// vertex is from there, 64 ulps of it. At the end of a line 10 long that is 1.4e-13, so at delta
// 1e-9 the line is still one segment.
TEST("greedy-stabbing/the margin grows along a run, and a long line is still one segment") {
  const ssk::frechet::Frechet<2> frechet;
  Curve<2> line;
  for (int i = 0; i <= 1000; ++i) {
    line.push_back({0.01 * i, 0.0});
  }
  const double delta = 1e-9;
  for (const double eps : {0.5, 0.1}) {
    for (const auto& s : stabbers({{"delta", delta}, {"eps", eps}})) {
      const Curve<2> out = s->simplify(line);
      CHECK_MSG(out.size() == 2, s->name() + ": " + std::to_string(out.size()) + " vertices");
      CHECK_MSG(frechet.within(line, out, (1.0 + eps) * delta),
                s->name() + " at eps " + std::to_string(eps));
    }
  }
}

// At delta 1e-12 a vertex 100 from its run's first one has a margin of 64 ulps of 100, more than
// delta: no point fits within delta of it in that run, so it is output itself, and the run before
// it is one vertex, which a no-decomp stabber, emitting from a shadow, cannot end. So the curve
// comes back whole.
TEST("greedy-stabbing/a vertex whose margin is delta or more is output itself") {
  const Curve<2> steps{{0, 0}, {100, 0}, {200, 0}, {200, 100}};
  for (const auto& s : stabbers({{"delta", 1e-12}, {"eps", 0.5}})) {
    CHECK_MSG(s->simplify(steps) == steps, s->name());
  }
}

// A run is written back from its own first vertex, so a last run of one vertex, after a turn back
// no segment can follow, is that vertex to the bit; from the curve's first it was q - o + o.
TEST("greedy-stabbing/a last run of one vertex ends the output on that vertex exactly") {
  const Curve<2> jump{{0.1, 0.7}, {10.1, 0.7}, {20.1, 0.7}, {-0.3, -0.3}};
  for (const auto& s : stabbers({{"delta", 0.05}, {"eps", 0.3}})) {
    const Curve<2> out = s->simplify(jump);
    CHECK_MSG(out.back() == jump.back(), s->name());
  }
}

// The margin grows with the coordinates of a run's first vertex. Past what delta resolves no point
// fits within delta of a vertex, so from there each vertex is output itself, while the line before
// them, near the origin, is still one segment: the margin of the whole curve kept all 70.
TEST("greedy-stabbing/where a stream's coordinates outgrow delta it keeps its vertices from there") {
  const ssk::frechet::Frechet<2> frechet;
  Curve<2> curve;
  for (int i = 0; i < 50; ++i) {
    curve.push_back({0.5 * i, 0.25 * i});
  }
  for (int i = 0; i < 20; ++i) {
    curve.push_back({1e9 + 3.0 * i, 1e9 - 2.0 * (i % 3)});
  }
  const double delta = 1e-7;
  for (const double eps : {0.5, 0.1}) {
    for (const auto& s : stabbers({{"delta", delta}, {"eps", eps}})) {
      const Curve<2> out = s->simplify(curve);
      CHECK_MSG(out.size() == 2 + 20, s->name() + " at eps " + std::to_string(eps) + ": " +
                                          std::to_string(out.size()) + " vertices");
      CHECK_MSG(std::equal(curve.end() - 20, curve.end(), out.end() - 20), s->name());
      CHECK_MSG(frechet.within(curve, out, (1.0 + eps) * delta),
                s->name() + " at eps " + std::to_string(eps));
    }
  }
}

// A spiral out from 1 to 2^20: the margin grows with the coordinates, so the regions shrink from
// vertex to vertex and the starts from run to run. Every stabber must take them alike: the bound
// holds, shared is naive to the bit, and each no-decomp one is its naive twin.
TEST("greedy-stabbing/a stream that moves out over many binades keeps every guarantee") {
  const ssk::frechet::Frechet<2> frechet;
  Curve<2> spiral;
  for (int i = 0; i < 200; ++i) {
    const double r = std::exp2(i / 10.0);
    spiral.push_back({r * std::cos(0.7 * i), r * std::sin(0.7 * i)});
  }
  const double delta = 0.3;
  const double eps = 0.5;
  const auto all = stabbers({{"delta", delta}, {"eps", eps}});
  std::vector<Curve<2>> out;
  for (const auto& s : all) {
    out.push_back(s->simplify(spiral));
    CHECK_MSG(frechet.within(spiral, out.back(), (1.0 + eps) * delta), s->name());
  }
  CHECK(out[3] == out[0]);
  for (const auto& [twin, no_decomp] : {std::pair{0, 4}, std::pair{1, 5}, std::pair{2, 6}}) {
    const Curve<2>& a = out[static_cast<std::size_t>(twin)];
    const Curve<2>& b = out[static_cast<std::size_t>(no_decomp)];
    bool close = a.size() == b.size();
    for (std::size_t i = 0; close && i < a.size(); ++i) {
      close = std::hypot(a[i][0] - b[i][0], a[i][1] - b[i][1]) <= 1e-9 * delta;
    }
    CHECK_MSG(close, all[static_cast<std::size_t>(no_decomp)]->name() + ": " +
                         std::to_string(a.size()) + " and " + std::to_string(b.size()));
  }
}

// Legs 50 delta long and right-angled corners: each corner needs a vertex near it, the corners do,
// so kappa is the corners and the two ends, within delta or within (1 + eps) delta alike. Then
// Theorem 3 bounds ours from both sides by numbers known exactly.
TEST("greedy-stabbing/an L, a U and a staircase take from kappa to 2 kappa - 2 vertices") {
  const auto walk = [](std::vector<std::pair<double, double>> corners) {
    Curve<2> c;
    for (std::size_t k = 0; k + 1 < corners.size(); ++k) {
      for (int i = 0; i < 100; ++i) {
        const double t = i / 100.0;
        c.push_back({corners[k].first + t * (corners[k + 1].first - corners[k].first),
                     corners[k].second + t * (corners[k + 1].second - corners[k].second)});
      }
    }
    c.push_back({corners.back().first, corners.back().second});
    return c;
  };
  struct Shape {
    const char* name;
    Curve<2> curve;
    std::size_t kappa;
  };
  for (const Shape& shape :
       {Shape{"L", walk({{0, 0}, {50, 0}, {50, 50}}), 3},
        Shape{"U", walk({{0, 0}, {50, 0}, {50, 50}, {0, 50}}), 4},
        Shape{"staircase",
              walk({{0, 0}, {50, 0}, {50, 50}, {100, 50}, {100, 100}, {150, 100}, {150, 150}}),
              7}}) {
    for (const double eps : {0.5, 0.2, 0.1}) {
      for (const auto& s : stabbers({{"delta", 1.0}, {"eps", eps}})) {
        if (eps < 0.2 && s->name().find("paper") != std::string::npos) {
          continue;
        }
        const std::size_t k = s->simplify(shape.curve).size();
        CHECK_MSG(k >= shape.kappa && k <= 2 * shape.kappa - 2,
                  s->name() + " on the " + shape.name + " at eps " + std::to_string(eps) + ": " +
                      std::to_string(k) + " vertices");
      }
    }
  }
}

TEST("greedy-stabbing/eps near either end of (0, 1) keeps the bound") {
  const ssk::frechet::Frechet<2> frechet;
  for (Curve<2> curve : curves()) {
    curve.resize(100);
    for (const double delta : {0.05, 0.4}) {
      for (const double eps : {0.01, 0.99}) {
        for (const auto& s : stabbers({{"delta", delta}, {"eps", eps}})) {
          if (eps < 0.05 && s->name().find("paper") != std::string::npos) {
            continue;
          }
          CHECK_MSG(frechet.within(curve, s->simplify(curve), (1.0 + eps) * delta),
                    s->name() + " at eps " + std::to_string(eps) + ", delta " +
                        std::to_string(delta));
        }
      }
    }
  }
}

// The shared stabber marks a start where the sides it is beyond change, found by angle and
// wrapped round the circle. Turning the curve about its first vertex in steps of a fifth of
// the starts' spacing carries every crossing past every start, the last and first among them.
TEST("greedy-stabbing/sharing changes nothing with the curve turned to every start's angle") {
  std::mt19937 rng(127);
  const Curve<2> walk = lattice_walk(rng, 60);
  for (const double eps : {0.5, 0.2}) {
    const auto starts = static_cast<int>(std::ceil(2.0 * std::acos(-1.0) / eps));
    for (int step = 0; step < 5 * starts; ++step) {
      const double turn = 2.0 * std::acos(-1.0) * step / (5.0 * starts) + 1e-3;
      Curve<2> turned;
      for (const auto& p : walk) {
        turned.push_back({p[0] * std::cos(turn) - p[1] * std::sin(turn),
                          p[0] * std::sin(turn) + p[1] * std::cos(turn)});
      }
      const Params params{{"delta", 1.1}, {"eps", eps}};
      CHECK_MSG(ssk::simplify::GreedyStabbingShared(params).simplify(turned) ==
                    ssk::simplify::GreedyStabbingNaive(params).simplify(turned),
                "eps " + std::to_string(eps) + ", turned " + std::to_string(step) + " fifths");
    }
  }
}

// Curves of the kinds the corpus is made of and the kinds that have broken a stabber: walks,
// a lattice like a GPS grid, steps along the axes, a zig-zag at the tangency, runs of repeats,
// a circle, spikes, some at State Plane coordinates.
Curve<2> family(std::mt19937& rng, int kind, std::size_t n, double step, double off) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  Curve<2> c;
  double x = off;
  double y = 0.6 * off;
  double heading = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double k = static_cast<double>(i);
    switch (kind) {
      case 0: x += step * (u(rng) - 0.5); y += step * (u(rng) - 0.5); break;
      case 1: heading += 0.3 * (u(rng) - 0.5); x += step * std::cos(heading);
              y += step * std::sin(heading); break;
      case 2: x = off + step * std::round(10 * u(rng));
              y = 0.6 * off + step * std::round(10 * u(rng)); break;
      case 3: (u(rng) < 0.5 ? x : y) += step * (u(rng) < 0.5 ? 1 : -1); break;
      case 4: x = off + step * static_cast<double>(i % 2);
              y = 0.6 * off + 1e-3 * step * u(rng); break;
      case 5: if (u(rng) < 0.4) { x += step * u(rng); y += step * u(rng); } break;
      case 6: x = off + step * std::cos(0.2 * k); y = 0.6 * off + step * std::sin(0.2 * k); break;
      default: x += step; y = 0.6 * off + (i % 7 == 3 ? 3 * step : 0.0); break;
    }
    c.push_back({x, y});
  }
  return c;
}

// Theorem 3 in both halves. Ours at delta / (1 + eps) is within delta, so kappa(delta) is at most
// its size, and ours at delta has at most twice that less two. Nothing else in the tests checks
// the size bound; the corpus does, after the run.
TEST("greedy-stabbing/within (1 + eps) delta and 2 kappa - 2 vertices, on curves of every kind") {
  const ssk::frechet::Frechet<2> frechet;
  std::mt19937 rng(89);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int trial = 0; trial < 160; ++trial) {
    const double eps = std::vector<double>{0.1, 0.2, 0.3, 0.5}[static_cast<std::size_t>(trial % 4)];
    const double off = trial % 3 == 0 ? 2.2e6 : 0.0;
    const double step = std::pow(10.0, -2.0 + 4.0 * u(rng));
    const Curve<2> curve =
        family(rng, trial % 8, 3 + static_cast<std::size_t>(40 * u(rng)), step, off);
    double delta = step * std::pow(10.0, -1.5 + 2.0 * u(rng));
    if (trial % 8 == 4 && trial % 16 == 4) {
      delta = step / (2.0 * (1.0 + eps)) * (1.0 - std::pow(10.0, -4.0 - 8.0 * u(rng)));
    }
    const std::string where = "trial " + std::to_string(trial) + ", eps " + std::to_string(eps);
    const auto outer = stabbers({{"delta", delta}, {"eps", eps}});
    const auto inner = stabbers({{"delta", delta / (1.0 + eps)}, {"eps", eps}});
    for (std::size_t k = 0; k < outer.size(); ++k) {
      if (eps < 0.2 && outer[k]->name().find("paper") != std::string::npos) {
        continue;
      }
      const Curve<2> out = outer[k]->simplify(curve);
      CHECK_MSG(frechet.within(curve, out, (1.0 + eps) * delta),
                outer[k]->name() + " over the bound, " + where);
      const Curve<2> witness = inner[k]->simplify(curve);
      CHECK_MSG(out.size() <= 2 * std::max<std::size_t>(witness.size(), 2) - 2,
                outer[k]->name() + ": " + std::to_string(out.size()) + " vertices against " +
                    std::to_string(witness.size()) + " within delta, " + where);
    }
    CHECK_MSG(outer[3]->simplify(curve) == outer[0]->simplify(curve),
              "shared is not naive, " + where);
  }
}

}  // namespace
