// Every greedy stabber keeps the guarantee, d_F <= (1 + eps) * delta, decided exactly, and
// Algorithm 2 only shares Algorithm 1's computations, so its output is Algorithm 1's to the bit.
#include "frechet/frechet.hpp"
#include "simplify/greedy_stabbing.hpp"
#include "test_main.hpp"

#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
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
}

}  // namespace
