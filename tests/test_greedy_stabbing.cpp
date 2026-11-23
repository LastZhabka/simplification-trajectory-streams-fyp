// Every greedy stabber keeps the guarantee, d_F <= (1 + eps) * delta, decided exactly, and
// Algorithm 2 only shares Algorithm 1's computations, so its output is Algorithm 1's to the bit.
#include "frechet/frechet.hpp"
#include "simplify/greedy_stabbing.hpp"
#include "test_main.hpp"

#include <cmath>
#include <memory>
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
