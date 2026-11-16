// Both implementations of ours keep the paper's guarantee, d_F <= (1 + eps) * delta, decided
// exactly, and agree on how many vertices that takes.
#include "frechet/frechet.hpp"
#include "simplify/convex_shadow.hpp"
#include "simplify/convex_shadow_far_points.hpp"
#include "test_main.hpp"

#include <cmath>
#include <vector>

namespace {

using ssk::simplify::Curve;
using ssk::simplify::Params;

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

TEST("convex-shadow/both implementations stay within (1 + eps) * delta and agree on size") {
  const ssk::frechet::Frechet<2> frechet;
  for (const Curve<2>& curve : curves()) {
    for (const double eps : {0.5, 0.1}) {
      for (const double delta : {0.05, 0.4, 2.0}) {
        const Params params{{"delta", delta}, {"eps", eps}};
        const Curve<2> edges = ssk::simplify::ConvexShadow(params).simplify(curve);
        const Curve<2> far = ssk::simplify::ConvexShadowFarPoints(params).simplify(curve);
        CHECK(frechet.within(curve, edges, (1.0 + eps) * delta));
        CHECK(frechet.within(curve, far, (1.0 + eps) * delta));
        CHECK(edges.size() == far.size());
      }
    }
  }
}

TEST("convex-shadow/each implementation names itself") {
  const Params params{{"delta", 1.0}};
  CHECK(ssk::simplify::ConvexShadow(params).name() == "convex-shadow");
  CHECK(ssk::simplify::ConvexShadowFarPoints(params).name() == "convex-shadow-far-points");
}

}  // namespace
