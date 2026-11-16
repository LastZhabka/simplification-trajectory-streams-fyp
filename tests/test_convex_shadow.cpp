// Every implementation of ours keeps the paper's guarantee, d_F <= (1 + eps) * delta, decided
// exactly. The two that share conv(G_v) agree on how many vertices that takes.
#include "frechet/frechet.hpp"
#include "simplify/convex_shadow.hpp"
#include "simplify/convex_shadow_32.hpp"
#include "simplify/convex_shadow_far_points.hpp"
#include "test_main.hpp"
#include "validate/validate.hpp"

#include <cmath>
#include <numbers>
#include <string>
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

TEST("convex-shadow/convex-shadow-32 stays within (1 + eps) * delta") {
  const ssk::frechet::Frechet<2> frechet;
  for (const Curve<2>& curve : curves()) {
    for (const double eps : {0.5, 0.1}) {
      for (const double delta : {0.05, 0.4, 2.0}) {
        const Curve<2> out =
            ssk::simplify::ConvexShadow32({{"delta", delta}, {"eps", eps}}).simplify(curve);
        CHECK(frechet.within(curve, out, (1.0 + eps) * delta));
      }
    }
  }
}

TEST("convex-shadow/the m-gon has the fewest sides that reach the inner ball") {
  // Worked by hand in docs/convex-shadow-32.md; conv(G_v) has 12, 16, 28 and 44 corners here.
  CHECK(ssk::simplify::ConvexShadow32::sides(0.5) == 6);
  CHECK(ssk::simplify::ConvexShadow32::sides(0.3) == 7);
  CHECK(ssk::simplify::ConvexShadow32::sides(0.1) == 11);
  CHECK(ssk::simplify::ConvexShadow32::sides(0.05) == 15);
  for (const double eps : {0.9, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001}) {
    const int m = ssk::simplify::ConvexShadow32::sides(eps);
    const double ratio = (1.0 + eps / 2.0) / (1.0 + eps);
    CHECK_MSG(std::cos(std::numbers::pi / m) >= ratio, "m sides reach the inner ball");
    CHECK_MSG(std::cos(std::numbers::pi / (m - 1)) < ratio, "m - 1 do not");
    CHECK_MSG(m <= std::ceil(std::numbers::pi * std::sqrt((1.0 + eps) / eps)),
              "within the O(eps^-1/2) bound at eps = " + std::to_string(eps));
  }
}

TEST("convex-shadow/the m-gon passes the region checks and has no grid hull to check") {
  for (const Curve<2>& curve : curves()) {
    ssk::validate::Validator v;
    const Curve<2> out =
        ssk::simplify::ConvexShadow32({{"delta", 0.4}, {"eps", 0.1}}, &v).simplify(curve);
    CHECK(!out.empty());
    using ssk::validate::Check;
    CHECK(v.stats()[static_cast<std::size_t>(Check::GridHull)].checked == 0);
    CHECK(v.stats()[static_cast<std::size_t>(Check::RegionCoversBall)].checked == 1);
    CHECK(ssk::validate::failed(v.stats()) == 0);
  }
}

TEST("convex-shadow/each implementation names itself") {
  const Params params{{"delta", 1.0}};
  CHECK(ssk::simplify::ConvexShadow(params).name() == "convex-shadow");
  CHECK(ssk::simplify::ConvexShadowFarPoints(params).name() == "convex-shadow-far-points");
  CHECK(ssk::simplify::ConvexShadow32(params).name() == "convex-shadow-32");
}

}  // namespace
