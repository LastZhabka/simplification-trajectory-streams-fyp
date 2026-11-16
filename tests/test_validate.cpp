// The validator has to fail on what is wrong as well as pass what is right, or a count of
// zero failures says nothing.
#include "simplify/convex_shadow.hpp"
#include "test_main.hpp"
#include "validate/validate.hpp"

#include <cmath>
#include <numbers>
#include <vector>

namespace {

using ssk::validate::Check;
using ssk::validate::Polygon;
using ssk::validate::Validator;
using ssk::validate::Vec;

std::size_t failures(const Validator& v, Check check) {
  return v.stats()[static_cast<std::size_t>(check)].failed;
}

std::size_t checks(const Validator& v, Check check) {
  return v.stats()[static_cast<std::size_t>(check)].checked;
}

const Polygon kSquare{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
const Polygon kPoints{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, 0}, {0, 1}};

TEST("validate/a correct grid hull passes the grid check and the three region checks") {
  Validator v;
  v.grid(kPoints, kSquare, 1.0, 1.5);
  CHECK(ssk::validate::checked(v.stats()) == 4);
  CHECK(ssk::validate::failed(v.stats()) == 0);
}

TEST("validate/a region with a vertex in the middle of an edge fails") {
  // What the hull of scaled grid points keeps when rounding hides that three are in a line.
  const Polygon ring{{-1, -1}, {0, -1}, {1, -1}, {1, 1}, {-1, 1}};
  Validator v;
  v.grid(ring, ring, 0.1, 2.0);
  CHECK(failures(v, Check::RegionConvex) == 1);
  CHECK_MSG(failures(v, Check::GridHull) == 0, "it is still the hull of the points");
}

TEST("validate/a region needs no grid") {
  // A regular hexagon of circumradius 1: inradius cos 30 deg = 0.866.
  Polygon hexagon;
  for (int i = 0; i < 6; ++i) {
    const double a = i * std::numbers::pi / 3;
    hexagon.push_back({std::cos(a), std::sin(a)});
  }
  Validator fits;
  fits.region(hexagon, 0.86, 1.0);
  CHECK(ssk::validate::checked(fits.stats()) == 3);
  CHECK(checks(fits, Check::GridHull) == 0);
  CHECK(ssk::validate::failed(fits.stats()) == 0);

  Validator thin;
  thin.region(hexagon, 0.87, 1.0);
  CHECK_MSG(failures(thin, Check::RegionCoversBall) == 1, "the inradius is 0.866");
}

TEST("validate/a hull that misses a point, turns clockwise or invents a vertex fails") {
  Validator missing;
  missing.grid(kPoints, {{-1, -1}, {1, -1}, {-1, 1}}, 0.1, 2.0);
  CHECK(failures(missing, Check::GridHull) == 1);

  Validator clockwise;
  clockwise.grid(kPoints, {{-1, 1}, {1, 1}, {1, -1}, {-1, -1}}, 0.1, 2.0);
  CHECK(failures(clockwise, Check::RegionConvex) == 1);

  Validator invented;
  invented.grid(kPoints, {{-2, -2}, {2, -2}, {2, 2}, {-2, 2}}, 0.1, 3.0);
  CHECK(failures(invented, Check::GridHull) == 1);
}

TEST("validate/a region too small for the ball or too large for the disk fails") {
  Validator v;
  v.grid(kPoints, kSquare, 1.1, 1.2);
  CHECK(failures(v, Check::GridHull) == 0);
  CHECK_MSG(failures(v, Check::RegionCoversBall) == 1, "the square's sides are 1 from the centre");
  CHECK_MSG(failures(v, Check::RegionInsideDisk) == 1, "its corners are sqrt 2 from the centre");
}

TEST("validate/a shadow seen through the previous one passes") {
  const Polygon grid{{8, -3}, {14, -3}, {14, 3}, {8, 3}};
  const Polygon next{{10, -1}, {12, -1}, {12, 1}, {10, 1}};
  const Polygon previous{{4, -1}, {6, -1}, {6, 1}, {4, 1}};
  Validator v;
  v.shadow(next, grid, {0, 0}, previous);
  CHECK(ssk::validate::checked(v.stats()) == 3);
  CHECK(ssk::validate::failed(v.stats()) == 0);
}

TEST("validate/a shadow not behind the previous one, outside the region or not convex fails") {
  const Polygon grid{{8, -3}, {14, -3}, {14, 3}, {8, 3}};
  const Polygon next{{10, -1}, {12, -1}, {12, 1}, {10, 1}};

  Validator aside;
  aside.shadow(next, grid, {0, 0}, {{4, 5}, {6, 5}, {6, 7}, {4, 7}});
  CHECK(failures(aside, Check::ShadowBehind) == 1);
  CHECK(failures(aside, Check::ShadowConvex) == 0);

  Validator outside;
  outside.shadow({{10, -1}, {16, -1}, {16, 1}, {10, 1}}, grid, {0, 0},
                 {{4, -1}, {6, -1}, {6, 1}, {4, 1}});
  CHECK(failures(outside, Check::ShadowInsideRegion) == 1);

  Validator dented;
  dented.shadow({{10, -1}, {12, -1}, {11, 0}, {12, 1}, {10, 1}}, grid, {0, 0},
                {{4, -2}, {6, -2}, {6, 2}, {4, 2}});
  CHECK(failures(dented, Check::ShadowConvex) == 1);
}

TEST("validate/an edge a few ulps long is not believed") {
  // What clipping leaves behind: a vertex repeated 1e-20 away, so the edge between the two has
  // no meaningful direction. The polygon is still the square, and is treated as one.
  const Polygon grid{{8, -3}, {14, -3}, {14, 3}, {8, 3}};
  const Polygon previous{{4, -1}, {6, -1}, {6, 1}, {6 + 1e-20, 1 - 1e-20}, {4, 1}};
  const Polygon next{{10, -1}, {12, -1}, {12, 1}, {12 - 1e-20, 1 + 1e-20}, {10, 1}};
  Validator v;
  v.shadow(next, grid, {0, 0}, previous);
  CHECK(ssk::validate::failed(v.stats()) == 0);
}

TEST("validate/a start inside its own shadow sees everything") {
  const Polygon grid{{8, -3}, {14, -3}, {14, 3}, {8, 3}};
  Validator v;
  v.shadow(grid, grid, {0, 0}, {{0, 0}});
  CHECK(failures(v, Check::ShadowBehind) == 0);
}

TEST("validate/a segment must meet every disk, and in order") {
  const std::vector<Vec> line{{0, 0.5}, {3, -0.5}, {6, 0.5}, {10, 0}};
  Validator passes;
  passes.segment({0, 0}, {10, 0}, line, 1.0);
  CHECK(failures(passes, Check::SegmentStabs) == 0);

  Validator misses;
  misses.segment({0, 0}, {10, 0}, std::vector<Vec>{{0, 0}, {5, 2}, {10, 0}}, 1.0);
  CHECK_MSG(failures(misses, Check::SegmentStabs) == 1, "the middle vertex is 2 from the segment");

  const std::vector<Vec> back{{0, 0}, {10, 0}, {0, 0}, {10, 0}};
  Validator backwards;
  backwards.segment({0, 0}, {10, 0}, back, 1.0);
  CHECK_MSG(failures(backwards, Check::SegmentStabs) == 1,
            "every disk is met, but not in the order the vertices come");

  Validator wide;
  wide.segment({0, 0}, {10, 0}, back, 11.0);
  CHECK(failures(wide, Check::SegmentStabs) == 0);
}

TEST("validate/convex-shadow passes its own invariants and is unchanged by being watched") {
  using ssk::simplify::Curve;
  std::vector<Curve<2>> curves(3);
  for (int i = 0; i < 300; ++i) {
    const double t = 0.05 * i;
    curves[0].push_back({t * std::cos(t), t * std::sin(t)});
    curves[1].push_back({0.5 * i, (i % 2 == 0) ? 0.0 : 1.2});
    curves[2].push_back({10.0 * std::sin(0.11 * i) + 0.01 * i, 7.0 * std::sin(0.07 * i)});
  }
  for (const Curve<2>& curve : curves) {
    for (const double eps : {0.5, 0.1}) {
      for (const double delta : {0.05, 0.4, 2.0}) {
        const ssk::simplify::Params params{{"delta", delta}, {"eps", eps}};
        Validator v;
        const Curve<2> watched = ssk::simplify::ConvexShadow(params, &v).simplify(curve);
        CHECK(watched == ssk::simplify::ConvexShadow(params).simplify(curve));
        CHECK(checks(v, Check::GridHull) == 1);
        CHECK(checks(v, Check::ShadowBehind) > 0);
        CHECK(checks(v, Check::SegmentStabs) > 0);
        for (std::size_t i = 0; i < ssk::validate::kChecks; ++i) {
          CHECK_MSG(v.stats()[i].failed == 0,
                    std::string(ssk::validate::name(static_cast<Check>(i))) + " failed");
        }
      }
    }
  }
}

}  // namespace
