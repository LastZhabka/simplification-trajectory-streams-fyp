#include "validate/validate.hpp"

#include <algorithm>
#include <cmath>

namespace ssk::validate {
namespace {

Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y}; }
double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y; }
double cross(Vec a, Vec b) { return a.x * b.y - a.y * b.x; }
double length(Vec a) { return std::sqrt(dot(a, a)); }

Vec next_of(const Polygon& poly, std::size_t i) { return poly[(i + 1) % poly.size()]; }

bool convex_ccw(const Polygon& poly) {
  if (poly.size() < 3) {
    return false;
  }
  // A dent is measured against the whole polygon, not the two edges at it: clipping leaves
  // edges a few ulps long, and the turn at one of those is rounding noise.
  const auto perimeter = [&poly] {
    double sum = 0.0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
      sum += length(next_of(poly, i) - poly[i]);
    }
    return sum;
  };
  double area = 0.0;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec b = next_of(poly, i);
    const Vec c = next_of(poly, i + 1);
    area += cross(b - poly[0], c - poly[0]);
    const double turn = cross(b - a, c - b);
    if (turn < 0.0 && turn < -kSlack * perimeter() * std::max(length(b - a), length(c - b))) {
      return false;
    }
  }
  return area > 0.0;
}

// No vertex lies in the line of its neighbours. A start on such a vertex would be in line with
// the edge beside it, and which side of that edge it is on would be rounding noise.
bool all_corners(const Polygon& poly) {
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec b = next_of(poly, i);
    const Vec c = next_of(poly, i + 1);
    if (cross(b - a, c - b) <= kSlack * length(b - a) * length(c - b)) {
      return false;
    }
  }
  return true;
}

// Clipping leaves edges a few ulps long whose direction is rounding noise. Such an edge bounds
// nothing its neighbours do not, so a test against it is skipped rather than believed.
bool negligible(Vec edge, double reach) { return length(edge) < kSlack * reach; }

bool inside_or_on(const Polygon& poly, Vec q) {
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec b = next_of(poly, i);
    const double side = cross(b - a, q - a);
    if (side >= 0.0) {
      continue;
    }
    const double reach = length(q - a) + length(q - b);
    if (!negligible(b - a, reach) && side < -kSlack * length(b - a) * reach) {
      return false;
    }
  }
  return true;
}

double distance_to_segment(Vec q, Vec a, Vec b) {
  const Vec ab = b - a;
  const double len2 = dot(ab, ab);
  const double t = len2 > 0.0 ? std::clamp(dot(q - a, ab) / len2, 0.0, 1.0) : 0.0;
  return length(q - Vec{a.x + t * ab.x, a.y + t * ab.y});
}

// Whether the segment from p to q meets the convex polygon: the part of [0, 1] on the inner
// side of every edge is not empty.
bool segment_meets(Vec p, Vec q, const Polygon& poly) {
  if (poly.size() < 3) {
    return std::any_of(poly.begin(), poly.end(), [&](Vec s) {
      return distance_to_segment(s, p, q) <= kSlack * (length(q - p) + length(s - p));
    });
  }
  double lo = 0.0;
  double hi = 1.0;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec edge = next_of(poly, i) - a;
    double at_p = cross(edge, p - a);
    double at_q = cross(edge, q - a);
    if (at_p >= 0.0 && at_q >= 0.0) {
      continue;
    }
    const double reach = length(p - a) + length(q - a);
    if (negligible(edge, reach)) {
      continue;
    }
    const double slack = kSlack * length(edge) * reach;
    at_p += slack;
    at_q += slack;
    if (at_p < 0.0 && at_q < 0.0) {
      return false;
    }
    if (at_p < 0.0) {
      lo = std::max(lo, at_p / (at_p - at_q));
    } else if (at_q < 0.0) {
      hi = std::min(hi, at_p / (at_p - at_q));
    }
  }
  return lo <= hi;
}

}  // namespace

const char* name(Check check) {
  switch (check) {
    case Check::GridHull: return "grid_hull";
    case Check::RegionConvex: return "region_convex";
    case Check::RegionCoversBall: return "region_covers_ball";
    case Check::RegionInsideDisk: return "region_inside_disk";
    case Check::ShadowConvex: return "shadow_convex";
    case Check::ShadowInsideRegion: return "shadow_inside_region";
    case Check::ShadowBehind: return "shadow_behind";
    case Check::SegmentStabs: return "segment_stabs";
  }
  return "";
}

std::size_t checked(const Stats& stats) {
  std::size_t total = 0;
  for (const Count& count : stats) {
    total += count.checked;
  }
  return total;
}

std::size_t failed(const Stats& stats) {
  std::size_t total = 0;
  for (const Count& count : stats) {
    total += count.failed;
  }
  return total;
}

void add(Stats& into, const Stats& from) {
  for (std::size_t i = 0; i < kChecks; ++i) {
    into[i].checked += from[i].checked;
    into[i].failed += from[i].failed;
  }
}

void Validator::record(Check check, bool ok) {
  Count& count = stats_[static_cast<std::size_t>(check)];
  ++count.checked;
  count.failed += !ok;
}

void Validator::grid(const Polygon& points, const Polygon& hull, double rho, double r) {
  record(Check::GridHull,
         convex_ccw(hull) &&
             std::all_of(points.begin(), points.end(),
                         [&](Vec q) { return inside_or_on(hull, q); }) &&
             std::all_of(hull.begin(), hull.end(), [&](Vec q) {
               return std::find(points.begin(), points.end(), q) != points.end();
             }));
  region(hull, rho, r);
}

void Validator::region(const Polygon& poly, double rho, double r) {
  const bool convex = convex_ccw(poly) && all_corners(poly);
  record(Check::RegionConvex, convex);

  bool covers = convex;
  for (std::size_t i = 0; covers && i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec b = next_of(poly, i);
    covers = cross(a, b) >= rho * (1.0 - kSlack) * length(b - a);
  }
  record(Check::RegionCoversBall, covers);

  record(Check::RegionInsideDisk, std::all_of(poly.begin(), poly.end(), [&](Vec q) {
    return length(q) <= r * (1.0 + kSlack);
  }));
}

void Validator::shadow(const Polygon& next, const Polygon& region, Vec p,
                       const Polygon& previous) {
  record(Check::ShadowConvex, convex_ccw(next));
  record(Check::ShadowInsideRegion, std::all_of(next.begin(), next.end(), [&](Vec q) {
    return inside_or_on(region, q);
  }));
  record(Check::ShadowBehind, std::all_of(next.begin(), next.end(), [&](Vec q) {
    return segment_meets(p, q, previous);
  }));
}

void Validator::segment(Vec a, Vec b, std::span<const Vec> stabbed, double r) {
  const double reach = r * (1.0 + kSlack);
  const Vec ab = b - a;
  const double len2 = dot(ab, ab);
  double at = 0.0;
  bool ok = true;
  for (const Vec v : stabbed) {
    if (len2 <= 0.0) {
      ok = length(v - a) <= reach;
    } else {
      const double closest = dot(v - a, ab) / len2;
      const double off = length(v - Vec{a.x + closest * ab.x, a.y + closest * ab.y});
      const double half = off <= reach ? std::sqrt((reach * reach - off * off) / len2) : -1.0;
      const double enter = std::max(closest - half, 0.0);
      const double leave = std::min(closest + half, 1.0);
      ok = half >= 0.0 && enter <= leave && leave >= at;
      at = std::max(at, enter);
    }
    if (!ok) {
      break;
    }
  }
  record(Check::SegmentStabs, ok);
}

}  // namespace ssk::validate
