#include "simplify/convex_shadow_32_sidenorms.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "simplify/convex_shadow_32.hpp"

namespace ssk::simplify {
namespace {

using validate::Polygon;
using validate::Vec;

Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.y + b.y}; }
Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y}; }
Vec operator*(Vec a, double s) { return {a.x * s, a.y * s}; }
double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y; }
double cross(Vec a, Vec b) { return a.x * b.y - a.y * b.x; }
double norm(Vec a) { return std::sqrt(dot(a, a)); }

// The m-gon around the origin: the edge from corners[k] is on dot(normals[k], x) = offsets[k].
struct Region {
  Polygon corners;
  std::vector<Vec> normals;
  std::vector<double> offsets;
};

// S[p]: dot(normals[k], x) <= offset[k] for every k, and between the rays from p along
// rays[0] (right) and rays[1] (left) once p is beyond a side. `tangents` touch S[p] and are the
// rays of the next one. Rebuilding S[p] needs the rays that cut it: the tangents describe the same
// set, but a tangent can run almost along an edge, and cutting by it again moves the corner there.
struct Candidate {
  Vec p;
  std::vector<double> offset;
  std::array<Vec, 2> rays;
  std::array<Vec, 2> tangents;
};

Polygon clip(const Polygon& poly, Vec n, double c) {
  Polygon out;
  out.reserve(poly.size() + 1);
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec cur = poly[i];
    const Vec nxt = poly[(i + 1) % poly.size()];
    const double sc = c - dot(n, cur);
    const double sn = c - dot(n, nxt);
    if (sc >= 0.0) {
      out.push_back(cur);
    }
    if ((sc > 0.0 && sn < 0.0) || (sc < 0.0 && sn > 0.0)) {
      out.push_back(cur + (nxt - cur) * (sc / (sc - sn)));
    }
  }
  return out;
}

Region moved(Region region, Vec v) {
  for (std::size_t k = 0; k < region.corners.size(); ++k) {
    region.corners[k] = v + region.corners[k];
    region.offsets[k] += dot(region.normals[k], v);
  }
  return region;
}

// The corners of S[p], whose offsets were last taken against `region`. A stricter offset is a
// cut, not just an edge: it can pass the corners of the sides next to it.
Polygon corners(const Candidate& c, const Region& region) {
  Polygon poly = region.corners;
  for (std::size_t k = 0; k < region.offsets.size(); ++k) {
    if (c.offset[k] < region.offsets[k]) {
      poly = clip(poly, region.normals[k], c.offset[k]);
    }
  }
  if (c.rays[0] != Vec{}) {
    for (const Vec d : {c.rays[0], c.rays[1] * -1.0}) {
      const Vec n{d.y, -d.x};
      poly = clip(poly, n, dot(n, c.p));
    }
  }
  return poly;
}

// The rays through the corners furthest right and left seen from p, once p is beyond a side.
// Each side of the inward normal of the side p is furthest beyond is ordered on its own: across
// it the order is noise when p is an ulp beyond.
void aim(Candidate& c, const Polygon& shadow, const std::vector<Vec>& normals) {
  Vec ahead;
  double furthest = 0.0;
  for (std::size_t k = 0; k < normals.size(); ++k) {
    if (dot(normals[k], c.p) - c.offset[k] > furthest) {
      furthest = dot(normals[k], c.p) - c.offset[k];
      ahead = normals[k] * -1.0;
    }
  }
  if (furthest > 0.0) {
    const auto further_left = [p = c.p, ahead](Vec a, Vec b) {
      const bool a_left = cross(ahead, a - p) > 0.0;
      const bool b_left = cross(ahead, b - p) > 0.0;
      return a_left == b_left ? cross(a - p, b - p) > 0.0 : b_left;
    };
    const auto [right, left] = std::minmax_element(shadow.begin(), shadow.end(), further_left);
    c.tangents = {*right - c.p, *left - c.p};
  }
}

// S[p] <- P_v ∩ F(S[p], p). Going away from p only goes deeper into a side p is beyond or on, so
// F keeps that offset and the strictest holds; every other offset is P_v's. Once beyond a side,
// p stays beyond it: the offset only falls.
std::vector<Candidate> advance(const std::vector<Candidate>& cands, const Region& before,
                               const Region& at, validate::Validator* validator) {
  std::vector<Candidate> out;
  for (const Candidate& c : cands) {
    Candidate next{c.p, at.offsets, c.tangents, {}};
    for (std::size_t k = 0; k < at.offsets.size(); ++k) {
      if (dot(at.normals[k], c.p) >= c.offset[k]) {
        next.offset[k] = std::min(c.offset[k], at.offsets[k]);
      }
    }
    const Polygon shadow = corners(next, at);
    // Under three corners is a tangency, without area. Lemma 2's witness always has area.
    if (shadow.size() < 3) {
      continue;
    }
    if (validator != nullptr) {
      validator->shadow(shadow, at.corners, c.p,
                        std::isinf(c.offset[0]) ? Polygon{c.p} : corners(c, before));
    }
    aim(next, shadow, at.normals);
    out.push_back(std::move(next));
  }
  return out;
}

// The point of poly nearest q. Only clearly inside counts, as that point becomes an output
// vertex; an edge a few ulps long, left where a cut passes a corner, points nowhere and is skipped.
Vec closest(const Polygon& poly, Vec q) {
  bool inside = true;
  Vec best = poly[0];
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec ab = poly[(i + 1) % poly.size()] - a;
    const double length = norm(ab);
    const double reach = norm(q - a) + norm(q - a - ab);
    inside = inside && (length < 1e-9 * reach || cross(ab, q - a) >= 1e-9 * length * reach);
    const double t = length > 0.0 ? std::clamp(dot(q - a, ab) / dot(ab, ab), 0.0, 1.0) : 0.0;
    const Vec on = a + ab * t;
    if (dot(q - on, q - on) < dot(q - best, q - best)) {
      best = on;
    }
  }
  return inside ? q : best;
}

void emit(Curve<2>& out, Vec q) {
  const Point<2> pt{q.x, q.y};
  if (out.empty() || out.back() != pt) {
    out.push_back(pt);
  }
}

// `last` is the m-gon around the last vertex stabbed, which the shadows were last cut from.
void emit_segment(Curve<2>& out, const std::vector<Candidate>& cands, const Region& last,
                  std::span<const Vec> stabbed, double r, validate::Validator* validator) {
  const Vec first = stabbed.front();
  const Candidate& best = *std::min_element(
      cands.begin(), cands.end(), [first](const Candidate& a, const Candidate& b) {
        return dot(a.p - first, a.p - first) < dot(b.p - first, b.p - first);
      });
  const Vec end = closest(corners(best, last), stabbed.back());
  if (validator != nullptr) {
    validator->segment(best.p, end, stabbed, r);
  }
  emit(out, best.p);
  emit(out, end);
}

// The grid points of the cells of G_v, side h, that meet the circle of radius rho around v.
Polygon boundary_points(double rho, double h) {
  const double radius = rho / h;
  const auto n = static_cast<long long>(std::ceil(radius));
  const auto near = [](long long k) { return std::clamp(0LL, k, k + 1); };
  const auto far = [](long long k) { return std::max(k * k, (k + 1) * (k + 1)); };
  std::vector<std::pair<long long, long long>> ids;
  for (long long i = -n; i < n; ++i) {
    for (long long j = -n; j < n; ++j) {
      const auto near2 = static_cast<double>(near(i) * near(i) + near(j) * near(j));
      const auto far2 = static_cast<double>(far(i) + far(j));
      if (near2 <= radius * radius && far2 > radius * radius) {
        ids.insert(ids.end(), {{i, j}, {i + 1, j}, {i, j + 1}, {i + 1, j + 1}});
      }
    }
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  Polygon out;
  for (const auto& [i, j] : ids) {
    out.push_back(Vec{static_cast<double>(i), static_cast<double>(j)} * h);
  }
  return out;
}

Region regular(int m, double r) {
  Region out;
  for (int k = 0; k < m; ++k) {
    const double corner = 2.0 * std::numbers::pi * k / m;
    const double normal = std::numbers::pi * (2 * k + 1) / m;
    out.corners.push_back({r * std::cos(corner), r * std::sin(corner)});
    out.normals.push_back({std::cos(normal), std::sin(normal)});
    out.offsets.push_back(dot(out.normals.back(), out.corners.back()));
  }
  return out;
}

std::vector<Candidate> start(Vec v, const Polygon& starts, std::size_t m) {
  std::vector<Candidate> out;
  for (const Vec o : starts) {
    out.push_back({v + o, std::vector<double>(m, std::numeric_limits<double>::infinity()), {}, {}});
  }
  return out;
}

}  // namespace

ConvexShadow32SideNorms::ConvexShadow32SideNorms(const Params& params,
                                                 validate::Validator* validator)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)), validator_(validator) {
  if (delta_ <= 0.0) {
    throw std::runtime_error("convex-shadow-32-sidenorms: delta must be positive");
  }
  if (eps_ <= 0.0 || eps_ >= 1.0) {
    throw std::runtime_error("convex-shadow-32-sidenorms: eps must lie in (0, 1)");
  }
}

Curve<2> ConvexShadow32SideNorms::run(const Curve<2>& curve, const Context&) const {
  // Ours reaches exactly (1 + eps) * delta, and writing a point back in the input's coordinates
  // can round it an ulp of them further, so aim a few ulps inside.
  double scale = 0.0;
  for (const Point<2>& q : curve) {
    scale = std::max({scale, std::abs(q[0]), std::abs(q[1])});
  }
  const double delta = delta_ - 16.0 * std::numeric_limits<double>::epsilon() * scale;
  if (delta <= 0.0) {
    return curve;  // delta is below the coordinates' precision: no invented point is safe
  }
  const double r = (1.0 + eps_) * delta;
  const double rho = (1.0 + eps_ / 2.0) * delta;
  const Polygon starts = boundary_points(rho, eps_ * delta / (2.0 * std::sqrt(2.0)));
  // A few ulps inside r, as every test near the error bound errs inward: cos and sin round.
  const Region region = regular(ConvexShadow32::sides(eps_),
                                r * (1.0 - 4.0 * std::numeric_limits<double>::epsilon()));
  if (validator_ != nullptr) {
    validator_->region(region.corners, rho, r);
  }

  // Relative to the first vertex, so that coordinates near zero keep more of delta's digits.
  const Vec origin{curve[0][0], curve[0][1]};
  std::vector<Vec> vertices;
  for (const Point<2>& q : curve) {
    vertices.push_back(Vec{q[0], q[1]} - origin);
  }

  Curve<2> out;
  std::size_t first = 0;
  std::vector<Candidate> cands = start(vertices[0], starts, region.normals.size());
  Region before = moved(region, vertices[0]);
  for (std::size_t i = 1; i < vertices.size(); ++i) {
    Region at = moved(region, vertices[i]);
    std::vector<Candidate> next = advance(cands, before, at, validator_);
    if (next.empty()) {
      emit_segment(out, cands, before, std::span(vertices).subspan(first, i - first), r,
                   validator_);
      first = i;
      next = start(vertices[i], starts, region.normals.size());
    }
    cands = std::move(next);
    before = std::move(at);
  }
  if (first + 1 == vertices.size()) {
    emit(out, vertices.back());
  } else {
    emit_segment(out, cands, before, std::span(vertices).subspan(first), r, validator_);
  }

  for (Point<2>& q : out) {
    q = {q[0] + origin.x, q[1] + origin.y};
  }
  return out;
}

}  // namespace ssk::simplify
