#include "simplify/convex_shadow_32.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

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
double dist2(Vec a, Vec b) { return dot(a - b, a - b); }

struct Candidate {
  Vec p;
  Polygon shadow;
};

Vec closest_on_segment(Vec q, Vec a, Vec b) {
  const Vec ab = b - a;
  const double len2 = dot(ab, ab);
  if (len2 <= 0.0) {
    return a;
  }
  return a + ab * std::clamp(dot(q - a, ab) / len2, 0.0, 1.0);
}

// Only clearly inside counts: a point taken as inside becomes an output vertex, so near the
// boundary the answer has to err inward.
bool contains(const Polygon& poly, Vec q) {
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec a = poly[i];
    const Vec b = poly[(i + 1) % poly.size()];
    const double scale = norm(b - a) * (norm(q - a) + norm(q - b)) + 1e-300;
    if (cross(b - a, q - a) < 1e-9 * scale) {
      return false;
    }
  }
  return true;
}

Vec closest(const Polygon& poly, Vec q) {
  if (poly.size() >= 3 && contains(poly, q)) {
    return q;
  }
  Vec best = poly[0];
  double best_d2 = dist2(q, best);
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec c = closest_on_segment(q, poly[i], poly[(i + 1) % poly.size()]);
    if (dist2(q, c) < best_d2) {
      best_d2 = dist2(q, c);
      best = c;
    }
  }
  return best;
}

Polygon clip_half_plane(const Polygon& poly, Vec a, Vec b) {
  Polygon out;
  out.reserve(poly.size() + 2);
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Vec cur = poly[i];
    const Vec nxt = poly[(i + 1) % poly.size()];
    const double sc = cross(b - a, cur - a);
    const double sn = cross(b - a, nxt - a);
    if (sc >= 0.0) {
      out.push_back(cur);
    }
    if ((sc > 0.0 && sn < 0.0) || (sc < 0.0 && sn > 0.0)) {
      out.push_back(cur + (nxt - cur) * (sc / (sc - sn)));
    }
  }
  return out;
}

Polygon hull(Polygon pts) {
  std::sort(pts.begin(), pts.end(),
            [](Vec a, Vec b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
  if (pts.size() <= 2) {
    return pts;
  }
  Polygon h;
  h.reserve(2 * pts.size());
  const auto add = [&h](Vec q, std::size_t keep) {
    while (h.size() >= keep && cross(h.back() - h[h.size() - 2], q - h[h.size() - 2]) <= 0.0) {
      h.pop_back();
    }
    h.push_back(q);
  };
  for (const Vec q : pts) {
    add(q, 2);
  }
  const std::size_t lower = h.size() + 1;
  for (auto it = pts.rbegin() + 1; it != pts.rend(); ++it) {
    add(*it, lower);
  }
  h.pop_back();
  return h;
}

// Whether edge a->b of a shadow faces p: p is beyond its line, or too near it to say. Too near
// is within 1e-9 of the lengths involved, where the sign of cross(b - a, p - a) is rounding's
// to decide, and a p in line with every edge of a thin shadow would otherwise be inside it.
bool faces(Vec a, Vec b, Vec p) {
  const double side = cross(b - a, p - a);
  return side < 0.0 || side * side < 1e-18 * dist2(a, b) * (dist2(p, a) + dist2(p, b));
}

// around ∩ F(shadow, p), where F is every point seen from p behind the shadow. F is bounded by
// the edges of the shadow that face p and by the two tangents from p, which touch the shadow
// where a facing edge meets one facing away; the edges facing away are the side it is open on.
// Every step here only cuts `around`, so taking an edge to face p when it does not costs area
// and never the bound.
Polygon behind(const Polygon& around, const Polygon& shadow, Vec p) {
  // A vertex of the shadow is in it, and then every point is behind it. The equality is exact,
  // which matters when a start is a vertex of the region it is cut from.
  if (std::find(shadow.begin(), shadow.end(), p) != shadow.end()) {
    return around;
  }
  Polygon out = around;
  std::size_t facing = 0;
  bool outside = false;
  Vec left;
  Vec right;
  bool previous = faces(shadow.back(), shadow.front(), p);
  for (std::size_t i = 0; i < shadow.size(); ++i) {
    const Vec a = shadow[i];
    const Vec b = shadow[(i + 1) % shadow.size()];
    const bool current = faces(a, b, p);
    if (current != previous) {
      (current ? left : right) = a;
    }
    previous = current;
    if (current) {
      ++facing;
      outside = outside || cross(b - a, p - a) < 0.0;
      out = clip_half_plane(out, a, b);
    }
  }
  if (facing == 0) {
    return out;
  }
  // From on an edge the tangents are its own line, so they are cut only from beyond one. With
  // every edge facing p there is nowhere for them to touch, and the cuts above have left no
  // more than the shadow, which is in F.
  if (outside && facing < shadow.size()) {
    out = clip_half_plane(out, p, right);
    out = clip_half_plane(out, left, p);
  }
  // Clipping leaves repeated points and edges a few ulps long whose direction is noise. Every
  // edge of a hull has the rest of it on its left, so the next step can believe each one.
  return hull(std::move(out));
}

// The grid points of the cells of G_v that meet its boundary circle, with v at the origin and
// cell side h. Only these start a segment: sliding a start backwards along its segment until
// it leaves the ball keeps every stab in order, so a boundary start reaches as far as an
// interior one.
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

// A regular m-gon of circumradius r around the origin, counter-clockwise.
Polygon regular(int m, double r) {
  Polygon out;
  for (int i = 0; i < m; ++i) {
    const double a = 2.0 * std::numbers::pi * i / m;
    out.push_back({r * std::cos(a), r * std::sin(a)});
  }
  return out;
}

std::vector<Candidate> start(Vec v, const Polygon& starts) {
  std::vector<Candidate> out;
  for (const Vec o : starts) {
    out.push_back({v + o, {v + o}});
  }
  return out;
}

// S[p] <- P_v ∩ F(S[p], p), where P_v is the m-gon around v.
std::vector<Candidate> advance(const std::vector<Candidate>& cands, Vec v, const Polygon& region,
                               validate::Validator* validator) {
  Polygon around;
  for (const Vec o : region) {
    around.push_back(v + o);
  }
  std::vector<Candidate> out;
  for (const Candidate& c : cands) {
    Polygon next = behind(around, c.shadow, c.p);
    // Under three vertices is a tangency, a shadow without area, whose edges cannot say where p
    // is. Lemma 2's witness stabs balls strictly inside P_v, so its shadow always has area.
    if (next.size() >= 3) {
      if (validator != nullptr) {
        validator->shadow(next, around, c.p, c.shadow);
      }
      out.push_back({c.p, std::move(next)});
    }
  }
  return out;
}

void emit(Curve<2>& out, Vec q) {
  const Point<2> pt{q.x, q.y};
  if (out.empty() || out.back() != pt) {
    out.push_back(pt);
  }
}

void emit_segment(Curve<2>& out, const std::vector<Candidate>& cands,
                  std::span<const Vec> stabbed, double r, validate::Validator* validator) {
  const Vec first = stabbed.front();
  const Candidate& best = *std::min_element(
      cands.begin(), cands.end(),
      [first](const Candidate& a, const Candidate& b) { return dist2(a.p, first) < dist2(b.p, first); });
  const Vec end = closest(best.shadow, stabbed.back());
  if (validator != nullptr) {
    validator->segment(best.p, end, stabbed, r);
  }
  emit(out, best.p);
  emit(out, end);
}

}  // namespace

ConvexShadow32::ConvexShadow32(const Params& params, validate::Validator* validator)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)), validator_(validator) {
  if (delta_ <= 0.0) {
    throw std::runtime_error("convex-shadow-32: delta must be positive");
  }
  if (eps_ <= 0.0 || eps_ >= 1.0) {
    throw std::runtime_error("convex-shadow-32: eps must lie in (0, 1)");
  }
}

int ConvexShadow32::sides(double eps) {
  const double ratio = (1.0 + eps / 2.0) / (1.0 + eps);  // rho / r
  int m = std::max(3, static_cast<int>(std::ceil(std::numbers::pi / std::acos(ratio))));
  while (std::cos(std::numbers::pi / m) < ratio) {
    ++m;
  }
  while (m > 3 && std::cos(std::numbers::pi / (m - 1)) >= ratio) {
    --m;
  }
  return m;
}

Curve<2> ConvexShadow32::run(const Curve<2>& curve, const Context&) const {
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
  const double h = eps_ * delta / (2.0 * std::sqrt(2.0));
  const double r = (1.0 + eps_) * delta;
  const double rho = (1.0 + eps_ / 2.0) * delta;
  const Polygon starts = boundary_points(rho, h);
  // A few ulps inside r, as every test near the error bound errs inward: cos and sin round.
  const Polygon region =
      regular(sides(eps_), r * (1.0 - 4.0 * std::numeric_limits<double>::epsilon()));
  if (validator_ != nullptr) {
    validator_->region(region, rho, r);
  }

  // The region moves with each vertex, so working relative to the first one changes nothing
  // but the rounding: coordinates near zero keep more of delta's digits.
  const Vec origin{curve[0][0], curve[0][1]};
  std::vector<Vec> vertices;
  for (const Point<2>& q : curve) {
    vertices.push_back(Vec{q[0], q[1]} - origin);
  }

  Curve<2> out;
  std::size_t first = 0;
  std::vector<Candidate> cands = start(vertices[0], starts);
  for (std::size_t i = 1; i < vertices.size(); ++i) {
    std::vector<Candidate> next = advance(cands, vertices[i], region, validator_);
    if (next.empty()) {
      emit_segment(out, cands, std::span(vertices).subspan(first, i - first), r, validator_);
      first = i;
      next = start(vertices[i], starts);
    }
    cands = std::move(next);
  }
  if (first + 1 == vertices.size()) {
    emit(out, vertices.back());
  } else {
    emit_segment(out, cands, std::span(vertices).subspan(first), r, validator_);
  }

  for (Point<2>& q : out) {
    q = {q[0] + origin.x, q[1] + origin.y};
  }
  return out;
}

}  // namespace ssk::simplify
