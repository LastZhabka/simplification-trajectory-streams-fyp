#include "simplify/convex_shadow.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ssk::simplify {
namespace {

constexpr double kPi = 3.14159265358979323846;

struct Vec {
  double x = 0.0;
  double y = 0.0;

  bool operator==(const Vec&) const = default;
};

Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.y + b.y}; }
Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y}; }
Vec operator*(Vec a, double s) { return {a.x * s, a.y * s}; }
double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y; }
double cross(Vec a, Vec b) { return a.x * b.y - a.y * b.x; }
double norm(Vec a) { return std::sqrt(dot(a, a)); }
double dist2(Vec a, Vec b) { return dot(a - b, a - b); }

using Polygon = std::vector<Vec>;

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

// Only clearly inside counts: a point taken as inside may become an output vertex, or widen a
// shadow to all of conv(G_v), so near the boundary the answer has to err inward.
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

Polygon clip(Polygon subject, const Polygon& clipper) {
  for (std::size_t i = 0; i < clipper.size() && !subject.empty(); ++i) {
    subject = clip_half_plane(subject, clipper[i], clipper[(i + 1) % clipper.size()]);
  }
  return subject;
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

// F(shadow, p): every point seen from p behind the shadow, cut off beyond the disk D(v, r).
// Rays at most 30 degrees apart and 1.3 times longer than the disk's far side keep the cut
// clear of the disk, so inside it this is F exactly.
Polygon cover(const Polygon& shadow, Vec p, Vec v, double r) {
  Vec centroid;
  for (const Vec s : shadow) {
    centroid = centroid + s;
  }
  centroid = centroid * (1.0 / static_cast<double>(shadow.size()));
  const double a0 = std::atan2(centroid.y - p.y, centroid.x - p.x);

  double lo = 0.0;
  double hi = 0.0;
  double far = 0.0;
  for (const Vec s : shadow) {
    const Vec d = s - p;
    const double a = std::remainder(std::atan2(d.y, d.x) - a0, 2.0 * kPi);
    lo = std::min(lo, a);
    hi = std::max(hi, a);
    far = std::max(far, norm(d));
  }

  const double reach = 1.3 * std::max(far, norm(p - v) + r) + r;
  const int k = std::max(2, static_cast<int>(std::ceil((hi - lo) / (kPi / 6.0))) + 1);
  Polygon pts = shadow;
  for (int i = 0; i < k; ++i) {
    const double a = a0 + lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(k - 1);
    pts.push_back({p.x + reach * std::cos(a), p.y + reach * std::sin(a)});
  }
  return hull(std::move(pts));
}

// The grid points of the cells of G_v that meet its boundary circle, with v at the origin.
// Only these start a segment: sliding a start backwards along its segment until it leaves
// the ball keeps every stab in order, so a boundary start reaches as far as an interior one.
Polygon boundary_points(double h, double rho) {
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
    out.push_back({static_cast<double>(i) * h, static_cast<double>(j) * h});
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

// S[p] <- conv(G_v) ∩ F(S[p], p), which is all of conv(G_v) once p is in S[p].
std::vector<Candidate> advance(const std::vector<Candidate>& cands, Vec v,
                               const Polygon& cells, double r) {
  Polygon around;
  for (const Vec o : cells) {
    around.push_back(v + o);
  }
  std::vector<Candidate> out;
  for (const Candidate& c : cands) {
    Polygon next = closest(c.shadow, c.p) == c.p ? around : clip(around, cover(c.shadow, c.p, v, r));
    // Under three vertices is a tangency, which cover() would widen to a whole line. Lemma 2's
    // witness stabs balls strictly inside conv(G_v), so its shadow always has area.
    if (next.size() >= 3) {
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

void emit_segment(Curve<2>& out, const std::vector<Candidate>& cands, Vec first, Vec last) {
  const Candidate& best = *std::min_element(
      cands.begin(), cands.end(),
      [first](const Candidate& a, const Candidate& b) { return dist2(a.p, first) < dist2(b.p, first); });
  emit(out, best.p);
  emit(out, closest(best.shadow, last));
}

}  // namespace

ConvexShadow::ConvexShadow(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  if (delta_ <= 0.0) {
    throw std::runtime_error("convex-shadow: delta must be positive");
  }
  if (eps_ <= 0.0 || eps_ >= 1.0) {
    throw std::runtime_error("convex-shadow: eps must lie in (0, 1)");
  }
}

Curve<2> ConvexShadow::run(const Curve<2>& curve, const Context&) const {
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
  const Polygon starts = boundary_points(h, (1.0 + eps_ / 2.0) * delta);
  const Polygon cells = hull(starts);

  // The grid moves with each vertex, so working relative to the first one changes nothing
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
    std::vector<Candidate> next = advance(cands, vertices[i], cells, r);
    if (next.empty()) {
      emit_segment(out, cands, vertices[first], vertices[i - 1]);
      first = i;
      next = start(vertices[i], starts);
    }
    cands = std::move(next);
  }
  if (first + 1 == vertices.size()) {
    emit(out, vertices.back());
  } else {
    emit_segment(out, cands, vertices[first], vertices.back());
  }

  for (Point<2>& q : out) {
    q = {q[0] + origin.x, q[1] + origin.y};
  }
  return out;
}

}  // namespace ssk::simplify
