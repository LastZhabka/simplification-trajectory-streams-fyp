#include "simplify/greedy_stabber/greedy_stabber_no_decomp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "geometry/geometry.hpp"
#include "simplify/greedy_stabber/greedy_stabber.hpp"

namespace ssk::simplify {
namespace {

using geometry::ConvexPolygon;
using geometry::Polygon;
using geometry::Vec;

// n . x <= b.
struct HalfPlane {
  Vec n;
  double b;
};

// The edge from `from` to the next corner, on the line of `h`, with the polygon inside `h`. An
// edge keeps the half-plane it came from rather than taking its line from the corners, whose
// direction is noise on an edge a few ulps long.
struct Edge {
  Vec from;
  HalfPlane h;
};

// A convex polygon, counter-clockwise.
using Shadow = std::vector<Edge>;

Polygon corners(const Shadow& s) {
  Polygon out;
  for (const Edge& e : s) {
    out.push_back(e.from);
  }
  return out;
}

// As in geometry::base, a shadow without area counts as empty.
bool has_area(const Shadow& s) {
  double area = 0.0;
  for (std::size_t j = 1; j + 1 < s.size(); ++j) {
    area += cross(s[j].from - s[0].from, s[j + 1].from - s[0].from);
  }
  return area > 0.0;
}

// a * b as the double nearest and what rounding left out: Dekker's product, as the build has no
// fused multiply-add.
std::pair<double, double> product(double a, double b) {
  const auto split = [](double x) {
    const double c = 134217729.0 * x;
    const double hi = c - (c - x);
    return std::pair{hi, x - hi};
  };
  const double p = a * b;
  const auto [ah, al] = split(a);
  const auto [bh, bl] = split(b);
  return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
}

// a * d - b * c, to a few ulps of itself even where the two products all but cancel.
double difference(double a, double d, double b, double c) {
  const auto [p, ep] = product(a, d);
  const auto [q, eq] = product(b, c);
  return (p - q) + (ep - eq);
}

// Where the lines of two half-planes meet, to a few ulps however near parallel they are. Two
// consecutive sides of a shadow can be: the tangents from a start a few ulps beyond a side, which
// run along it. geometry::meet's cross products cancel to noise there, and put the corner as far
// off as ulps of the coordinates over the angle between them.
Vec meet(const HalfPlane& a, const HalfPlane& b) {
  const double d = difference(a.n.x, b.n.y, a.n.y, b.n.x);
  return {difference(a.b, b.n.y, b.b, a.n.y) / d, difference(a.n.x, b.b, b.n.x, a.b) / d};
}

// The t for which a + t d is in the first `count` half-planes, or none. A half-plane missed by no
// more than `close`, what the coordinates resolve, counts as met, and so does a range shorter than
// nothing by no more than that, which is a corner where more than two sides meet.
std::optional<std::pair<double, double>> along(const std::vector<HalfPlane>& h, std::size_t count,
                                               Vec a, Vec d, double close) {
  double lo = -std::numeric_limits<double>::infinity();
  double hi = std::numeric_limits<double>::infinity();
  for (std::size_t j = 0; j < count; ++j) {
    const double toward = dot(h[j].n, d);
    const double room = h[j].b - dot(h[j].n, a);
    if (toward > 0.0) {
      hi = std::min(hi, room / toward);
    } else if (toward < 0.0) {
      lo = std::max(lo, room / toward);
    } else if (room < -close) {
      return std::nullopt;
    }
  }
  if (lo > hi + close) {
    return std::nullopt;
  }
  return lo > hi ? std::pair{(lo + hi) / 2.0, (lo + hi) / 2.0} : std::pair{lo, hi};
}

// The point of the half-planes' intersection furthest along u, or none when it is empty. Seidel's
// incremental linear program: the first four half-planes are a box, whose corner furthest along u
// starts it, and the optimum moves only when a half-plane cuts it off by more than `close`, onto
// that half-plane's line, to the furthest point there that the half-planes before it allow. It
// only finds a point to build the hull about, which polygon() checks exactly, so to give way by
// `close` costs nothing there.
std::optional<Vec> furthest(const std::vector<HalfPlane>& h, Vec u, double close) {
  Vec x{u.x >= 0.0 ? h[0].b : -h[1].b, u.y >= 0.0 ? h[2].b : -h[3].b};
  for (std::size_t i = 4; i < h.size(); ++i) {
    if (dot(h[i].n, x) <= h[i].b + close) {
      continue;
    }
    const Vec d{-h[i].n.y, h[i].n.x};
    const Vec a = h[i].n * (h[i].b / dot(h[i].n, h[i].n));
    const std::optional<std::pair<double, double>> t = along(h, i, a, d, close);
    if (!t) {
      return std::nullopt;
    }
    x = a + d * (dot(u, d) >= 0.0 ? t->second : t->first);
  }
  return x;
}

// A point strictly inside the half-planes' intersection, or none when it has no area: the middle
// of its chord halfway between its leftmost and rightmost points. A chord strictly between the
// two runs through the inside, as no side can lie along it.
std::optional<Vec> inside(const std::vector<HalfPlane>& h, double close) {
  const std::optional<Vec> left = furthest(h, {-1.0, 0.0}, close);
  const std::optional<Vec> right = furthest(h, {1.0, 0.0}, close);
  if (!left || !right || !(right->x > left->x)) {
    return std::nullopt;
  }
  const double x = (left->x + right->x) / 2.0;
  const std::optional<std::pair<double, double>> chord =
      along(h, h.size(), {x, 0.0}, {0.0, 1.0}, close);
  if (!chord || !(chord->second > chord->first)) {
    return std::nullopt;
  }
  const Vec c{x, (chord->first + chord->second) / 2.0};
  for (const HalfPlane& p : h) {
    if (!(dot(p.n, c) < p.b)) {
      return std::nullopt;
    }
  }
  return c;
}

// The half-planes, by index, whose sides bound their intersection, counter-clockwise. About a
// point c strictly inside it, n . x <= b is d . (x - c) <= 1 with d = n / (b - n . c), so the
// intersection is the polar of the convex hull of the d, and its sides are the hull's corners, in
// the hull's order: the half-plane intersection of the paper's Lemma 1(ii), by its convex hull
// [Chazelle 1993]. The hull is Andrew's monotone chain; a d inside it is a half-plane that does not
// bound the intersection, such as a side of S facing p that is no stricter than the region's.
std::vector<std::size_t> bounding(const std::vector<HalfPlane>& h, Vec c) {
  std::vector<Vec> d;
  for (const HalfPlane& p : h) {
    d.push_back(p.n * (1.0 / (p.b - dot(p.n, c))));
  }
  std::vector<std::size_t> order(h.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
    return d[i].x < d[j].x || (d[i].x == d[j].x && d[i].y < d[j].y);
  });
  const auto turns_left = [&](std::size_t a, std::size_t b, std::size_t q) {
    return cross(d[b] - d[a], d[q] - d[a]) > 0.0;
  };
  std::vector<std::size_t> hull(2 * order.size());
  std::size_t k = 0;
  for (const std::size_t i : order) {
    while (k >= 2 && !turns_left(hull[k - 2], hull[k - 1], i)) {
      --k;
    }
    hull[k++] = i;
  }
  const std::size_t lower = k + 1;
  for (std::size_t r = order.size() - 1; r-- > 0;) {
    while (k >= lower && !turns_left(hull[k - 2], hull[k - 1], order[r])) {
      --k;
    }
    hull[k++] = order[r];
  }
  hull.resize(k - 1);
  return hull;
}

// The polygon of the half-planes, the first four a box around it, or nothing when it has no area.
Shadow polygon(const std::vector<HalfPlane>& h, double close) {
  const std::optional<Vec> c = inside(h, close);
  if (!c) {
    return {};
  }
  const std::vector<std::size_t> sides = bounding(h, *c);
  if (sides.size() < 3) {
    return {};
  }
  Shadow out;
  for (std::size_t j = 0; j < sides.size(); ++j) {
    const HalfPlane& here = h[sides[j]];
    const HalfPlane& next = h[sides[(j + 1) % sides.size()]];
    out.push_back({meet(here, next), next});
  }
  return has_area(out) ? out : Shadow{};
}

// A_v as half-planes, the first four a box of half-size `reach` around v, and as a polygon; and
// `close`, 16 ulps of the coordinates, what they resolve. That is within the margin run() leaves
// inside delta.
struct Around {
  std::vector<HalfPlane> h;
  Shadow at;
  Vec v;
  double reach;
  double close;
};

struct Candidate {
  Vec p;
  Shadow shadow;
  // p was outside its shadow, so it is outside for good: p not in S puts p outside F(S, p). The
  // sides alone cannot say so on a sliver a few ulps wide between the two tangents, whose closing
  // side rounding has lost: none of its sides then faces p.
  bool outside = false;
};

// A_v ∩ F(S, p), where F(S, p) is the points y whose segment py meets S, as Lemma 1(ii) of the
// paper builds it. If p is in S, F is the plane and the step is A_v itself. Otherwise F is bounded
// by the sides of S that face p, the support lines that separate p from S, and the two tangents
// from p, through the corners where those sides stop and start; with the box and A_v's sides,
// their intersection is the new shadow. A tangent's b is n . p itself, so p is on it to the bit
// and it never faces p at a later step. An empty shadow is a start that dies.
Candidate step(const Candidate& c, const Around& around) {
  const Shadow& s = c.shadow;
  const Vec p = c.p;
  std::vector<bool> faces;
  std::size_t facing = 0;
  for (const Edge& e : s) {
    faces.push_back(dot(e.h.n, p) > e.h.b);
    facing += faces.back() ? 1 : 0;
  }
  if (facing == 0) {
    return {p, c.outside ? Shadow{} : around.at, c.outside};
  }
  if (facing == s.size()) {
    return {p, {}, true};
  }
  Vec right;
  Vec left;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const std::size_t next = (i + 1) % s.size();
    if (faces[i] && !faces[next]) {
      right = s[next].from;
    }
    if (!faces[i] && faces[next]) {
      left = s[next].from;
    }
  }
  const Vec r = right - p;
  const Vec l = left - p;
  // A cone of zero width is a tangency, as in geometry::intersect: the witness of the size bound
  // always leaves area.
  if (cross(r, l) <= 0.0) {
    return {p, {}, true};
  }
  // Two half-planes on one line, as far as the coordinates resolve it anywhere in the box, are one:
  // the one stricter at v is kept. They are common, a side of S and of A_v where a vertex repeats,
  // a tangent and a side of A_v where a start of the grid lies on that side, and kept apart they
  // leave `meet` a corner far off. The box is left as it is.
  std::vector<HalfPlane> h = around.h;
  const Vec v = around.v;
  const double close = around.close;
  const auto add = [&](HalfPlane q) {
    const auto same = std::find_if(h.begin() + 4, h.end(), [&](const HalfPlane& o) {
      const Vec dn = q.n - o.n;
      const double tilt = around.reach * (std::abs(dn.x) + std::abs(dn.y));
      return std::abs(dot(dn, v) - (q.b - o.b)) + tilt <= close;
    });
    if (same == h.end()) {
      h.push_back(q);
    } else if (q.b - dot(q.n, v) < same->b - dot(same->n, v)) {
      *same = q;
    }
  };
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (faces[i]) {
      add(s[i].h);
    }
  }
  const Vec nr = Vec{r.y, -r.x} * (1.0 / norm(r));
  const Vec nl = Vec{-l.y, l.x} * (1.0 / norm(l));
  add({nr, dot(nr, p)});
  add({nl, dot(nl, p)});
  return {p, polygon(h, close), true};
}

class NoDecompStabber : public GreedyStabber<Candidate> {
 protected:
  // The region around a vertex and the starts around a run's first one, with it at the origin.
  [[nodiscard]] virtual ConvexPolygon region(double delta, double eps) const {
    return mgon(delta, eps);
  }
  [[nodiscard]] virtual std::vector<Vec> start_points(double delta, double eps) const {
    return circle(delta, eps);
  }

  void init(double delta, double eps) override {
    region_ = region(delta, eps);
    starts_ = start_points(delta, eps);
    reach_ = 2.0 * (1.0 + eps) * delta;
  }

  // S[p] = {p}: a shadow with no sides, none of which faces p, so its first step is the region.
  std::vector<Candidate> start(Vec v) override {
    std::vector<Candidate> out;
    for (const Vec o : starts_) {
      out.push_back({v + o, {}, false});
    }
    return out;
  }

  // Every region lies within (1 + eps) delta of its vertex, so the box at twice that bounds every
  // shadow, and none of its sides bounds one.
  std::vector<Candidate> advance(const std::vector<Candidate>& cands, std::span<const Vec>,
                                 Vec v) override {
    Around around{{{{1.0, 0.0}, v.x + reach_},
                   {{-1.0, 0.0}, reach_ - v.x},
                   {{0.0, 1.0}, v.y + reach_},
                   {{0.0, -1.0}, reach_ - v.y}},
                  {},
                  v,
                  reach_,
                  16.0 * std::numeric_limits<double>::epsilon() *
                      (std::abs(v.x) + std::abs(v.y) + reach_)};
    const ConvexPolygon a = translated(region_, v);
    for (std::size_t k = 0; k < a.normals.size(); ++k) {
      around.h.push_back({a.normals[k], a.offsets[k]});
    }
    around.at = polygon(around.h, around.close);
    std::vector<Candidate> out;
    for (const Candidate& c : cands) {
      Candidate next = step(c, around);
      if (!next.shadow.empty()) {
        out.push_back(std::move(next));
      }
    }
    return out;
  }

  void emit_segment(Curve<2>& out, std::span<const Vec> stabbed) override {
    for (const Candidate& c : candidates_) {
      if (!c.shadow.empty()) {
        emit(out, c.p);
        emit(out, closest(corners(c.shadow), stabbed.back()));
        return;
      }
    }
    throw std::logic_error("no-decomp greedy stabber: no surviving start has a shadow");
  }

 private:
  ConvexPolygon region_;
  std::vector<Vec> starts_;
  double reach_ = 0.0;
};

class NoDecompStabberPaper : public NoDecompStabber {
 protected:
  [[nodiscard]] ConvexPolygon region(double delta, double eps) const override {
    return grid_hull(delta, eps);
  }
  [[nodiscard]] std::vector<Vec> start_points(double delta, double eps) const override {
    return grid_points(delta, eps, false);
  }
};

class NoDecompStabberPaperBoundary final : public NoDecompStabberPaper {
 protected:
  [[nodiscard]] std::vector<Vec> start_points(double delta, double eps) const override {
    return grid_points(delta, eps, true);
  }
};

}  // namespace

GreedyStabbingNoDecomp::GreedyStabbingNoDecomp(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingNoDecomp::run(const Curve<2>& curve, const Context&) const {
  return NoDecompStabber().run(curve, delta_, eps_);
}

GreedyStabbingNoDecompPaper::GreedyStabbingNoDecompPaper(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingNoDecompPaper::run(const Curve<2>& curve, const Context&) const {
  return NoDecompStabberPaper().run(curve, delta_, eps_);
}

GreedyStabbingNoDecompPaperBoundary::GreedyStabbingNoDecompPaperBoundary(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingNoDecompPaperBoundary::run(const Curve<2>& curve, const Context&) const {
  return NoDecompStabberPaperBoundary().run(curve, delta_, eps_);
}

}  // namespace ssk::simplify
