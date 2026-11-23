#include "simplify/greedy_stabbing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "geometry/geometry.hpp"
#include "simplify/greedy_stabber.hpp"

namespace ssk::simplify {
namespace {

using geometry::Base;
using geometry::ConvexPolygon;
using geometry::Polygon;
using geometry::Vec;
using geometry::Wedge;

void check_params(const std::string& name, double delta, double eps) {
  if (delta <= 0.0) {
    throw std::runtime_error(name + ": delta must be positive");
  }
  if (eps <= 0.0 || eps >= 1.0) {
    throw std::runtime_error(name + ": eps must lie in (0, 1)");
  }
}

// The fewest sides of a regular polygon between the circles of radius (1 + eps/2) and (1 + eps),
// moved off ceil(pi / acos) either way where acos rounds.
std::size_t sides(double eps) {
  const double ratio = (1.0 + eps / 2.0) / (1.0 + eps);
  const auto reaches = [ratio](std::size_t m) {
    return std::cos(std::numbers::pi / static_cast<double>(m)) >= ratio;
  };
  auto m = std::max<std::size_t>(
      3, static_cast<std::size_t>(std::ceil(std::numbers::pi / std::acos(ratio))));
  while (!reaches(m)) {
    ++m;
  }
  while (m > 3 && reaches(m - 1)) {
    --m;
  }
  return m;
}

// A few ulps inside (1 + eps) * delta, as every test near the error bound errs inward: cos and
// sin round.
ConvexPolygon mgon(double delta, double eps) {
  return geometry::regular_polygon(
      sides(eps), (1.0 + eps) * delta * (1.0 - 4.0 * std::numeric_limits<double>::epsilon()));
}

std::vector<Vec> circle(double delta, double eps) {
  const auto count = static_cast<std::size_t>(std::ceil(2.0 * std::numbers::pi / eps));
  std::vector<Vec> out;
  for (std::size_t s = 0; s < count; ++s) {
    const double angle =
        2.0 * std::numbers::pi * static_cast<double>(s) / static_cast<double>(count);
    out.push_back({delta * std::cos(angle), delta * std::sin(angle)});
  }
  return out;
}

// The paper's grid G_v has cells of side eps delta / (2 sqrt 2), with v at a corner.
double cell_side(double delta, double eps) { return eps * delta / (2.0 * std::sqrt(2.0)); }

// The corners of the cells of G_v, in cells around v: of the cells that meet the disk of radius
// (1 + eps/2) delta, or with `boundary` of only those that meet its circle. Nearest v first, so
// the start emitted is the nearest that survived.
std::vector<Vec> grid(double eps, bool boundary) {
  const double radius = (1.0 + eps / 2.0) * 2.0 * std::sqrt(2.0) / eps;
  const auto n = static_cast<long long>(std::ceil(radius));
  const auto near = [](long long k) { return std::clamp(0LL, k, k + 1); };
  const auto far = [](long long k) { return std::max(k * k, (k + 1) * (k + 1)); };
  std::vector<std::array<long long, 3>> corners;
  for (long long i = -n; i < n; ++i) {
    for (long long j = -n; j < n; ++j) {
      const auto near2 = static_cast<double>(near(i) * near(i) + near(j) * near(j));
      const auto far2 = static_cast<double>(far(i) + far(j));
      if (near2 <= radius * radius && (!boundary || far2 > radius * radius)) {
        for (const long long a : {i, i + 1}) {
          for (const long long b : {j, j + 1}) {
            corners.push_back({a * a + b * b, a, b});
          }
        }
      }
    }
  }
  std::sort(corners.begin(), corners.end());
  corners.erase(std::unique(corners.begin(), corners.end()), corners.end());
  std::vector<Vec> out;
  for (const auto& [distance, i, j] : corners) {
    out.push_back({static_cast<double>(i), static_cast<double>(j)});
  }
  return out;
}

std::vector<Vec> scaled(std::vector<Vec> points, double s) {
  for (Vec& q : points) {
    q = q * s;
  }
  return points;
}

// conv(G_v), the paper's region. The hull is taken in cells, where grid points in a line are
// exactly in a line, and scaled after.
ConvexPolygon grid_hull(double delta, double eps) {
  return geometry::from_corners(
      scaled(geometry::full_convex_hull(grid(eps, false), false), cell_side(delta, eps)));
}

// Q[p]: the run's strictest offset where p is beyond it, the region's around v elsewhere.
ConvexPolygon bounds(Vec p, const std::vector<double>& strictest, const ConvexPolygon& at) {
  ConvexPolygon out{at.normals, {}};
  for (std::size_t k = 0; k < at.normals.size(); ++k) {
    const bool beyond = dot(at.normals[k], p) > strictest[k];
    out.offsets.push_back(beyond ? strictest[k] : at.offsets[k]);
  }
  return out;
}

// S[p] = R(S[p], p) ∩ Q[p], with the offsets last taken against the region `at`.
Polygon shadow(Vec p, const Wedge& wedge, const std::vector<double>& strictest,
               const ConvexPolygon& at) {
  return cut(base(bounds(p, strictest, at), p).corners, p, wedge);
}

struct NaiveCandidate {
  Vec p;
  Wedge wedge;
  std::vector<double> offsets;
};

class NaiveStabber : public GreedyStabber<NaiveCandidate> {
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
  }

  std::vector<NaiveCandidate> start(Vec v) override {
    std::vector<NaiveCandidate> out;
    for (const Vec o : starts_) {
      out.push_back({v + o, Wedge{},
                     std::vector<double>(region_.normals.size(),
                                         std::numeric_limits<double>::infinity())});
    }
    return out;
  }

  std::vector<NaiveCandidate> advance(const std::vector<NaiveCandidate>& cands,
                                      std::span<const Vec>, Vec v) override {
    const ConvexPolygon at = translated(region_, v);
    std::vector<NaiveCandidate> out;
    for (const NaiveCandidate& c : cands) {
      NaiveCandidate next = c;
      for (std::size_t k = 0; k < at.offsets.size(); ++k) {
        next.offsets[k] = std::min(c.offsets[k], at.offsets[k]);
      }
      const Base q = base(bounds(c.p, next.offsets, at), c.p);
      if (q.corners.empty()) {
        continue;
      }
      const std::optional<Wedge> wedge = intersect(c.wedge, cone(q, c.p));
      if (!wedge) {
        continue;
      }
      next.wedge = *wedge;
      out.push_back(std::move(next));
    }
    return out;
  }

  void emit_segment(Curve<2>& out, std::span<const Vec> stabbed) override {
    const ConvexPolygon at = translated(region_, stabbed.back());
    for (const NaiveCandidate& c : candidates_) {
      const Polygon s = shadow(c.p, c.wedge, c.offsets, at);
      if (!s.empty()) {
        emit(out, c.p);
        emit(out, closest(s, stabbed.back()));
        return;
      }
    }
    throw std::logic_error("naive greedy stabber: no surviving start has a shadow to end in");
  }

 private:
  ConvexPolygon region_;
  std::vector<Vec> starts_;
};

class NaiveStabberPaper : public NaiveStabber {
 protected:
  [[nodiscard]] ConvexPolygon region(double delta, double eps) const override {
    return grid_hull(delta, eps);
  }
  [[nodiscard]] std::vector<Vec> start_points(double delta, double eps) const override {
    return scaled(grid(eps, false), cell_side(delta, eps));
  }
};

class NaiveStabberPaperBoundary final : public NaiveStabberPaper {
 protected:
  [[nodiscard]] std::vector<Vec> start_points(double delta, double eps) const override {
    return scaled(grid(eps, true), cell_side(delta, eps));
  }
};

struct SharedCandidate {
  std::size_t start;
  Vec p;
  Wedge wedge;
};

class SharedStabber final : public GreedyStabber<SharedCandidate> {
 protected:
  void init(double delta, double eps) override {
    delta_ = delta;
    region_ = mgon(delta, eps);
    starts_ = circle(delta, eps);
  }

  std::vector<SharedCandidate> start(Vec v) override {
    center_ = v;
    strictest_.assign(region_.normals.size(), std::numeric_limits<double>::infinity());
    std::vector<SharedCandidate> out;
    for (std::size_t s = 0; s < starts_.size(); ++s) {
      out.push_back({s, v + starts_[s], Wedge{}});
    }
    return out;
  }

  // Between two marked starts the sides a start is beyond are the same, so are Q[p] and its two
  // corners: they are built at the first start after a mark and kept until the next.
  std::vector<SharedCandidate> advance(const std::vector<SharedCandidate>& cands,
                                       std::span<const Vec>, Vec v) override {
    const ConvexPolygon at = translated(region_, v);
    std::vector<double> strictest(at.offsets.size());
    for (std::size_t k = 0; k < at.offsets.size(); ++k) {
      strictest[k] = std::min(strictest_[k], at.offsets[k]);
    }
    const std::vector<bool> marked = crossings(strictest);
    std::vector<SharedCandidate> out;
    std::optional<Base> shared;
    std::size_t scanned = 0;
    for (const SharedCandidate& c : cands) {
      bool crossed = false;
      for (; scanned <= c.start; ++scanned) {
        crossed = crossed || marked[scanned];
      }
      if (!shared || crossed) {
        shared = base(bounds(c.p, strictest, at), c.p);
      }
      if (shared->corners.empty()) {
        continue;
      }
      const std::optional<Wedge> wedge = intersect(c.wedge, cone(*shared, c.p));
      if (wedge) {
        out.push_back({c.start, c.p, *wedge});
      }
    }
    if (!out.empty()) {
      strictest_ = std::move(strictest);
    }
    return out;
  }

  void emit_segment(Curve<2>& out, std::span<const Vec> stabbed) override {
    const ConvexPolygon at = translated(region_, stabbed.back());
    for (const SharedCandidate& c : candidates_) {
      const Polygon s = shadow(c.p, c.wedge, strictest_, at);
      if (!s.empty()) {
        emit(out, c.p);
        emit(out, closest(s, stabbed.back()));
        return;
      }
    }
    throw std::logic_error("greedy-stabbing-shared: no surviving start has a shadow to end in");
  }

 private:
  // The starts just after a crossing of the circle with a line n_k . x = strictest_k, where the
  // sides they are beyond can change. A crossing is found by angle and a side by dot product, which
  // round differently, so each marks the starts within a quarter spacing of it and the next.
  [[nodiscard]] std::vector<bool> crossings(const std::vector<double>& strictest) const {
    const auto count = static_cast<long long>(starts_.size());
    const double spacing = 2.0 * std::numbers::pi / static_cast<double>(count);
    std::vector<bool> marked(starts_.size(), false);
    for (std::size_t k = 0; k < strictest.size(); ++k) {
      const Vec n = region_.normals[k];
      const double t = (strictest[k] - dot(n, center_)) / delta_;
      if (std::abs(t) > 1.0 + 1e-9) {
        continue;
      }
      const double half = std::acos(std::clamp(t, -1.0, 1.0));
      for (const double angle : {std::atan2(n.y, n.x) - half, std::atan2(n.y, n.x) + half}) {
        const auto first = static_cast<long long>(std::ceil(angle / spacing - 0.25));
        const auto last = static_cast<long long>(std::floor(angle / spacing + 0.25)) + 1;
        for (long long s = first; s <= last; ++s) {
          marked[static_cast<std::size_t>((s % count + count) % count)] = true;
        }
      }
    }
    return marked;
  }

  double delta_ = 0.0;
  ConvexPolygon region_;
  std::vector<Vec> starts_;
  Vec center_;
  std::vector<double> strictest_;
};

}  // namespace

GreedyStabbingNaive::GreedyStabbingNaive(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingNaive::run(const Curve<2>& curve, const Context&) const {
  return NaiveStabber().run(curve, delta_, eps_);
}

GreedyStabbingShared::GreedyStabbingShared(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingShared::run(const Curve<2>& curve, const Context&) const {
  return SharedStabber().run(curve, delta_, eps_);
}

GreedyStabbingNaivePaper::GreedyStabbingNaivePaper(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingNaivePaper::run(const Curve<2>& curve, const Context&) const {
  return NaiveStabberPaper().run(curve, delta_, eps_);
}

GreedyStabbingNaivePaperBoundary::GreedyStabbingNaivePaperBoundary(const Params& params)
    : delta_(params.get("delta")), eps_(params.get("eps", 0.1)) {
  check_params(name(), delta_, eps_);
}

Curve<2> GreedyStabbingNaivePaperBoundary::run(const Curve<2>& curve, const Context&) const {
  return NaiveStabberPaperBoundary().run(curve, delta_, eps_);
}

}  // namespace ssk::simplify
