#include "simplify/greedy_stabber/greedy_stabber_decomp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "geometry/geometry.hpp"
#include "simplify/greedy_stabber/greedy_stabber.hpp"

namespace ssk::simplify {
namespace {

using geometry::Base;
using geometry::ConvexPolygon;
using geometry::Polygon;
using geometry::Vec;
using geometry::Wedge;

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

  std::vector<NaiveCandidate> start(double scale) override {
    std::vector<NaiveCandidate> out;
    for (const Vec o : starts_) {
      out.push_back({o * scale, Wedge{},
                     std::vector<double>(region_.normals.size(),
                                         std::numeric_limits<double>::infinity())});
    }
    return out;
  }

  std::vector<NaiveCandidate> advance(const std::vector<NaiveCandidate>& cands, Vec v,
                                      double scale) override {
    const ConvexPolygon at = placed(region_, scale, v);
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

  void emit_segment(Curve<2>& out, Vec last, double scale) override {
    const ConvexPolygon at = placed(region_, scale, last);
    for (const NaiveCandidate& c : candidates_) {
      const Polygon s = shadow(c.p, c.wedge, c.offsets, at);
      if (!s.empty()) {
        emit(out, c.p);
        emit(out, closest(s, last));
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
    return grid_points(delta, eps, false);
  }
};

class NaiveStabberPaperBoundary final : public NaiveStabberPaper {
 protected:
  [[nodiscard]] std::vector<Vec> start_points(double delta, double eps) const override {
    return grid_points(delta, eps, true);
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

  std::vector<SharedCandidate> start(double scale) override {
    scale_ = scale;
    strictest_.assign(region_.normals.size(), std::numeric_limits<double>::infinity());
    std::vector<SharedCandidate> out;
    for (std::size_t s = 0; s < starts_.size(); ++s) {
      out.push_back({s, starts_[s] * scale, Wedge{}});
    }
    return out;
  }

  // Between two marked starts the sides a start is beyond are the same, so are Q[p] and its two
  // corners: they are built at the first start after a mark and kept until the next.
  std::vector<SharedCandidate> advance(const std::vector<SharedCandidate>& cands, Vec v,
                                       double scale) override {
    const ConvexPolygon at = placed(region_, scale, v);
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

  void emit_segment(Curve<2>& out, Vec last, double scale) override {
    const ConvexPolygon at = placed(region_, scale, last);
    for (const SharedCandidate& c : candidates_) {
      const Polygon s = shadow(c.p, c.wedge, strictest_, at);
      if (!s.empty()) {
        emit(out, c.p);
        emit(out, closest(s, last));
        return;
      }
    }
    throw std::logic_error("greedy-stabbing-shared: no surviving start has a shadow to end in");
  }

 private:
  // The starts where the sides they are beyond change from the start before: one past each
  // crossing of the circle with a line n_k . x = strictest_k. The crossing's angle only says where
  // to look. It rounds differently from the dot product that decides beyond in bounds() and
  // base(), so that test itself is asked of the starts a spacing either side, and a start is
  // marked only where its answer differs from the one before. A mark costs a build of Q[p], a
  // test only a dot product.
  [[nodiscard]] std::vector<bool> crossings(const std::vector<double>& strictest) const {
    const auto count = static_cast<long long>(starts_.size());
    const double spacing = 2.0 * std::numbers::pi / static_cast<double>(count);
    const auto wrap = [count](long long s) {
      return static_cast<std::size_t>((s % count + count) % count);
    };
    std::vector<bool> marked(starts_.size(), false);
    for (std::size_t k = 0; k < strictest.size(); ++k) {
      const Vec n = region_.normals[k];
      const double t = strictest[k] / (delta_ * scale_);
      if (std::abs(t) > 1.0 + 1e-9) {
        continue;
      }
      const auto beyond = [&](long long s) {
        return dot(n, starts_[wrap(s)] * scale_) > strictest[k];
      };
      const double half = std::acos(std::clamp(t, -1.0, 1.0));
      for (const double angle : {std::atan2(n.y, n.x) - half, std::atan2(n.y, n.x) + half}) {
        const auto first = static_cast<long long>(std::floor(angle / spacing)) - 1;
        bool before = beyond(first - 1);
        for (long long s = first; s <= first + 3; ++s) {
          const bool here = beyond(s);
          if (here != before) {
            marked[wrap(s)] = true;
          }
          before = here;
        }
      }
    }
    return marked;
  }

  double delta_ = 0.0;
  ConvexPolygon region_;
  std::vector<Vec> starts_;
  double scale_ = 1.0;
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
