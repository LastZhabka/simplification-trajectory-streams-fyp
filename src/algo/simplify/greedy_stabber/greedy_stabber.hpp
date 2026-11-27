// The greedy of ordered stabbing, the same however a stabber keeps its candidates: advance them by
// each vertex, and when none survives, emit a segment from the ones that reached the vertex before
// and start again from this one. A stabber says what a candidate is and how it starts, advances
// and ends a segment; `init` fixes what depends on delta and eps.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "geometry/geometry.hpp"
#include "simplify/simplifier.hpp"

namespace ssk::simplify {

void check_params(const std::string& name, double delta, double eps);

// The write-up's region and starts around a vertex at the origin: the regular m-gon with the
// fewest sides between the circles of radius (1 + eps/2) delta and (1 + eps) delta, and
// M = ceil(2 pi / eps) points on the circle of radius delta.
[[nodiscard]] geometry::ConvexPolygon mgon(double delta, double eps);
[[nodiscard]] std::vector<geometry::Vec> circle(double delta, double eps);

// The paper's G_v around a vertex at the origin: the cells of side eps delta / (2 sqrt 2), with v
// at a corner, that meet the disk of radius (1 + eps/2) delta. Its region conv(G_v), and its grid
// points, or with `boundary` those of only the cells that meet the circle, nearest v first.
[[nodiscard]] geometry::ConvexPolygon grid_hull(double delta, double eps);
[[nodiscard]] std::vector<geometry::Vec> grid_points(double delta, double eps, bool boundary);

// The region built for delta, around v and for scale * delta: the same normals, so only the
// offsets change, one multiplication each.
[[nodiscard]] geometry::ConvexPolygon placed(const geometry::ConvexPolygon& region, double scale,
                                             geometry::Vec v);

template <class Candidate>
class GreedyStabber {
 public:
  virtual ~GreedyStabber() = default;

  Curve<2> run(const Curve<2>& curve, double delta, double eps);

 protected:
  using Vec = geometry::Vec;

  virtual void init(double delta, double eps) = 0;
  // Every size below is the one built for delta times `scale`, which is (delta - Delta) / delta
  // for the margin Delta at that vertex. The run's first vertex is the origin.
  virtual std::vector<Candidate> start(double scale) = 0;
  virtual std::vector<Candidate> advance(const std::vector<Candidate>& cands, Vec v,
                                         double scale) = 0;
  // Called while `candidates_` still holds the ones that reached `last`, the open run's last
  // vertex, whose region had `scale`.
  virtual void emit_segment(Curve<2>& out, Vec last, double scale) = 0;

  // `q` is relative to the open run's first vertex.
  void emit(Curve<2>& out, Vec q) const { keep(out, {q.x + origin_[0], q.y + origin_[1]}); }

  std::vector<Candidate> candidates_;

 private:
  static void keep(Curve<2>& out, const Point<2>& q) {
    if (out.empty() || out.back() != q) {
      out.push_back(q);
    }
  }

  Point<2> origin_{};
};

template <class Candidate>
Curve<2> GreedyStabber<Candidate>::run(const Curve<2>& curve, double delta, double eps) {
  // Ours reaches exactly (1 + eps) * delta, so each vertex's region is built for delta less what
  // the arithmetic can err there, from what has arrived. A run works relative to its first vertex
  // o, so its numbers are as large as the run: writing p and q back rounds by ulps of |o| + 2
  // delta and how far the vertex is from o, and every step errs by ulps of that distance. 64 is a
  // safety factor, not a derived bound; the output is decided exactly by the tests and pipelines.
  init(delta, eps);
  const double ulp = std::numeric_limits<double>::epsilon();
  const auto scale_at = [&](Vec v) {
    const double far = std::max(std::abs(v.x), std::abs(v.y));
    const double from = std::max(std::abs(origin_[0]), std::abs(origin_[1]));
    return 1.0 - ulp * (from + 2.0 * delta + 64.0 * far) / delta;
  };
  Curve<2> out;
  bool open = false;
  // Of the open run only its last vertex is needed, its scale, and whether it is the first.
  Vec last{};
  double last_scale = 1.0;
  bool alone = true;
  // A run of one vertex is that vertex: no start has stepped yet, so a stabber that keeps S[p] has
  // none to emit from.
  const auto close = [&] {
    if (alone) {
      keep(out, origin_);
    } else {
      emit_segment(out, last, last_scale);
    }
    open = false;
  };

  for (std::size_t i = 0; i < curve.size(); ++i) {
    if (open) {
      const Vec v{curve[i][0] - origin_[0], curve[i][1] - origin_[1]};
      const double scale = scale_at(v);
      std::vector<Candidate> next;
      if (scale > 0.0) {
        next = advance(candidates_, v, scale);
      }
      if (!next.empty()) {
        candidates_ = std::move(next);
        last = v;
        last_scale = scale;
        alone = false;
        continue;
      }
      close();
      if (scale <= 0.0) {
        // No point fits within delta of this vertex, but the vertex itself does, exactly.
        keep(out, curve[i]);
        continue;
      }
    }
    origin_ = curve[i];
    const double scale = scale_at({0.0, 0.0});
    if (scale <= 0.0) {
      keep(out, curve[i]);
      continue;
    }
    candidates_ = start(scale);
    last = {0.0, 0.0};
    last_scale = scale;
    alone = true;
    open = true;
  }
  if (open) {
    close();
  }
  // Exact repeats are dropped, so a curve that never moves would be one point.
  if (out.size() == 1 && curve.size() > 1) {
    out.push_back(curve.back());
  }
  return out;
}

}  // namespace ssk::simplify
