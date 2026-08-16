// Douglas-Peucker, the point-count variant (DPn), ported from psimpl
// (psimpl.sourceforge.net) v7, MPL 1.1: `simplify_douglas_peucker_n` and
// `DPHelper::ApproximateN`.
//
// The knob is a vertex budget, not a tolerance. Instead of splitting one sub-polyline at a
// time until every point is within tol, all sub-polylines sit in a max-heap keyed by the
// squared distance of their furthest point, and the globally furthest point is promoted to a
// key at each step -- until `count` keys exist. The result is exactly `count` points.
//
// Upstream's tolerance form is *not* ported: with no tolerance there is no radial-distance
// pre-pass either, and DPn does not use one. Distances stay point-to-*segment*
// (psimpl::math::segment_distance2), as upstream.
//
// psimpl works on a flat iterator range of interleaved coordinates and writes points; this
// works on Curve<D> and returns indices into the curve, so upstream's coordinate indices
// (`index / DIM`) are point indices here.
#pragma once

#include <cstddef>
#include <queue>
#include <string>
#include <vector>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

template <std::size_t D>
class DouglasPeucker : public SubsetSimplifier<D> {
 public:
  explicit DouglasPeucker(const Params& params)
      : count_(static_cast<int>(params.get("count"))) {}

  [[nodiscard]] std::string name() const override { return "douglas-peucker-n"; }

 protected:
  [[nodiscard]] std::vector<std::size_t> run_indices(const Curve<D>& curve,
                                                     const Context&) const override {
    // Upstream copies the input through when it cannot deliver the budget, rather than
    // failing.
    if (count_ < 2 || curve.size() <= static_cast<std::size_t>(count_)) {
      std::vector<std::size_t> all(curve.size());
      for (std::size_t i = 0; i < all.size(); ++i) {
        all[i] = i;
      }
      return all;
    }

    const std::size_t n = curve.size();
    std::vector<bool> keys(n, false);
    keys.front() = true;
    keys.back() = true;

    if (count_ > 2) {
      std::priority_queue<SubPoly> queue;
      queue.push(with_key(curve, {0, n - 1}));

      int key_count = 2;
      while (!queue.empty()) {
        const SubPoly poly = queue.top();
        queue.pop();
        keys[poly.key_index] = true;
        if (++key_count == count_) {
          break;
        }
        const SubPoly left = with_key(curve, {poly.first, poly.key_index});
        if (left.key_index != 0) {
          queue.push(left);
        }
        const SubPoly right = with_key(curve, {poly.key_index, poly.last});
        if (right.key_index != 0) {
          queue.push(right);
        }
      }
    }

    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < n; ++i) {
      if (keys[i]) {
        out.push_back(i);
      }
    }
    return out;
  }

 private:
  struct SubPoly {
    std::size_t first = 0;
    std::size_t last = 0;
    std::size_t key_index = 0;  // 0 means "no key found"
    double key_dist2 = 0.0;

    bool operator<(const SubPoly& other) const { return key_dist2 < other.key_dist2; }
  };

  // Furthest point from the segment, over the open range; index 0 means "none found".
  [[nodiscard]] static SubPoly with_key(const Curve<D>& curve, SubPoly poly) {
    for (std::size_t i = poly.first + 1; i < poly.last; ++i) {
      const double d2 = segment_distance2(curve[poly.first], curve[poly.last], curve[i]);
      if (d2 < poly.key_dist2) {
        continue;
      }
      poly.key_index = i;
      poly.key_dist2 = d2;
    }
    return poly;
  }

  [[nodiscard]] static double point_distance2(const Point<D>& a, const Point<D>& b) {
    double sum = 0.0;
    for (std::size_t k = 0; k < D; ++k) {
      const double d = a[k] - b[k];
      sum += d * d;
    }
    return sum;
  }

  [[nodiscard]] static double segment_distance2(const Point<D>& s1, const Point<D>& s2,
                                                const Point<D>& p) {
    double cw = 0.0;
    double cv = 0.0;
    for (std::size_t k = 0; k < D; ++k) {
      const double v = s2[k] - s1[k];
      cw += (p[k] - s1[k]) * v;
      cv += v * v;
    }
    if (cw <= 0.0) {
      return point_distance2(p, s1);
    }
    if (cv <= cw) {
      return point_distance2(p, s2);
    }
    const double fraction = (cv == 0.0) ? 0.0 : cw / cv;
    Point<D> proj{};
    for (std::size_t k = 0; k < D; ++k) {
      proj[k] = s1[k] + (s2[k] - s1[k]) * fraction;
    }
    return point_distance2(p, proj);
  }

  int count_;
};

}  // namespace ssk::simplify
