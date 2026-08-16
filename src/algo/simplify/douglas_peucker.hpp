// Douglas-Peucker, ported from psimpl (psimpl.sourceforge.net), MIT.
//
// Faithful to psimpl's `simplify_douglas_peucker`, including the part people forget: it runs
// a radial-distance pass at the same tolerance first, and applies DP to that reduced
// polyline. Distances are point-to-*segment* (psimpl::math::segment_distance2), not
// point-to-line, and the recursion is an explicit stack, as upstream.
//
// psimpl works on a flat iterator range of interleaved coordinates and writes points; this
// works on Curve<D> and returns indices into the original curve, so the radial pass has to
// carry the surviving positions rather than copy coordinates.
#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

template <std::size_t D>
class DouglasPeucker : public SubsetSimplifier<D> {
 public:
  explicit DouglasPeucker(const Params& params) : tol_(params.get("tol")) {}

  [[nodiscard]] std::string name() const override { return "douglas-peucker"; }

 protected:
  [[nodiscard]] std::vector<std::size_t> run_indices(const Curve<D>& curve,
                                                     const Context&) const override {
    if (tol_ == 0.0) {
      std::vector<std::size_t> all(curve.size());
      for (std::size_t i = 0; i < all.size(); ++i) {
        all[i] = i;
      }
      return all;
    }

    const std::vector<std::size_t> reduced = radial_distance(curve);
    const std::size_t n = reduced.size();
    std::vector<bool> keys(n, false);
    keys.front() = true;
    keys.back() = true;

    const double tol2 = tol_ * tol_;
    std::vector<std::pair<std::size_t, std::size_t>> stack{{0, n - 1}};
    while (!stack.empty()) {
      const auto [first, last] = stack.back();
      stack.pop_back();
      const auto [index, dist2] = find_key(curve, reduced, first, last);
      if (index != 0 && tol2 < dist2) {
        keys[index] = true;
        stack.emplace_back(index, last);
        stack.emplace_back(first, index);
      }
    }

    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < n; ++i) {
      if (keys[i]) {
        out.push_back(reduced[i]);
      }
    }
    return out;
  }

 private:
  // psimpl runs this at the same tolerance before DP proper.
  [[nodiscard]] std::vector<std::size_t> radial_distance(const Curve<D>& curve) const {
    const double tol2 = tol_ * tol_;
    std::vector<std::size_t> out{0};
    std::size_t current = 0;
    for (std::size_t i = 1; i + 1 < curve.size(); ++i) {
      if (point_distance2(curve[current], curve[i]) < tol2) {
        continue;
      }
      current = i;
      out.push_back(i);
    }
    out.push_back(curve.size() - 1);
    return out;
  }

  // Furthest point from the segment, over the open range; index 0 means "none found".
  [[nodiscard]] static std::pair<std::size_t, double> find_key(
      const Curve<D>& curve, const std::vector<std::size_t>& reduced, std::size_t first,
      std::size_t last) {
    std::pair<std::size_t, double> key{0, 0.0};
    for (std::size_t i = first + 1; i < last; ++i) {
      const double d2 =
          segment_distance2(curve[reduced[first]], curve[reduced[last]], curve[reduced[i]]);
      if (d2 < key.second) {
        continue;
      }
      key = {i, d2};
    }
    return key;
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

  double tol_;
};

}  // namespace ssk::simplify
