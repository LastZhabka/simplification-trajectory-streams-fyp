// Continuous Frechet distance, ported from the Dutsch-Vahrenhold implementation in
// mwernerds/frechetrange (MIT) -- `detail/dv/frechet_distance.hpp` and `distance_sqr.hpp`.
// That code won the ACM SIGSPATIAL GIS Cup 2017, whose task was exactly this decision
// problem, and it needs nothing outside the standard library, so it ports cleanly here.
//
// The primitive is the decision `within(a, b, r)`: is d_F(a, b) <= r? It is exact, and it is
// what the algorithm actually computes. The value comes from bisecting it, so it costs a
// tolerance and roughly nine times a single decision -- see docs/frechet-distance.md.
//
// The frontline buffer is a member rather than a local so a driver measuring a whole corpus
// reuses one allocation across millions of calls, as upstream intends.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace ssk::frechet {

template <std::size_t D>
using Point = std::array<double, D>;

template <std::size_t D>
using Curve = std::vector<Point<D>>;

template <std::size_t D>
class Frechet {
 public:
  [[nodiscard]] bool within(const Curve<D>& a, const Curve<D>& b, double r) const {
    if (a.empty() || b.empty()) {
      return a.empty() && b.empty();
    }
    const double r2 = r * r;

    if (dist2(a.front(), b.front()) > r2 || dist2(a.back(), b.back()) > r2) {
      return false;
    }

    const bool a_smaller = a.size() <= b.size();
    const Curve<D>& small = a_smaller ? a : b;
    const Curve<D>& big = a_smaller ? b : a;

    if (small.size() == 1) {
      return point_covers(small[0], big, r2);
    }
    if (!monotone_match(a, b, r2) || !monotone_match(b, a, r2)) {
      return false;
    }
    return traverse(small, big, r2);
  }

  // Bisection on `within`. Upstream provides no value routine; this one is ours.
  [[nodiscard]] double distance(const Curve<D>& a, const Curve<D>& b, double tol) const {
    if (a.empty() || b.empty()) {
      return 0.0;
    }
    double hi = std::max(std::sqrt(dist2(a.front(), b.front())),
                         std::sqrt(dist2(a.back(), b.back())));
    if (hi <= 0.0) {
      hi = tol;
    }
    while (!within(a, b, hi)) {
      hi *= 2.0;
    }
    double lo = 0.0;
    while (hi - lo > tol) {
      const double mid = 0.5 * (lo + hi);
      if (within(a, b, mid)) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    return hi;
  }

 private:
  // Upstream's sentinel: a segment is reachable when its beginning is <= 1.
  static constexpr double kUnreachable = 2.0;

  mutable std::vector<double> left_begins_;

  // Coordinates are summed from D-1 down to 0, as upstream's template recursion does, so
  // the floating-point result is identical.
  [[nodiscard]] static double dist2(const Point<D>& p, const Point<D>& q) {
    double sum = 0.0;
    for (std::size_t k = D; k-- > 0;) {
      const double d = p[k] - q[k];
      sum += d * d;
    }
    return sum;
  }

  [[nodiscard]] static bool reachable(double begin) { return begin <= 1.0; }

  [[nodiscard]] static bool point_covers(const Point<D>& p, const Curve<D>& c, double r2) {
    for (std::size_t i = 1; i < c.size(); ++i) {
      if (dist2(p, c[i]) > r2) {
        return false;
      }
    }
    return true;
  }

  // Where the line through p1,p2 meets the ball of radius^2 around cp, as the scalars s with
  // p1 + s*(p2-p1) on the boundary. `begin` becomes kUnreachable when they do not meet.
  static void intersections(const Point<D>& p1, const Point<D>& p2, const Point<D>& cp,
                            double r2, double& begin, double& end) {
    double a = 0.0;
    double b = 0.0;
    double c = -r2;
    for (std::size_t k = D; k-- > 0;) {
      const double u = p2[k] - p1[k];
      const double v = p1[k] - cp[k];
      a += u * u;
      b += u * v;
      c += v * v;
    }
    b *= 2.0;

    const double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
      begin = kUnreachable;
      return;
    }
    if (a == 0.0) {
      if (dist2(p1, cp) <= r2) {
        begin = 0.0;
        end = 1.0;
      } else {
        begin = end = kUnreachable;
      }
      return;
    }
    const double root = std::sqrt(discriminant);
    begin = (-b - root) / (2.0 * a);
    end = (-b + root) / (2.0 * a);
  }

  // Every interior point of `points` must match a point on `segments`, monotonically. A cheap
  // one-dimensional necessary condition, run both ways before the diagram is built.
  [[nodiscard]] static bool monotone_match(const Curve<D>& points, const Curve<D>& segments,
                                           double r2) {
    const std::size_t last_point = points.size() - 1;
    const std::size_t segment_count = segments.size() - 1;
    if (last_point <= 1 || segment_count == 0) {
      return true;
    }

    std::size_t point_idx = 1;
    std::size_t seg_idx = 0;
    double consumed = 0.0;
    double begin = 0.0;
    double end = 1.0;
    while (true) {
      intersections(segments[seg_idx], segments[seg_idx + 1], points[point_idx], r2, begin,
                    end);
      if (begin <= 1.0 && end >= consumed) {
        consumed = std::max(consumed, begin);
        if (++point_idx == last_point) {
          return true;
        }
      } else {
        consumed = 0.0;
        if (++seg_idx == segment_count) {
          return false;
        }
      }
    }
  }

  [[nodiscard]] static double reachable_begin(double begin, double end, double parallel,
                                              double orthogonal) {
    if (end < 0.0) {
      return kUnreachable;
    }
    if (begin > 1.0) {
      return begin;
    }
    if (reachable(orthogonal)) {
      return begin;
    }
    if (begin < parallel) {
      return (end >= parallel) ? parallel : kUnreachable;
    }
    return begin;
  }

  // The free-space diagram, one column at a time, keeping only the frontline: the reachable
  // beginnings of the previous column's right segments.
  [[nodiscard]] bool traverse(const Curve<D>& p1, const Curve<D>& p2, double r2) const {
    const std::size_t rows = p1.size() - 1;
    if (left_begins_.size() < rows) {
      left_begins_.resize(rows, kUnreachable);
    }
    std::fill_n(left_begins_.begin(), rows, kUnreachable);
    left_begins_[0] = 0.0;

    double bottom_begin = 0.0;
    std::size_t bottom_row = 0;
    const std::size_t last_row = p1.size() - 2;
    const std::size_t last_column = p2.size() - 2;

    for (std::size_t col = 0; col < last_column; ++col) {
      if (bottom_begin == 0.0 && dist2(p1[0], p2[col]) > r2) {
        bottom_begin = kUnreachable;
      }

      double cur_bottom_begin = bottom_begin;
      double cur_bottom_end = cur_bottom_begin;
      double right_begin = kUnreachable;
      double right_end = kUnreachable;

      for (std::size_t row = bottom_row; row < last_row; ++row) {
        if (reachable(left_begins_[row]) || reachable(cur_bottom_begin)) {
          const bool top_right_free = dist2(p1[row + 1], p2[col + 1]) <= r2;

          if (top_right_free && cur_bottom_end >= 1.0 && cur_bottom_begin <= 1.0) {
            right_begin = 0.0;
            right_end = 1.0;
          } else {
            intersections(p1[row], p1[row + 1], p2[col + 1], r2, right_begin, right_end);
            right_begin =
                reachable_begin(right_begin, right_end, left_begins_[row], cur_bottom_begin);
          }

          double top_begin = 0.0;
          double top_end = 1.0;
          if (!(top_right_free && left_begins_[row + 1] <= 0.0)) {
            intersections(p2[col], p2[col + 1], p1[row + 1], r2, top_begin, top_end);
            top_begin =
                reachable_begin(top_begin, top_end, cur_bottom_begin, left_begins_[row]);
          }

          left_begins_[row] = right_begin;
          cur_bottom_begin = top_begin;
          cur_bottom_end = top_end;
        }

        if (bottom_row == row && !reachable(left_begins_[row])) {
          ++bottom_row;
        }
      }

      if (reachable(left_begins_[last_row]) || reachable(cur_bottom_begin)) {
        if (dist2(p1[last_row + 1], p2[col + 1]) <= r2 && cur_bottom_end >= 1.0 &&
            cur_bottom_begin <= 1.0) {
          right_begin = 0.0;
        } else {
          intersections(p1[last_row], p1[last_row + 1], p2[col + 1], r2, right_begin,
                        right_end);
          right_begin = reachable_begin(right_begin, right_end, left_begins_[last_row],
                                        cur_bottom_begin);
        }
        left_begins_[last_row] = right_begin;
      }

      if (bottom_row == last_row && !reachable(left_begins_[last_row])) {
        return false;
      }
    }

    if (reachable(left_begins_[last_row])) {
      return true;
    }

    double bottom = kUnreachable;
    for (std::size_t row = bottom_row; row < last_row; ++row) {
      if (reachable(left_begins_[row]) || reachable(bottom)) {
        if (dist2(p1[row + 1], p2[last_column + 1]) <= r2 && left_begins_[row + 1] <= 0.0) {
          bottom = 0.0;
        } else {
          double top_begin = 0.0;
          double top_end = 1.0;
          intersections(p2[last_column], p2[last_column + 1], p1[row + 1], r2, top_begin,
                        top_end);
          bottom = reachable_begin(top_begin, top_end, bottom, left_begins_[row]);
        }
      }
    }
    return reachable(bottom);
  }
};

}  // namespace ssk::frechet
