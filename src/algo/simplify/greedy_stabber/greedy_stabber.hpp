// The greedy of ordered stabbing, the same however a stabber keeps its candidates: advance them by
// each vertex, and when none survives, emit a segment from the ones that reached the vertex before
// and start again from this one. A stabber says what a candidate is and how it starts, advances
// and ends a segment. One is built per curve, as `init` fixes what depends on delta and eps.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
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

template <class Candidate>
class GreedyStabber {
 public:
  virtual ~GreedyStabber() = default;

  Curve<2> run(const Curve<2>& curve, double delta, double eps);

 protected:
  using Vec = geometry::Vec;

  virtual void init(double delta, double eps) = 0;
  virtual std::vector<Candidate> start(Vec v) = 0;
  // `stabbed` is the open run, the vertices the candidates have reached.
  virtual std::vector<Candidate> advance(const std::vector<Candidate>& cands,
                                         std::span<const Vec> stabbed, Vec v) = 0;
  // Called while `candidates_` still holds the ones that reached `stabbed.back()`.
  virtual void emit_segment(Curve<2>& out, std::span<const Vec> stabbed) = 0;

  static void emit(Curve<2>& out, Vec q) {
    const Point<2> pt{q.x, q.y};
    if (out.empty() || out.back() != pt) {
      out.push_back(pt);
    }
  }

  std::vector<Candidate> candidates_;
};

template <class Candidate>
Curve<2> GreedyStabber<Candidate>::run(const Curve<2>& curve, double delta, double eps) {
  // Ours reaches exactly (1 + eps) * delta, so aim inside by what the arithmetic can err. It runs
  // relative to the first vertex, erring by ulps of the curve's extent; only writing a point back
  // rounds by half an ulp of its own coordinates, at most scale + 2 * delta. Ulps of the
  // coordinates for all of it, 16 of them, put NGSIM's floor at 2.3e-8 ft where 1.4e-9 does.
  const Point<2> origin = curve[0];
  double scale = 0.0;
  double extent = 0.0;
  for (const Point<2>& q : curve) {
    scale = std::max({scale, std::abs(q[0]), std::abs(q[1])});
    extent = std::max({extent, std::abs(q[0] - origin[0]), std::abs(q[1] - origin[1])});
  }
  const double ulp = std::numeric_limits<double>::epsilon();
  const double inside = delta - ulp * (scale + 2.0 * delta) - 64.0 * ulp * extent;
  if (inside <= 0.0) {
    // Below what the coordinates resolve no new point can be placed within delta, but dropping
    // a repeat of the point before is exact. Returning the input whole kept every repeat and
    // broke the size bound on Mopsi's and NGSIM's runs of identical fixes.
    Curve<2> out;
    for (const Point<2>& q : curve) {
      if (out.empty() || out.back() != q) {
        out.push_back(q);
      }
    }
    if (out.size() == 1) {
      out.push_back(curve.back());
    }
    return out;
  }
  init(inside, eps);

  std::vector<Vec> vertices;
  for (const Point<2>& q : curve) {
    vertices.push_back({q[0] - origin[0], q[1] - origin[1]});
  }

  Curve<2> out;
  std::size_t first = 0;
  candidates_ = start(vertices[0]);
  for (std::size_t i = 1; i < vertices.size(); ++i) {
    const std::span<const Vec> stabbed = std::span(vertices).subspan(first, i - first);
    std::vector<Candidate> next = advance(candidates_, stabbed, vertices[i]);
    if (next.empty()) {
      emit_segment(out, stabbed);
      first = i;
      next = start(vertices[i]);
    }
    candidates_ = std::move(next);
  }
  if (first + 1 == vertices.size()) {
    emit(out, vertices.back());
  } else {
    emit_segment(out, std::span(vertices).subspan(first));
  }

  for (Point<2>& q : out) {
    q = {q[0] + origin[0], q[1] + origin[1]};
  }
  return out;
}

}  // namespace ssk::simplify
