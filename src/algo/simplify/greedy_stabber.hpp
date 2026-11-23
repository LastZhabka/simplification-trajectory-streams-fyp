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
#include <utility>
#include <vector>

#include "geometry/geometry.hpp"
#include "simplify/simplifier.hpp"

namespace ssk::simplify {

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
  // Ours reaches exactly (1 + eps) * delta, and writing a point back in the input's coordinates
  // can round it an ulp of them further, so aim a few ulps inside.
  double scale = 0.0;
  for (const Point<2>& q : curve) {
    scale = std::max({scale, std::abs(q[0]), std::abs(q[1])});
  }
  const double inside = delta - 16.0 * std::numeric_limits<double>::epsilon() * scale;
  if (inside <= 0.0) {
    return curve;
  }
  init(inside, eps);

  // Relative to the first vertex, so that coordinates near zero keep more of delta's digits.
  const Point<2> origin = curve[0];
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
