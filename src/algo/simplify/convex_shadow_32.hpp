// Streaming delta-simplification, Algorithm 1 (Simplify) of Cheng, Huang and Jiang,
// "Simplification of Trajectory Streams", SoCG 2025 (arXiv:2503.23025). Polygon code from
// archive/streaming-simplification.
//
// G_v is the grid of side eps * delta / (2 sqrt 2), with v as a grid point, cut to the cells
// meeting the ball of radius (1 + eps/2) * delta around v. Each candidate start p keeps a
// convex shadow S[p], the points where a segment from p can be after stabbing conv(G_v) of
// every vertex so far in order. When every shadow is empty the last segment is emitted and a
// new one starts at the vertex that emptied them. The output is within (1 + eps) * delta and
// has at most 2 * opt - 2 vertices (Theorem 3).
//
// One deviation: only the grid points of G_v's boundary cells start a segment, O(1/eps) of
// them rather than O(1/eps^2). Lemma 2's bound still holds; see boundary_points.
//
// ConvexShadow32 is ConvexShadow with one change: the region around each vertex v is not
// conv(G_v) but P_v, the regular m-gon around v with circumradius (1 + eps) * delta and the
// fewest sides whose inradius reaches (1 + eps/2) * delta, the two radii conv(G_v) lies
// between. That is m = O(eps^-1/2) sides against conv(G_v)'s O(eps^-2/3) corners, and no
// convex polygon between those radii has fewer: docs/convex-shadow-32.md. The starts are still
// the boundary grid points, and both guarantees hold unchanged.
//
// The output points are grid and shadow points, not input vertices, so this is a Simplifier
// and not a SubsetSimplifier.
#pragma once

#include <string>

#include "simplify/simplifier.hpp"
#include "validate/validate.hpp"

namespace ssk::simplify {

class ConvexShadow32 : public Simplifier<2> {
 public:
  // With a validator, every grid, shadow and emitted segment is checked as it is made and
  // the outcome counted there. It costs time, so it is off unless one is given.
  explicit ConvexShadow32(const Params& params, validate::Validator* validator = nullptr);

  [[nodiscard]] std::string name() const override { return "convex-shadow-32"; }

  // m: the fewest sides with cos(pi/m) >= (1 + eps/2) / (1 + eps). It depends on eps alone.
  [[nodiscard]] static int sides(double eps);

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
  validate::Validator* validator_;
};

}  // namespace ssk::simplify
