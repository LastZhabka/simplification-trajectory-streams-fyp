// The planar greedy stabbers of "Faster δ-simplification of Trajectory Streams", a Simplifier
// each. Unless said otherwise every region is the same regular m-gon moved to a vertex, and the
// starts are M = ⌈2π/ε⌉ points on the circle of radius δ. docs/convex-shadow-32.md §10–11.
#pragma once

#include <string>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

// Algorithm 1: each start keeps its own m offsets and the cone R(S[p], p), the whole plane or two
// rays. O(Mm) = O(ε^-3/2) per vertex.
class GreedyStabbingNaive : public Simplifier<2> {
 public:
  explicit GreedyStabbingNaive(const Params& params);

  [[nodiscard]] std::string name() const override { return "greedy-stabbing-naive"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

// Algorithm 1 as in the paper: the region is conv(G_v), the hull of the cells of side εδ/(2√2)
// that meet the disk of radius (1 + ε/2)δ, and the starts are every grid point of G_v, O(ε^-2).
class GreedyStabbingNaivePaper : public Simplifier<2> {
 public:
  explicit GreedyStabbingNaivePaper(const Params& params);

  [[nodiscard]] std::string name() const override { return "greedy-stabbing-naive-paper"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

// The same with the grid points of only the cells of G_v that meet the circle, O(ε^-1).
class GreedyStabbingNaivePaperBoundary : public Simplifier<2> {
 public:
  explicit GreedyStabbingNaivePaperBoundary(const Params& params);

  [[nodiscard]] std::string name() const override {
    return "greedy-stabbing-naive-paper-boundary";
  }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

// Algorithm 2: the m offsets are shared by every start, and the starts between two crossings of
// the circle with the offsets' lines share one polygon and its two corners. O(M + m²) = O(ε^-1)
// per vertex.
class GreedyStabbingShared : public Simplifier<2> {
 public:
  explicit GreedyStabbingShared(const Params& params);

  [[nodiscard]] std::string name() const override { return "greedy-stabbing-shared"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

}  // namespace ssk::simplify
