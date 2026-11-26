// Algorithm 1 as Cheng, Huang and Jiang write it, a Simplifier each: every start keeps its shadow
// S[p] itself, a polygon, and steps it as S[p] <- A_v ∩ F(S[p], p), with no offsets or cone kept
// beside it, intersecting the half-planes as their Lemma 1(ii) does, by a convex hull in the dual.
// Each has the region A_v and the starts of the decomposed stabber named with naive in place of
// no-decomp, so only the step differs. docs/convex-shadow-32.md §12.
#pragma once

#include <string>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

// The regular m-gon and M = ⌈2π/ε⌉ starts on the circle of radius δ, as GreedyStabbingNaive.
class GreedyStabbingNoDecomp : public Simplifier<2> {
 public:
  explicit GreedyStabbingNoDecomp(const Params& params);

  [[nodiscard]] std::string name() const override { return "greedy-stabbing-no-decomp"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

// The paper's own: conv(G_v), and every grid point of G_v, as GreedyStabbingNaivePaper.
class GreedyStabbingNoDecompPaper : public Simplifier<2> {
 public:
  explicit GreedyStabbingNoDecompPaper(const Params& params);

  [[nodiscard]] std::string name() const override { return "greedy-stabbing-no-decomp-paper"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

// conv(G_v), and the grid points of only the cells of G_v that meet the circle, as
// GreedyStabbingNaivePaperBoundary.
class GreedyStabbingNoDecompPaperBoundary : public Simplifier<2> {
 public:
  explicit GreedyStabbingNoDecompPaperBoundary(const Params& params);

  [[nodiscard]] std::string name() const override {
    return "greedy-stabbing-no-decomp-paper-boundary";
  }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
};

}  // namespace ssk::simplify
