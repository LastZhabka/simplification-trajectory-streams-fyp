// Ours at each baseline's size -- step 11 of docs/running.md, the protocol "fix size, compare
// error" of docs/comparison.md.
//
// ssk_frechet measured each baseline at its budget, N / 2^m vertices or a point or two under, at
// some Frechet distance d*. Ours cannot be set by size, so a bisection on delta finds where it
// fits in the baseline's own size: the top of a bracket [lo, hi], ours at hi within the size and
// at lo not, to the tolerance ssk_frechet measured with. Ours at hi is then measured the same
// way, and the comparison is of error alone: ours' measured distance against d*.
//
// The search depends on the trajectory and the size, not on the baseline, so every size of one
// trajectory shares one cache of runs: two baselines at the same size get the same run of ours.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "frechet/frechet.hpp"
#include "io/json.hpp"
#include "io/trajectory.hpp"
#include "measure.hpp"
#include "simplify/simplifier.hpp"

namespace ssk::pipelines {

struct Fit {
  std::size_t budget = 0;  // the baseline's vertices
  double delta = 0.0;      // ours at delta has at most `budget` vertices
  double below = 0.0;      // the bracket's bottom, a tolerance under delta: ours there has more
  simplify::Curve<2> simplified;
  Measurement measured;
  bool within_bound = false;  // the exact decision at the guarantee, (1 + eps) * delta
  int runs = 0;               // runs of ours its search made, beyond those already cached
};

// The vertices of a baseline's frechet-distances document: what ours is fitted to.
[[nodiscard]] std::size_t budget_of(const json::Value& baseline);

// One per budget, in their order. Bisection assumes a larger delta never gives more vertices;
// where that fails it stops at a boundary, so delta is the top of a bracket, not proven the least.
[[nodiscard]] std::vector<Fit> fit(const frechet::Frechet<2>& computer, const std::string& ours,
                                   const simplify::Curve<2>& curve, double eps,
                                   const std::vector<std::size_t>& budgets);

// Ours as a result document, as ssk_match writes one, with the baseline document it was fitted to.
[[nodiscard]] json::Value fit_to_json(const io::Document& input, const std::string& dataset,
                                      const std::string& source, const std::string& ours,
                                      double eps, const Fit& fit, const json::Value& baseline);

struct Fitted {
  json::Value baseline;  // its frechet-distances document
  Fit fit;
};

// Every baseline beside ours at its size, in the shape of ssk_match's summary, so
// scripts/ratio_buckets.py reads both the same way; its delta_rule is "size".
[[nodiscard]] json::Value fit_summary_to_json(const io::Document& input,
                                              const std::string& dataset,
                                              const std::string& source, const std::string& ours,
                                              double eps, const std::vector<Fitted>& fitted);

}  // namespace ssk::pipelines
