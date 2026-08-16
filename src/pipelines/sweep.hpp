// Driving the baselines to a fixed compression rate.
//
// For a trajectory of N points, rate M asks for ceil(N / 2^M) vertices. DPn and SQUISH take
// that budget as their parameter and hit it exactly; DOTS bounds error rather than size, so
// its threshold has to be searched for -- see sweep_dots. The protocol and what it is for are
// in docs/comparison.md.
//
// A Sweep is one algorithm's whole set of rates for one trajectory, which is what lands in
// one output document. Indices, not points: all three algorithms return a subsequence, so
// indices are lossless and a third of the size.
#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "io/json.hpp"
#include "io/trajectory.hpp"
#include "simplify/simplifier.hpp"

namespace ssk::pipelines {

using simplify::Context;
using simplify::Curve;

struct Run {
  int m = 0;                  // compression rate is 2^m
  std::size_t target = 0;
  std::vector<std::pair<std::string, double>> params;
  std::vector<std::size_t> indices;
  int algorithm_runs = 1;     // more than one only where the budget was searched for
};

struct Sweep {
  std::string algorithm;
  std::vector<Run> runs;
};

// ceil(N / 2^m) for m = 1..max_m, stopping at the first budget below `floor` -- a budget
// under 2 cannot be met when both endpoints are always kept, and SQUISH needs more than 4.
std::vector<std::size_t> budget_targets(std::size_t n, int max_m, std::size_t floor);

Sweep sweep_dpn(const Curve<2>& curve, int max_m);
Sweep sweep_squish(const Curve<2>& curve, const Context& in, int max_m);
Sweep sweep_dots(const Curve<2>& curve, const Context& in, int max_m);

// The output document: the sweep, plus enough provenance to rehydrate points from `source`.
json::Value to_json(const io::Document& doc, const std::string& source, const Sweep& sweep);

}  // namespace ssk::pipelines
