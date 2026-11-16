// Our algorithm at the baselines' operating points -- step 3 of docs/comparison.md.
//
// ssk_frechet measured each baseline simplification at some Frechet distance d*. Ours runs at
// delta = d* / (1 + eps): its guarantee is (1 + eps) * delta, so it is within d* by
// construction, and what is left to compare is the protocol's question -- at that error, how
// many vertices. eps has no baseline value to match, so it is swept, one run per value.
//
// Ours is measured in the same pass, with the computer and tolerance ssk_frechet used on the
// baselines: the protocol reports measured error on both sides, and measuring here saves
// reading every result back.
//
// Ours has more than one implementation. Each is run by the name it gives itself, and its
// results are kept apart from the others', so two implementations can be compared on the same
// baseline operating points.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "io/json.hpp"
#include "io/trajectory.hpp"
#include "measure.hpp"
#include "simplify/simplifier.hpp"
#include "validate/validate.hpp"

namespace ssk::pipelines {

// The first is the default.
inline const std::vector<std::string> kOurs{"convex-shadow", "convex-shadow-far-points",
                                            "convex-shadow-32"};

struct Match {
  std::string ours;      // which implementation, one of kOurs
  json::Value baseline;  // its frechet-distances document
  double delta = 0.0;
  simplify::Curve<2> simplified;
  Measurement measured;
  bool within_bound = false;  // the exact decision at the guarantee, (1 + eps) * delta
  std::optional<validate::Stats> validation;  // set when the run was validated
};

[[nodiscard]] Match match(const frechet::Frechet<2>& computer, const std::string& ours,
                          const simplify::Curve<2>& input, json::Value baseline, double eps,
                          double relative_tol, bool validate);

// {"checked", "failed", "checks": {<name>: {"checked", "failed"}}}
[[nodiscard]] json::Value validation_to_json(const validate::Stats& stats);

// Ours as a result document: a trajectory document of its points -- no clock, since its
// vertices are not input vertices -- with params, stats, its measurement and the baseline
// document it was matched to.
[[nodiscard]] json::Value match_to_json(const io::Document& input, const std::string& dataset,
                                        const std::string& source, double eps,
                                        const Match& match);

// Every match of one trajectory as one small document, so scripts/comparison_stats.py opens
// one file per trajectory rather than every result. Validation counts are summed over them.
[[nodiscard]] json::Value summary_to_json(const io::Document& input, const std::string& dataset,
                                          const std::string& source, const std::string& ours,
                                          double eps, const std::vector<Match>& matches);

}  // namespace ssk::pipelines
