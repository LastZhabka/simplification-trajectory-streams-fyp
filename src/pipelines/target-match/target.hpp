// The baselines and ours at a fixed Frechet error -- step 8 of docs/running.md, the protocol
// "fix error, compare compression" of docs/comparison.md.
//
// Each trajectory gets targets T = f * its bounding-box diagonal, for f in kFractions, so a
// level means the same in degrees, feet and pixels. A baseline at a target is the fewest
// vertices it was found to reach within T: a bisection on its own knob, deciding each step
// exactly with Frechet::within. Ours at a target runs at delta = T / (1 + eps), within T by
// its guarantee, and does not depend on which baseline it is compared with.
#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "frechet/frechet.hpp"
#include "io/json.hpp"
#include "io/trajectory.hpp"
#include "measure.hpp"
#include "simplify/simplifier.hpp"

namespace ssk::pipelines {

// Half a decade apart, around the medians of step 6's d*/diagonal at m2..m6 (8e-4 to 2.5e-2).
inline const std::vector<double> kFractions{3e-4, 1e-3, 3e-3, 1e-2, 3e-2};

inline const std::vector<std::string> kTargetBaselines{"douglas-peucker-n", "squish", "squish-e",
                                                       "dots"};

struct Target {
  std::string level;  // "f0.001", the directory its documents go under
  double fraction = 0.0;
  double diagonal = 0.0;
  double distance = 0.0;  // fraction * diagonal
};

// Ascending; none for a curve with no extent, where every target would be 0.
[[nodiscard]] std::vector<Target> targets(double diagonal);

struct Best {
  Target target;
  std::vector<std::pair<std::string, double>> params;  // the knob it was run at
  std::vector<std::size_t> indices;
  int algorithm_runs = 0;  // runs of the algorithm the search for this target made
  bool input = false;      // no setting tried was within the target: the input itself
};

// One per target, ascending. Bisection assumes a smaller output is never closer; where that
// fails it stops at a boundary, so the size is the fewest found, not proven the fewest.
[[nodiscard]] std::vector<Best> search(const std::string& algorithm,
                                       const frechet::Frechet<2>& computer,
                                       const simplify::Curve<2>& curve,
                                       const simplify::Context& in,
                                       const std::vector<Target>& targets);

[[nodiscard]] json::Value best_to_json(const io::Document& doc, const std::string& dataset,
                                       const std::string& source, const std::string& algorithm,
                                       const Best& best, const Measurement& measured);

struct Found {
  std::string algorithm;
  Best best;
  Measurement measured;
};

// Every baseline at every target of one trajectory, the targets with them. Written last, so with
// --skip-existing a trajectory that has one is done.
[[nodiscard]] json::Value search_summary_to_json(const io::Document& doc,
                                                 const std::string& dataset,
                                                 const std::string& source,
                                                 const std::vector<Target>& targets,
                                                 const std::vector<Found>& found);

[[nodiscard]] std::vector<Target> targets_of(const json::Value& search_summary);

struct OursRun {
  Target target;
  double delta = 0.0;
  simplify::Curve<2> simplified;
  Measurement measured;
  bool within = false;         // the exact decision at ours' guarantee, (1 + eps) * delta
  bool within_target = false;  // and at the target, the same unless `unscaled`
};

// At delta = T / (1 + eps), within T by the guarantee; with `unscaled` at delta = T, within
// (1 + eps) T and its size bound against kappa(T) -- step 10 of docs/running.md.
[[nodiscard]] OursRun run_ours(const frechet::Frechet<2>& computer, const std::string& ours,
                               const simplify::Curve<2>& curve, const Target& target, double eps,
                               bool unscaled);

[[nodiscard]] json::Value ours_to_json(const io::Document& doc, const std::string& dataset,
                                       const std::string& source, const std::string& ours,
                                       double eps, bool unscaled, const OursRun& run);

// Ours beside every baseline at each target, in the shape of ssk_match's summary with the level
// where it has m: {"algorithm", "level", "fraction", "target", "baseline": {...}, "ours": {...}}.
[[nodiscard]] json::Value match_summary_to_json(const io::Document& doc,
                                                const std::string& dataset,
                                                const std::string& source,
                                                const std::string& ours, double eps,
                                                bool unscaled, const json::Value& search_summary,
                                                const std::vector<OursRun>& runs);

}  // namespace ssk::pipelines
