// The baselines and ours at a fixed Frechet error, step 8 of docs/running.md.
#include "target.hpp"
#include "measure.hpp"
#include "simplify/douglas_peucker.hpp"
#include "test_main.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace {

using ssk::pipelines::Best;
using ssk::pipelines::Target;
using ssk::simplify::Curve;

struct Case {
  ssk::io::Document doc;
  Curve<2> curve;
  ssk::simplify::Context in;
};

Case wiggle() {
  Case c;
  c.doc.dim = 2;
  c.doc.name = "wiggle";
  c.doc.t_unit = "ms";
  for (int i = 0; i < 400; ++i) {
    const double t = 0.03 * i;
    const double x = 10.0 * t + std::sin(7.0 * t);
    const double y = 3.0 * std::sin(t) + 0.2 * std::sin(41.0 * t);
    c.doc.points.push_back({x, y});
    c.doc.t.push_back(1000.0 * i);
    c.curve.push_back({x, y});
  }
  c.in = ssk::simplify::context_of(c.doc);
  return c;
}

Curve<2> kept(const Curve<2>& curve, const std::vector<std::size_t>& indices) {
  Curve<2> out;
  for (const std::size_t i : indices) {
    out.push_back(curve[i]);
  }
  return out;
}

TEST("target/levels are fractions of the diagonal, ascending, none without extent") {
  const auto ts = ssk::pipelines::targets(200.0);
  CHECK(ts.size() == ssk::pipelines::kFractions.size());
  CHECK(ts.front().level == "f0.0003");
  CHECK(ts.back().level == "f0.03");
  CHECK_CLOSE(ts[1].distance, 0.2, 1e-12);
  for (std::size_t i = 1; i < ts.size(); ++i) {
    CHECK(ts[i].distance > ts[i - 1].distance);
  }
  CHECK(ssk::pipelines::targets(0.0).empty());
}

TEST("target/every baseline is within each target, and fewer vertices as it loosens") {
  const Case c = wiggle();
  const ssk::frechet::Frechet<2> computer;
  const auto ts = ssk::pipelines::targets(ssk::pipelines::diagonal(c.curve));
  for (const std::string& algorithm : ssk::pipelines::kTargetBaselines) {
    const std::vector<Best> bests = ssk::pipelines::search(algorithm, computer, c.curve, c.in, ts);
    CHECK(bests.size() == ts.size());
    for (std::size_t i = 0; i < bests.size(); ++i) {
      CHECK_MSG(computer.within(c.curve, kept(c.curve, bests[i].indices), ts[i].distance),
                algorithm + " at " + ts[i].level);
      CHECK_MSG(bests[i].input == (bests[i].indices.size() == c.curve.size()), algorithm);
      if (i > 0) {
        CHECK_MSG(bests[i].indices.size() <= bests[i - 1].indices.size(),
                  algorithm + " at " + ts[i].level);
      }
    }
    // At the tightest level SQUISH-E keeps everything here: dropping even one point of the
    // 41-cycle wiggle costs 2.3 times the target.
    CHECK_MSG(bests.back().indices.size() < c.curve.size() / 4, algorithm + " compresses");
  }
}

TEST("target/the budget found is a boundary: one vertex fewer is not within") {
  const Case c = wiggle();
  const ssk::frechet::Frechet<2> computer;
  const auto ts = ssk::pipelines::targets(ssk::pipelines::diagonal(c.curve));
  for (const Best& best : ssk::pipelines::search("douglas-peucker-n", computer, c.curve, c.in, ts)) {
    const std::size_t k = best.indices.size();
    if (k <= 2) {
      continue;
    }
    const auto fewer = ssk::simplify::DouglasPeucker<2>({{"count", static_cast<double>(k - 1)}})
                           .indices(c.curve);
    CHECK_MSG(!computer.within(c.curve, kept(c.curve, fewer), best.target.distance),
              best.target.level);
  }
}

TEST("target/ours runs at T / (1 + eps) and is within T") {
  const Case c = wiggle();
  const ssk::frechet::Frechet<2> computer;
  for (const Target& t : ssk::pipelines::targets(ssk::pipelines::diagonal(c.curve))) {
    const auto run =
        ssk::pipelines::run_ours(computer, "greedy-stabbing-shared", c.curve, t, 0.25, false);
    CHECK_CLOSE(run.delta, t.distance / 1.25, 1e-15 * t.distance);
    CHECK(run.within && run.within_target);
    CHECK(run.measured.distance <= t.distance * (1.0 + 1e-9));
  }
}

TEST("target/unscaled, ours runs at T, within (1 + eps) T and no larger than scaled") {
  const Case c = wiggle();
  const ssk::frechet::Frechet<2> computer;
  for (const Target& t : ssk::pipelines::targets(ssk::pipelines::diagonal(c.curve))) {
    const auto scaled =
        ssk::pipelines::run_ours(computer, "greedy-stabbing-shared", c.curve, t, 0.5, false);
    const auto unscaled =
        ssk::pipelines::run_ours(computer, "greedy-stabbing-shared", c.curve, t, 0.5, true);
    CHECK(unscaled.delta == t.distance);
    CHECK(unscaled.within);
    CHECK(unscaled.measured.distance <= 1.5 * t.distance * (1.0 + 1e-9));
    CHECK_MSG(unscaled.simplified.size() <= scaled.simplified.size(), t.level);
  }
}

TEST("target/the match summary pairs ours with every baseline at its level") {
  const Case c = wiggle();
  const ssk::frechet::Frechet<2> computer;
  const auto ts = ssk::pipelines::targets(ssk::pipelines::diagonal(c.curve));
  std::vector<ssk::pipelines::Found> found;
  for (const std::string& algorithm : {std::string("douglas-peucker-n"), std::string("squish")}) {
    for (Best& best : ssk::pipelines::search(algorithm, computer, c.curve, c.in, ts)) {
      const auto measured = ssk::pipelines::measure(computer, c.curve, kept(c.curve, best.indices), 1e-9);
      found.push_back({algorithm, std::move(best), measured});
    }
  }
  const auto searched = ssk::pipelines::search_summary_to_json(c.doc, "test", "test/w.json", ts, found);
  const auto back = ssk::pipelines::targets_of(searched);
  CHECK(back.size() == ts.size());
  CHECK(back[2].level == ts[2].level && back[2].distance == ts[2].distance);

  std::vector<ssk::pipelines::OursRun> runs;
  for (const Target& t : back) {
    runs.push_back(
        ssk::pipelines::run_ours(computer, "greedy-stabbing-shared", c.curve, t, 0.5, false));
  }
  const auto matched = ssk::pipelines::match_summary_to_json(c.doc, "test", "test/w.json",
                                                             "greedy-stabbing-shared", 0.5, false,
                                                             searched, runs);
  CHECK(matched.find("delta_rule")->str() == "T/(1+eps)");
  const auto& pairs = matched.find("pairs")->arr();
  CHECK(pairs.size() == 2 * ts.size());
  for (const auto& p : pairs) {
    CHECK(p.find("ours")->find("within_bound")->boolean());
    CHECK(p.find("baseline")->find("distance")->num() <= p.find("target")->num() * (1 + 1e-9));
  }
}

}  // namespace
