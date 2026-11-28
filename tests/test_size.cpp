// ssk_size_match, step 9: ours bisected on delta to each baseline's size, its error measured.
#include "size.hpp"
#include "match.hpp"
#include "measure.hpp"
#include "simplify/douglas_peucker.hpp"
#include "test_main.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ssk::json::Array;
using ssk::json::Object;
using ssk::json::Value;
using ssk::pipelines::Fit;
using ssk::simplify::Curve;

const std::string kShared = "greedy-stabbing-shared";

Curve<2> wander(std::size_t n) {
  Curve<2> c;
  double x = 0.0;
  double y = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({x, y});
    x += 0.5 + 0.3 * std::sin(0.7 * static_cast<double>(i));
    y += std::cos(0.45 * static_cast<double>(i));
  }
  return c;
}

ssk::io::Document document(const Curve<2>& curve) {
  ssk::io::Document doc;
  doc.dim = 2;
  doc.name = "mopsi/7/1400";
  for (const auto& p : curve) {
    doc.points.push_back({p[0], p[1]});
  }
  return doc;
}

std::size_t size_at(const std::string& ours, const Curve<2>& curve, double delta, double eps) {
  return ssk::pipelines::make_ours(ours, {{"delta", delta}, {"eps", eps}})->simplify(curve).size();
}

// DPn at N / 2^m, measured, as ssk_frechet writes it.
Value baseline(const ssk::frechet::Frechet<2>& computer, const Curve<2>& curve, int m) {
  const auto k = static_cast<std::size_t>(std::ceil(curve.size() / std::pow(2.0, m)));
  const Curve<2> kept =
      ssk::simplify::DouglasPeucker<2>({{"count", static_cast<double>(k)}}).simplify(curve);
  const ssk::pipelines::Measurement measured = ssk::pipelines::measure(computer, curve, kept, 1e-9);
  return Value(Object{
      {"algorithm", Value("douglas-peucker-n")},
      {"mode", Value("budget")},
      {"params", Value(Object{{"count", Value(k)}})},
      {"stats", Value(Object{{"m", Value(m)},
                             {"output_size", Value(kept.size())},
                             {"target", Value(k)},
                             {"algorithm_runs", Value(1)}})},
      {"frechet", Value(Object{{"distance", Value(measured.distance)},
                               {"tolerance", Value(measured.tolerance)},
                               {"computer", Value("dv-gis-cup-2017")}})},
      {"simplified", Value("data/simplified-trajectories/douglas-peucker-n/mopsi/m2/x.json")},
      {"dataset", Value("mopsi")},
      {"source", Value("mopsi/x.json")},
  });
}

TEST("size/ours fits in each budget at delta, and not a tolerance under it") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(200);
  const double tol = 1e-9 * ssk::pipelines::diagonal(curve);
  const std::vector<std::size_t> budgets{2, 3, 7, 13, 25, 50, 100};
  for (const std::string& ours : {kShared, std::string("greedy-stabbing-naive")}) {
    for (const double eps : {0.1, 0.5}) {
      const std::vector<Fit> fits = ssk::pipelines::fit(computer, ours, curve, eps, budgets);
      CHECK(fits.size() == budgets.size());
      for (std::size_t i = 0; i < fits.size(); ++i) {
        const Fit& f = fits[i];
        const std::string where = ours + ", budget " + std::to_string(budgets[i]);
        CHECK_MSG(f.budget == budgets[i], where);
        CHECK_MSG(f.simplified.size() <= f.budget, where + ": over the budget");
        CHECK_MSG(f.simplified ==
                      ssk::pipelines::make_ours(ours, {{"delta", f.delta}, {"eps", eps}})
                          ->simplify(curve),
                  where + ": not ours at delta");
        CHECK_MSG(f.below < f.delta && f.delta - f.below <= tol, where + ": bracket too wide");
        CHECK_MSG(f.below == 0.0 || size_at(ours, curve, f.below, eps) > f.budget,
                  where + ": fits at the bracket's bottom");
      }
      CHECK(fits.front().runs > 0);
    }
  }
}

TEST("size/ours is measured as the baselines are, and within its guarantee") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(150);
  for (const Fit& f : ssk::pipelines::fit(computer, kShared, curve, 0.1, {2, 10, 38})) {
    CHECK(f.within_bound);
    CHECK(f.within_bound == computer.within(curve, f.simplified, 1.1 * f.delta));
    CHECK_CLOSE(f.measured.tolerance, 1e-9 * ssk::pipelines::diagonal(curve), 1e-24);
    CHECK(computer.within(curve, f.simplified, f.measured.distance));
    CHECK(!computer.within(curve, f.simplified, f.measured.distance - f.measured.tolerance * 1.0001));
  }
}

// The cache only saves runs: each search asks the same questions, so it gives the same answer.
TEST("size/budgets give the same fits in any order, a repeated one at no cost") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(120);
  const std::vector<Fit> up = ssk::pipelines::fit(computer, kShared, curve, 0.1, {4, 9, 30, 30});
  const std::vector<Fit> down = ssk::pipelines::fit(computer, kShared, curve, 0.1, {30, 9, 4});
  const std::vector<Fit> alone = ssk::pipelines::fit(computer, kShared, curve, 0.1, {9});
  CHECK(up[0].delta == down[2].delta && up[0].simplified == down[2].simplified);
  CHECK(up[1].delta == down[1].delta && up[1].delta == alone[0].delta);
  CHECK(up[2].delta == down[0].delta && up[2].simplified == down[0].simplified);
  CHECK(up[3].delta == up[2].delta && up[3].simplified == up[2].simplified);
  CHECK(up[1].runs < alone[0].runs);
}

// Not a theorem -- ours' size need not fall as delta grows -- but on a curve this smooth it does, and
// its measured error need not follow delta down at all.
TEST("size/a looser budget needs no larger delta, on a smooth curve") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(300);
  const std::vector<Fit> fits =
      ssk::pipelines::fit(computer, kShared, curve, 0.1, {2, 4, 8, 16, 32, 64, 128});
  for (std::size_t i = 1; i < fits.size(); ++i) {
    CHECK(fits[i].delta <= fits[i - 1].delta);
  }
}

TEST("size/a line fits in 2 at the tolerance, and the input's own size at no error") {
  const ssk::frechet::Frechet<2> computer;
  Curve<2> line;
  for (int i = 0; i < 50; ++i) {
    line.push_back({0.5 * i, 0.25 * i});
  }
  const double tol = 1e-9 * ssk::pipelines::diagonal(line);
  const Fit two = ssk::pipelines::fit(computer, kShared, line, 0.1, {2}).front();
  CHECK(two.simplified.size() == 2 && two.below == 0.0 && two.delta <= tol);
  CHECK(two.measured.distance <= two.measured.tolerance);

  const Curve<2> curve = wander(60);
  const Fit all = ssk::pipelines::fit(computer, kShared, curve, 0.1, {curve.size()}).front();
  CHECK(all.below == 0.0 && all.delta <= 1e-9 * ssk::pipelines::diagonal(curve));
  CHECK(all.simplified.size() <= curve.size());
}

TEST("size/ours at a baseline's size, against its error") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(256);
  for (const int m : {2, 3, 4, 5, 6}) {
    const Value base = baseline(computer, curve, m);
    const std::size_t k = ssk::pipelines::budget_of(base);
    const Fit f = ssk::pipelines::fit(computer, kShared, curve, 0.1, {k}).front();
    CHECK(f.simplified.size() <= k);
    // Theorem 3: at delta = d*, kappa is at most k, so ours has at most 2k - 2, and it fits in k
    // at a delta no larger than the one where 2k - 2 is all it is guaranteed.
    CHECK(size_at(kShared, curve, base.find("frechet")->find("distance")->num(), 0.1) <=
          2 * k - 2);
  }
}

TEST("size/a budget under 2, or a baseline without its size, is refused") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(30);
  bool refused = false;
  try {
    (void)ssk::pipelines::fit(computer, kShared, curve, 0.1, {1});
  } catch (const std::runtime_error&) {
    refused = true;
  }
  CHECK(refused);
  Value unsized = baseline(computer, curve, 2);
  unsized.obj()["stats"].obj().erase("output_size");
  refused = false;
  try {
    (void)ssk::pipelines::budget_of(unsized);
  } catch (const std::runtime_error&) {
    refused = true;
  }
  CHECK(refused);
}

TEST("size/the result document and the summary, in step 7's shape with delta_rule size") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(100);
  std::vector<ssk::pipelines::Fitted> fitted;
  std::vector<std::size_t> budgets;
  for (const int m : {2, 3, 4}) {
    fitted.push_back({baseline(computer, curve, m), {}});
    budgets.push_back(ssk::pipelines::budget_of(fitted.back().baseline));
  }
  const std::vector<Fit> fits = ssk::pipelines::fit(computer, kShared, curve, 0.1, budgets);
  for (std::size_t i = 0; i < fits.size(); ++i) {
    fitted[i].fit = fits[i];
  }

  const Value doc = ssk::json::parse(ssk::pipelines::fit_to_json(document(curve), "mopsi",
                                                                 "mopsi/x.json", kShared, 0.1,
                                                                 fits[1], fitted[1].baseline)
                                         .dump());
  CHECK(doc.find("algorithm")->str() == kShared);
  CHECK(doc.find("mode")->str() == "size-matched");
  CHECK(doc.find("delta_rule")->str() == "size");
  CHECK(doc.find("params")->find("delta")->num() == fits[1].delta);
  CHECK(doc.find("stats")->find("budget")->num() == double(budgets[1]));
  CHECK(doc.find("stats")->find("output_size")->num() == double(fits[1].simplified.size()));
  CHECK(doc.find("frechet")->find("distance")->num() == fits[1].measured.distance);
  CHECK(doc.find("frechet")->find("within_bound")->boolean());
  CHECK(doc.find("matched")->find("frechet")->dump() == fitted[1].baseline.find("frechet")->dump());
  CHECK(doc.find("points")->arr().size() == fits[1].simplified.size());
  CHECK(doc.find("name")->str() == "mopsi/7/1400");

  const Value summary = ssk::json::parse(
      ssk::pipelines::fit_summary_to_json(document(curve), "mopsi", "mopsi/x.json", kShared, 0.1,
                                          fitted)
          .dump());
  CHECK(summary.find("delta_rule")->str() == "size");
  CHECK(summary.find("eps")->num() == 0.1);
  CHECK(summary.find("input_size")->num() == 100.0);
  const Array& pairs = summary.find("pairs")->arr();
  CHECK(pairs.size() == 3);
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const Value& base = *pairs[i].find("baseline");
    const Value& ours = *pairs[i].find("ours");
    CHECK(pairs[i].find("algorithm")->str() == "douglas-peucker-n");
    CHECK(pairs[i].find("m")->num() == double(i + 2));
    CHECK(base.find("output_size")->num() == double(budgets[i]));
    CHECK(base.find("distance")->num() ==
          fitted[i].baseline.find("frechet")->find("distance")->num());
    CHECK(ours.find("delta")->num() == fits[i].delta);
    CHECK(ours.find("output_size")->num() <= double(budgets[i]));
    CHECK(ours.find("distance")->num() == fits[i].measured.distance);
    CHECK(ours.find("within_bound")->boolean());
  }
}

}  // namespace
