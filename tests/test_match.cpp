// ssk_match, steps 7 and 9: ours at a baseline's measured error, d*. Every number of the
// comparison passes through match(), and through the two documents it writes.
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
using ssk::pipelines::Match;
using ssk::simplify::Curve;

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

// A baseline's frechet-distances document as ssk_frechet writes it: DPn at N / 2^m, measured.
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

double distance_of(const Value& doc) { return doc.find("frechet")->find("distance")->num(); }

TEST("match/ours runs at d* / (1 + eps), or at d* unscaled, as the implementation named") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(120);
  for (const int m : {2, 4}) {
    const Value base = baseline(computer, curve, m);
    const double d = distance_of(base);
    for (const std::string& ours : ssk::pipelines::kOurs) {
      for (const bool unscaled : {false, true}) {
        const double eps = 0.25;
        const Match got = ssk::pipelines::match(computer, ours, curve, base, eps, 1e-9, unscaled);
        const std::string where = ours + (unscaled ? " unscaled" : "") + ", m" + std::to_string(m);
        CHECK_MSG(got.ours == ours && got.unscaled == unscaled, where);
        CHECK_MSG(got.delta == (unscaled ? d : d / (1.0 + eps)), where + ": delta");
        CHECK_MSG(got.simplified ==
                      ssk::pipelines::make_ours(ours, {{"delta", got.delta}, {"eps", eps}})
                          ->simplify(curve),
                  where + ": not the implementation's own output");
        CHECK_MSG(got.within_bound, where + ": over (1 + eps) delta");
        CHECK(got.within_bound == computer.within(curve, got.simplified, (1.0 + eps) * got.delta));
        CHECK_CLOSE(got.measured.tolerance, 1e-9 * ssk::pipelines::diagonal(curve), 1e-24);
        CHECK(computer.within(curve, got.simplified, got.measured.distance));
        CHECK(!computer.within(curve, got.simplified,
                               got.measured.distance - got.measured.tolerance * 1.0001));
        CHECK(distance_of(got.baseline) == d);
      }
    }
  }
}

// Scaled, ours is within d* by its guarantee: never worse than the baseline it is matched with.
TEST("match/scaled, ours is within the baseline's own measured error") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(200);
  for (const int m : {1, 2, 3, 4, 5}) {
    const Value base = baseline(computer, curve, m);
    for (const double eps : {0.5, 0.1}) {
      const Match got =
          ssk::pipelines::match(computer, "greedy-stabbing-shared", curve, base, eps, 1e-9, false);
      CHECK_MSG(computer.within(curve, got.simplified, distance_of(base)),
                "m" + std::to_string(m) + ", eps " + std::to_string(eps));
    }
  }
}

// Shared is naive to the bit, so only its name tells which of the two a run used.
TEST("match/each implementation of ours is built by its own name") {
  for (const std::string& ours : ssk::pipelines::kOurs) {
    const auto built = ssk::pipelines::make_ours(ours, {{"delta", 1.0}, {"eps", 0.5}});
    CHECK_MSG(built->name() == ours, ours);
  }
  CHECK(ssk::pipelines::kOurs.front() == "greedy-stabbing-shared");
}

TEST("match/an unknown implementation, or a baseline without its measurement, is refused") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(30);
  const auto raises = [&](const std::string& ours, const Value& base) {
    try {
      (void)ssk::pipelines::match(computer, ours, curve, base, 0.5, 1e-9, false);
    } catch (const std::runtime_error&) {
      return true;
    }
    return false;
  };
  CHECK(raises("greedy-stabbing-psychic", baseline(computer, curve, 2)));
  Value unmeasured = baseline(computer, curve, 2);
  unmeasured.obj().erase("frechet");
  CHECK(raises("greedy-stabbing-shared", unmeasured));
  Value no_distance = baseline(computer, curve, 2);
  no_distance.obj()["frechet"].obj().erase("distance");
  CHECK(raises("greedy-stabbing-shared", no_distance));
}

TEST("match/the result document: ours' points and numbers, beside the baseline it was matched to") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(80);
  const Value base = baseline(computer, curve, 3);
  for (const bool unscaled : {false, true}) {
    const Match got =
        ssk::pipelines::match(computer, "greedy-stabbing-naive", curve, base, 0.2, 1e-9, unscaled);
    const Value doc = ssk::json::parse(
        ssk::pipelines::match_to_json(document(curve), "mopsi", "mopsi/x.json", 0.2, got).dump());
    CHECK(doc.find("algorithm")->str() == "greedy-stabbing-naive");
    CHECK(doc.find("mode")->str() == "error-matched");
    CHECK(doc.find("delta_rule")->str() == (unscaled ? "d*" : "d*/(1+eps)"));
    CHECK(doc.find("params")->find("delta")->num() == got.delta);
    CHECK(doc.find("params")->find("eps")->num() == 0.2);
    CHECK(doc.find("stats")->find("input_size")->num() == 80.0);
    CHECK(doc.find("stats")->find("output_size")->num() == double(got.simplified.size()));
    CHECK_CLOSE(doc.find("stats")->find("compression")->num(), 80.0 / double(got.simplified.size()),
                1e-12);
    CHECK(doc.find("frechet")->find("distance")->num() == got.measured.distance);
    CHECK(doc.find("frechet")->find("bound")->num() == 1.2 * got.delta);
    CHECK(doc.find("frechet")->find("within_bound")->boolean());
    CHECK(doc.find("dataset")->str() == "mopsi" && doc.find("source")->str() == "mopsi/x.json");
    CHECK(doc.find("name")->str() == "mopsi/7/1400");
    const Array& points = doc.find("points")->arr();
    CHECK(points.size() == got.simplified.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      CHECK(points[i].arr()[0].num() == got.simplified[i][0]);
      CHECK(points[i].arr()[1].num() == got.simplified[i][1]);
    }
    const Value* matched = doc.find("matched");
    for (const char* key : {"algorithm", "mode", "params", "stats", "frechet", "simplified"}) {
      CHECK_MSG(matched->find(key) != nullptr &&
                    matched->find(key)->dump() == base.find(key)->dump(),
                std::string("matched.") + key);
    }
    CHECK(doc.find("t") == nullptr);
  }
}

TEST("match/the summary: one pair per baseline run, each side's size and error") {
  const ssk::frechet::Frechet<2> computer;
  const Curve<2> curve = wander(150);
  std::vector<Match> matches;
  for (const int m : {2, 3, 4}) {
    matches.push_back(ssk::pipelines::match(computer, "greedy-stabbing-shared", curve,
                                            baseline(computer, curve, m), 0.1, 1e-9, false));
  }
  const Value doc = ssk::json::parse(ssk::pipelines::summary_to_json(
      document(curve), "mopsi", "mopsi/x.json", "greedy-stabbing-shared", 0.1, false, matches)
                                         .dump());
  CHECK(doc.find("ours")->str() == "greedy-stabbing-shared");
  CHECK(doc.find("eps")->num() == 0.1);
  CHECK(doc.find("delta_rule")->str() == "d*/(1+eps)");
  CHECK(doc.find("input_size")->num() == 150.0);
  CHECK(doc.find("name")->str() == "mopsi/7/1400");
  const Array& pairs = doc.find("pairs")->arr();
  CHECK(pairs.size() == 3);
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const Value& pair = pairs[i];
    const Value& base = matches[i].baseline;
    CHECK(pair.find("algorithm")->str() == "douglas-peucker-n");
    CHECK(pair.find("m")->num() == double(i + 2));
    CHECK(pair.find("baseline")->find("output_size")->num() ==
          base.find("stats")->find("output_size")->num());
    CHECK(pair.find("baseline")->find("distance")->num() == distance_of(base));
    CHECK(pair.find("ours")->find("delta")->num() == matches[i].delta);
    CHECK(pair.find("ours")->find("output_size")->num() == double(matches[i].simplified.size()));
    CHECK(pair.find("ours")->find("distance")->num() == matches[i].measured.distance);
    CHECK(pair.find("ours")->find("within_bound")->boolean());
  }
  const std::vector<Match> unscaled{ssk::pipelines::match(
      computer, "greedy-stabbing-shared", curve, baseline(computer, curve, 2), 0.1, 1e-9, true)};
  const Value doc_unscaled = ssk::json::parse(
      ssk::pipelines::summary_to_json(document(curve), "mopsi", "mopsi/x.json",
                                      "greedy-stabbing-shared", 0.1, true, unscaled)
          .dump());
  CHECK(doc_unscaled.find("delta_rule")->str() == "d*");
  CHECK(doc_unscaled.find("pairs")->arr()[0].find("ours")->find("delta")->num() ==
        distance_of(unscaled[0].baseline));
}

}  // namespace
