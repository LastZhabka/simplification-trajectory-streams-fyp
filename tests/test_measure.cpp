// Measuring a simplification against its input, and the document that records it.
#include "measure.hpp"
#include "test_main.hpp"

#include <cmath>
#include <string>

namespace {

using ssk::frechet::Curve;
using ssk::frechet::Frechet;
using ssk::pipelines::measure;

Curve<2> zigzag(std::size_t n, double amp) {
  Curve<2> c;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({static_cast<double>(i), (i % 2 == 0) ? 0.0 : amp});
  }
  return c;
}

TEST("measure/the tolerance scales with the curve, not the units") {
  const Frechet<2> f;
  Curve<2> small, big;
  for (std::size_t i = 0; i < 20; ++i) {
    const auto t = static_cast<double>(i);
    small.push_back({t * 1e-5, 0.0});
    big.push_back({t * 1e5, 0.0});
  }
  const Curve<2> small_ends{small.front(), small.back()};
  const Curve<2> big_ends{big.front(), big.back()};

  const auto a = measure(f, small, small_ends, 1e-9);
  const auto b = measure(f, big, big_ends, 1e-9);
  CHECK_MSG(b.tolerance > a.tolerance * 1e9,
            "a curve ten orders of magnitude larger must get a proportionally larger tolerance");
  CHECK_CLOSE(ssk::pipelines::diagonal(big), 19e5, 1.0);
}

TEST("measure/a known distance comes back") {
  const Frechet<2> f;
  // A zigzag of amplitude 2 against its own endpoints: the far corner is 2 away.
  const Curve<2> c = zigzag(5, 2.0);
  const Curve<2> ends{c.front(), c.back()};
  const auto m = measure(f, c, ends, 1e-12);
  CHECK(m.input_points == 5);
  CHECK(m.output_points == 2);
  CHECK_CLOSE(m.distance, 2.0, 1e-6);
}

TEST("measure/a degenerate curve still terminates") {
  const Frechet<2> f;
  const Curve<2> point{{3.0, 3.0}, {3.0, 3.0}, {3.0, 3.0}};
  const auto m = measure(f, point, Curve<2>{{3.0, 3.0}, {3.0, 3.0}}, 1e-9);
  CHECK_MSG(m.tolerance > 0.0, "a curve with no extent needs a positive tolerance");
  CHECK_CLOSE(m.distance, 0.0, 1e-6);
}

TEST("measure/the document carries the measurement and what produced it") {
  const auto input = ssk::io::read_trajectory(ssk::json::parse(
      R"({"dim": 2, "name": "toy/1", "points": [[0,0],[1,2],[2,0],[3,2],[4,0]]})"));
  const auto simplified = ssk::json::parse(
      R"({"dim": 2, "algorithm": "dots", "mode": "budget",
          "params": {"lssd_threshold": 0.25},
          "stats": {"input_size": 5, "output_size": 2, "m": 2, "target": 2},
          "points": [[0,0],[4,0]]})");

  ssk::pipelines::Measurement m;
  m.distance = 2.0;
  m.tolerance = 1e-6;
  const auto doc = ssk::json::parse(
      ssk::pipelines::to_json(input, simplified, "toy", "toy/toy-1.json",
                              "dots/toy/m2/toy-1.json", m)
          .dump());

  CHECK(doc.find("dataset")->str() == "toy");
  CHECK(doc.find("source")->str() == "toy/toy-1.json");
  CHECK(doc.find("simplified")->str() == "dots/toy/m2/toy-1.json");
  CHECK(doc.find("name")->str() == "toy/1");
  CHECK_CLOSE(doc.find("frechet")->find("distance")->num(), 2.0, 1e-12);
  CHECK_CLOSE(doc.find("frechet")->find("tolerance")->num(), 1e-6, 1e-18);
  CHECK(doc.find("frechet")->find("computer")->str() == "dv-gis-cup-2017");

  CHECK_MSG(doc.find("algorithm")->str() == "dots", "the algorithm must travel with the result");
  CHECK(doc.find("mode")->str() == "budget");
  CHECK_MSG(doc.find("params")->find("lssd_threshold")->num() == 0.25,
            "the hyper-parameters must travel too, so a table needs only this document");
  CHECK(doc.find("stats")->find("m")->num() == 2.0);
  CHECK_MSG(doc.find("points") == nullptr, "the geometry stays in the simplification");
}

}  // namespace
