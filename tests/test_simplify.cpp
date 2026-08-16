// The three algorithms share one contract: endpoints kept, output a subsequence of the
// input, output in order, and a tighter parameter never keeps more points.
#include "simplify/dots.hpp"
#include "simplify/douglas_peucker.hpp"
#include "simplify/squish.hpp"
#include "test_main.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ssk::simplify::Context;
using ssk::simplify::Curve;
using ssk::simplify::Params;

// A zigzag on a rising line: the spikes are what a simplifier has to decide about.
Curve<2> zigzag(std::size_t n, double amp) {
  Curve<2> c;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({static_cast<double>(i), (i % 2 == 0) ? 0.0 : amp});
  }
  return c;
}

// Coordinates in sixteenths, so every value is an exact double and the reference indices
// below are reproducible anywhere.
Curve<2> wander(std::size_t n) {
  Curve<2> c;
  std::uint32_t s = 12345;
  double x = 0.0;
  double y = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    c.push_back({x, y});
    s = s * 1664525u + 1013904223u;
    x += static_cast<double>((s >> 16) % 32u) / 16.0;
    s = s * 1664525u + 1013904223u;
    y += static_cast<double>((s >> 16) % 32u) / 16.0 - 1.0;
  }
  return c;
}

Context uniform_time(std::size_t n) {
  Context in;
  for (std::size_t i = 0; i < n; ++i) {
    in.t.push_back(static_cast<double>(i));
  }
  return in;
}

void check_contract(const std::vector<std::size_t>& out, std::size_t n,
                    const std::string& who) {
  ::ssk::test::check(out.size() >= 2, who + ": must keep at least both endpoints");
  ::ssk::test::check(out.front() == 0, who + ": must keep the first point");
  ::ssk::test::check(out.back() == n - 1, who + ": must keep the last point");
  for (std::size_t i = 1; i < out.size(); ++i) {
    ::ssk::test::check(out[i - 1] < out[i], who + ": indices must be strictly increasing");
  }
  ::ssk::test::check(out.back() < n, who + ": index out of range");
}

// --- Douglas-Peucker ---------------------------------------------------------------

TEST("simplify/dp returns exactly the requested vertex count") {
  const Curve<2> c = wander(200);
  for (const double count : {2.0, 3.0, 17.0, 100.0, 199.0}) {
    const ssk::simplify::DouglasPeucker<2> dp(Params{{"count", count}});
    const auto out = dp.indices(c);
    check_contract(out, c.size(), "dp");
    CHECK_MSG(out.size() == static_cast<std::size_t>(count),
              "DPn delivers the budget exactly, not at most");
  }
}

TEST("simplify/dp passes a curve through when the budget cannot bind") {
  const Curve<2> c = wander(20);
  for (const double count : {20.0, 40.0, 1.0, 0.0}) {
    const ssk::simplify::DouglasPeucker<2> dp(Params{{"count", count}});
    const auto out = dp.indices(c);
    check_contract(out, c.size(), "dp");
    CHECK_MSG(out.size() == c.size(),
              "upstream copies the input through rather than failing");
  }
}

TEST("simplify/dp simplifications are nested as the budget grows") {
  const Curve<2> c = wander(120);
  std::vector<std::size_t> previous;
  for (const double count : {2.0, 5.0, 13.0, 40.0, 90.0}) {
    const ssk::simplify::DouglasPeucker<2> dp(Params{{"count", count}});
    const auto out = dp.indices(c);
    check_contract(out, c.size(), "dp");
    for (const std::size_t kept : previous) {
      bool still_there = false;
      for (const std::size_t i : out) {
        still_there = still_there || i == kept;
      }
      CHECK_MSG(still_there, "a larger budget must keep every point a smaller one kept");
    }
    previous = out;
  }
}

TEST("simplify/dp works in 3D") {
  Curve<3> c;
  for (std::size_t i = 0; i < 40; ++i) {
    const auto d = static_cast<double>(i);
    c.push_back({d, 0.0, 0.0});
  }
  c[20] = {20.0, 9.0, 9.0};
  const ssk::simplify::DouglasPeucker<3> dp(Params{{"count", 5.0}});
  const auto out = dp.indices(c);
  check_contract(out, c.size(), "dp3");
  CHECK(out.size() == 5);
  bool kept_spike = false;
  for (const std::size_t i : out) {
    kept_spike = kept_spike || i == 20;
  }
  CHECK_MSG(kept_spike, "the point furthest from the chord must survive");
}

TEST("simplify/dp simplify() returns the points the indices name") {
  const Curve<2> c = zigzag(21, 5.0);
  const ssk::simplify::DouglasPeucker<2> dp(Params{{"count", 8.0}});
  const auto idx = dp.indices(c);
  const auto pts = dp.simplify(c);
  CHECK(pts.size() == idx.size());
  for (std::size_t i = 0; i < idx.size(); ++i) {
    CHECK_CLOSE(pts[i][0], c[idx[i]][0], 1e-12);
    CHECK_CLOSE(pts[i][1], c[idx[i]][1], 1e-12);
  }
}

// --- SQUISH ------------------------------------------------------------------------

TEST("simplify/squish honours its point budget") {
  const Curve<2> c = zigzag(60, 1.0);
  const Context in = uniform_time(c.size());
  const ssk::simplify::Squish sq(Params{{"buffer_size", 10}});
  const auto out = sq.indices(c, in);
  check_contract(out, c.size(), "squish");
  CHECK_MSG(out.size() <= 10, "output must stay within the buffer size");
}

TEST("simplify/squish is monotone in the buffer size") {
  const Curve<2> c = zigzag(80, 1.0);
  const Context in = uniform_time(c.size());
  std::size_t previous = 0;
  for (const double budget : {6.0, 12.0, 25.0, 50.0}) {
    const ssk::simplify::Squish sq(Params{{"buffer_size", budget}});
    const auto out = sq.indices(c, in);
    check_contract(out, c.size(), "squish");
    CHECK_MSG(out.size() >= previous, "a bigger budget must not keep fewer points");
    previous = out.size();
  }
}

TEST("simplify/squish rejects a buffer size of 4 or less") {
  bool threw = false;
  try {
    const ssk::simplify::Squish sq(Params{{"buffer_size", 4}});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK_MSG(threw, "upstream requires buffer_size > 4");
}

// --- DOTS --------------------------------------------------------------------------

TEST("simplify/dots keeps a straight line's endpoints only") {
  Curve<2> line;
  for (std::size_t i = 0; i < 50; ++i) {
    line.push_back({static_cast<double>(i), 0.0});
  }
  const ssk::simplify::Dots dots(Params{{"lssd_threshold", 1.0}});
  const auto out = dots.indices(line, uniform_time(line.size()));
  check_contract(out, line.size(), "dots");
  CHECK(out.size() == 2);
}

TEST("simplify/dots is monotone in the LSSD threshold") {
  const Curve<2> c = zigzag(60, 1.0);
  const Context in = uniform_time(c.size());
  std::size_t previous = c.size() + 1;
  for (const double th : {0.01, 0.5, 5.0, 100.0}) {
    const ssk::simplify::Dots dots(Params{{"lssd_threshold", th}});
    const auto out = dots.indices(c, in);
    check_contract(out, c.size(), "dots");
    CHECK_MSG(out.size() <= previous, "a looser threshold must not keep more points");
    previous = out.size();
  }
}

TEST("simplify/dots normalises so epoch-ms timestamps are accepted") {
  const Curve<2> c = zigzag(30, 1.0);
  Context in;
  for (std::size_t i = 0; i < c.size(); ++i) {
    in.t.push_back(1216461503656.0 + 100.0 * static_cast<double>(i));
  }
  const ssk::simplify::Dots dots(Params{{"lssd_threshold", 1.0}});
  const auto out = dots.indices(c, in);
  check_contract(out, c.size(), "dots");
}

TEST("simplify/dots raw feed still rejects an un-normalised timestamp") {
  ssk::simplify::Dots dots(Params{{"lssd_threshold", 1.0}});
  bool threw = false;
  try {
    dots.feed(0.0, 0.0, 1216461503656.0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK_MSG(threw, "upstream's numerical guard must survive the port");
}

// --- the shared contract -----------------------------------------------------------

TEST("simplify/all reject a missing clock when they need one") {
  const Curve<2> c = zigzag(20, 1.0);
  for (int which = 0; which < 2; ++which) {
    bool threw = false;
    try {
      const auto out = (which == 0)
                           ? ssk::simplify::Dots(Params{{"lssd_threshold", 1.0}}).indices(c)
                           : ssk::simplify::Squish(Params{{"buffer_size", 8}}).indices(c);
      CHECK(out.empty());
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK_MSG(threw, "an SED-based algorithm must refuse to run without timestamps");
  }
}

TEST("simplify/all pass through a curve too short to simplify") {
  Curve<2> two{{0.0, 0.0}, {1.0, 1.0}};
  const Context in = uniform_time(2);
  CHECK(ssk::simplify::DouglasPeucker<2>(Params{{"count", 2.0}}).indices(two).size() == 2);
  CHECK(ssk::simplify::Dots(Params{{"lssd_threshold", 1.0}}).indices(two, in).size() == 2);
  CHECK(ssk::simplify::Squish(Params{{"buffer_size", 8}}).indices(two, in).size() == 2);
}

TEST("simplify/params reports a missing required parameter") {
  bool threw = false;
  try {
    const ssk::simplify::DouglasPeucker<2> dp{Params{}};
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK_MSG(threw, "a required parameter must be reported, not defaulted");
  const Params some{{"a", 2.0}};
  CHECK_CLOSE(some.get("b", 7.0), 7.0, 1e-12);
  CHECK_CLOSE(some.get("a"), 2.0, 1e-12);
}

// An algorithm that emits vertices which are nowhere in the input -- the case the interface
// has to allow, and the reason `simplify` returns points rather than positions.
class Midpoints : public ssk::simplify::Simplifier<2> {
 public:
  [[nodiscard]] std::string name() const override { return "midpoints"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& c, const Context&) const override {
    Curve<2> out{c.front()};
    for (std::size_t i = 1; i < c.size(); ++i) {
      out.push_back({(c[i - 1][0] + c[i][0]) / 2.0, (c[i - 1][1] + c[i][1]) / 2.0});
    }
    out.push_back(c.back());
    return out;
  }
};

TEST("simplify/a non-subset algorithm needs only simplify()") {
  const Curve<2> c = zigzag(9, 2.0);
  const Midpoints m;
  const auto out = m.simplify(c);
  CHECK(out.size() == c.size() + 1);
  CHECK_CLOSE(out[1][0], 0.5, 1e-12);
  CHECK_MSG(out[1][1] != c[0][1] && out[1][1] != c[1][1],
            "an invented vertex must not have to match an input point");
  const ssk::simplify::Simplifier<2>& base = m;
  CHECK(base.simplify(c).size() == out.size());
}

TEST("simplify/subset algorithms agree between simplify and indices") {
  const Curve<2> c = zigzag(41, 1.0);
  const Context in = uniform_time(c.size());
  const ssk::simplify::DouglasPeucker<2> dp(Params{{"count", 9.0}});
  const ssk::simplify::Squish sq(Params{{"buffer_size", 12}});
  const ssk::simplify::Dots dt(Params{{"lssd_threshold", 0.5}});

  for (const ssk::simplify::SubsetSimplifier<2>* algo :
       {static_cast<const ssk::simplify::SubsetSimplifier<2>*>(&dp),
        static_cast<const ssk::simplify::SubsetSimplifier<2>*>(&sq),
        static_cast<const ssk::simplify::SubsetSimplifier<2>*>(&dt)}) {
    const auto idx = algo->indices(c, in);
    const auto pts = algo->simplify(c, in);
    ::ssk::test::check(pts.size() == idx.size(), algo->name() + ": size mismatch");
    for (std::size_t i = 0; i < idx.size(); ++i) {
      ::ssk::test::check_close(pts[i][0], c[idx[i]][0], 1e-12, algo->name() + ": x");
      ::ssk::test::check_close(pts[i][1], c[idx[i]][1], 1e-12, algo->name() + ": y");
    }
  }
}

// Pinned against psimpl and the dots repo, run on this curve through the Qt shim described
// in docs/simplification.md. A port that drifts from upstream fails here.
TEST("simplify/all reproduce the upstream reference output") {
  const Curve<2> c = wander(40);
  const Context in = uniform_time(c.size());

  using Idx = std::vector<std::size_t>;
  CHECK(ssk::simplify::DouglasPeucker<2>(Params{{"count", 5.0}}).indices(c) ==
        Idx({0, 13, 17, 21, 39}));
  CHECK(ssk::simplify::DouglasPeucker<2>(Params{{"count", 12.0}}).indices(c) ==
        Idx({0, 1, 6, 8, 11, 13, 17, 21, 27, 32, 34, 39}));
  CHECK(ssk::simplify::Squish(Params{{"buffer_size", 5}}).indices(c, in) ==
        Idx({0, 9, 14, 27, 39}));
  CHECK(ssk::simplify::Squish(Params{{"buffer_size", 9}}).indices(c, in) ==
        Idx({0, 6, 10, 13, 18, 21, 27, 32, 39}));
  CHECK(ssk::simplify::Dots(Params{{"lssd_threshold", 1.0}}).indices(c, in) ==
        Idx({0, 1, 6, 8, 10, 12, 14, 17, 21, 23, 27, 29, 32, 34, 39}));
  CHECK(ssk::simplify::Dots(Params{{"lssd_threshold", 10.0}}).indices(c, in) ==
        Idx({0, 10, 14, 28, 39}));
}

TEST("simplify/curve_of bridges a trajectory document") {
  const auto doc = ssk::io::read_trajectory(ssk::json::parse(
      R"({"dim": 2, "t": [0, 1, 2], "points": [[0, 0], [1, 9], [2, 0]]})"));
  const auto c = ssk::simplify::curve_of<2>(doc);
  const auto in = ssk::simplify::context_of(doc);
  CHECK(c.size() == 3);
  CHECK_CLOSE(c[1][1], 9.0, 1e-12);
  CHECK(in.t.size() == 3);

  bool threw = false;
  try {
    ssk::simplify::curve_of<3>(doc);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK_MSG(threw, "a 2D document must not load as a 3D curve");
}

}  // namespace
