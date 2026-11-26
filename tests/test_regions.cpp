// The region around a vertex and the starts of a run, as the greedy stabbers build them. Both
// guarantees rest here: the error bound on the region lying within (1 + eps) delta, the size
// bound on its holding the ball of radius rho = (1 + eps/2) delta and on the starts being
// within a snap of every point Lemma 2's witness could begin at.
#include "geometry/geometry.hpp"
#include "simplify/greedy_stabber/greedy_stabber.hpp"
#include "test_main.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace geo = ssk::geometry;
using geo::Vec;
using ssk::simplify::circle;
using ssk::simplify::grid_hull;
using ssk::simplify::grid_points;
using ssk::simplify::mgon;

constexpr double kPi = std::numbers::pi;
constexpr double kUlp = std::numeric_limits<double>::epsilon();
const std::vector<double> kCorpusEps{0.5, 0.25, 0.2, 0.1, 0.05};
// Beside them, eps where rho is not a whole number of cell diagonals, as it is at every eps of
// the corpus, so cells meet the circle of radius rho in other ways.
const std::vector<double> kGridEps{0.5, 0.25, 0.2, 0.1, 0.05, 0.37, 0.3, 0.15, 0.07};

double cell_side(double delta, double eps) { return eps * delta / (2.0 * std::sqrt(2.0)); }

std::string at(double delta, double eps) {
  return "delta " + std::to_string(delta) + ", eps " + std::to_string(eps);
}

// Every eps of the corpus, and a sweep across (0, 1).
std::vector<double> sweep(double from) {
  std::vector<double> out = kCorpusEps;
  for (double eps = from; eps < 1.0; eps += 0.0137) {
    out.push_back(eps);
  }
  return out;
}

// The grid points as integer cell coordinates; empty when one is off the lattice.
std::set<std::pair<long long, long long>> cells(const std::vector<Vec>& points, double side) {
  std::set<std::pair<long long, long long>> out;
  for (const Vec q : points) {
    const double i = std::round(q.x / side);
    const double j = std::round(q.y / side);
    if (std::abs(q.x / side - i) > 1e-9 || std::abs(q.y / side - j) > 1e-9) {
      return {};
    }
    out.insert({static_cast<long long>(i), static_cast<long long>(j)});
  }
  return out;
}

// Along one axis, in cells, the coordinate of cell [k, k + 1] nearest 0 and the one furthest.
double nearest(long long k) { return std::clamp(0.0, double(k), double(k + 1)); }
double furthest(long long k) { return std::max(std::abs(double(k)), std::abs(double(k + 1))); }

bool corners_of(const std::set<std::pair<long long, long long>>& set, long long i, long long j) {
  return set.count({i, j}) && set.count({i + 1, j}) && set.count({i, j + 1}) &&
         set.count({i + 1, j + 1});
}

// cos(pi / m) >= (1 + eps/2) / (1 + eps) is what puts the ball of radius rho inside a regular
// m-gon whose corners are at (1 + eps) delta.
TEST("regions/the m-gon has the fewest sides that fit between rho and (1 + eps) delta") {
  for (const auto& [eps, m] : {std::pair{0.5, 6u}, std::pair{0.3, 7u}, std::pair{0.25, 7u},
                               std::pair{0.2, 8u}, std::pair{0.1, 11u}, std::pair{0.05, 15u}}) {
    CHECK_MSG(mgon(1.0, eps).normals.size() == m, "eps " + std::to_string(eps));
  }
  for (const double eps : sweep(0.01)) {
    const std::size_t m = mgon(1.0, eps).normals.size();
    const double ratio = (1.0 + eps / 2.0) / (1.0 + eps);
    CHECK_MSG(m >= 3 && std::cos(kPi / static_cast<double>(m)) >= ratio,
              "eps " + std::to_string(eps) + ": " + std::to_string(m) + " sides do not fit");
    CHECK_MSG(m == 3 || std::cos(kPi / static_cast<double>(m - 1)) < ratio,
              "eps " + std::to_string(eps) + ": " + std::to_string(m - 1) + " sides would");
  }
}

// The m-gon is built 4 ulps inside (1 + eps) delta, which keeps every corner within it from eps
// 0.02 up. Below, at 28 sides and more, a corner can still round past by under an ulp, which
// run()'s margin under delta, 2 ulps of it at least, covers.
TEST("regions/the m-gon holds the ball of radius rho and has its corners within (1 + eps) delta") {
  std::mt19937 rng(97);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (const double eps : sweep(0.01)) {
    for (int i = 0; i < 6; ++i) {
      const double delta = std::pow(10.0, -6.0 + 12.0 * u(rng));
      const geo::ConvexPolygon gon = mgon(delta, eps);
      const std::vector<geo::Vec> corners = geo::base(gon, {0, 0}).corners;
      CHECK(corners.size() == gon.normals.size());
      for (std::size_t k = 0; k < gon.normals.size(); ++k) {
        CHECK_CLOSE(geo::norm(gon.normals[k]), 1.0, 1e-15);
        CHECK_MSG(gon.offsets[k] >= (1.0 + eps / 2.0) * delta,
                  at(delta, eps) + ": a side cuts rho");
        const double reach = (1.0 + eps) * delta * (eps >= 0.02 ? 1.0 : 1.0 + kUlp);
        CHECK_MSG(geo::norm(corners[k]) <= reach,
                  at(delta, eps) + ": a corner past (1 + eps) delta");
      }
    }
  }
}

// Every point of the circle of radius delta is within eps delta / 2 of a start, so a witness
// starting anywhere on it is a snap away from one.
TEST("regions/M = ceil(2 pi / eps) starts on the circle of radius delta, eps delta apart at most") {
  for (const double eps : sweep(0.01)) {
    const double delta = 0.37;
    const std::vector<Vec> starts = circle(delta, eps);
    const auto m = static_cast<std::size_t>(std::ceil(2.0 * kPi / eps));
    CHECK_MSG(starts.size() == m, "eps " + std::to_string(eps));
    CHECK((starts.front() == Vec{delta, 0.0}));
    for (std::size_t s = 0; s < starts.size(); ++s) {
      const Vec next = starts[(s + 1) % starts.size()];
      CHECK_CLOSE(geo::norm(starts[s]), delta, 1e-15 * delta);
      CHECK(geo::cross(starts[s], next) > 0.0);
      CHECK_MSG(geo::norm(next - starts[s]) <= eps * delta, "eps " + std::to_string(eps));
    }
  }
}

TEST("regions/the grid starts: on the lattice, once each, v first, nearest first") {
  for (const double eps : kGridEps) {
    for (const bool boundary : {false, true}) {
      const double delta = 2.5;
      const std::vector<Vec> points = grid_points(delta, eps, boundary);
      const auto lattice = cells(points, cell_side(delta, eps));
      CHECK_MSG(lattice.size() == points.size(), at(delta, eps) + ": off the lattice or twice");
      if (!boundary) {
        CHECK((points.front() == Vec{0.0, 0.0}));
      }
      for (std::size_t k = 1; k < points.size(); ++k) {
        CHECK(geo::dot(points[k - 1], points[k - 1]) <=
              geo::dot(points[k], points[k]) * (1 + 4 * kUlp));
      }
      for (const Vec q : points) {
        CHECK(geo::norm(q) <= (1.0 + eps) * delta * (1.0 + 4.0 * kUlp));
      }
    }
  }
}

// The snap of Lemma 2: a point anywhere in the ball of radius rho lies in a cell whose four
// corners are all starts, so some start is within a cell diagonal, eps delta / 2, of it.
TEST("regions/every point of the ball of radius rho lies in a cell of four grid starts") {
  std::mt19937 rng(101);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (const double eps : kGridEps) {
    const double delta = 1.75;
    const double side = cell_side(delta, eps);
    const double rho = (1.0 + eps / 2.0) * delta;
    const auto all = cells(grid_points(delta, eps, false), side);
    for (int i = 0; i < 4000; ++i) {
      const Vec x = geo::Vec{std::cos(2 * kPi * u(rng)), std::sin(2 * kPi * u(rng))} *
                    (rho * std::sqrt(u(rng)));
      if (geo::norm(x) > rho) {
        continue;
      }
      const auto ci = static_cast<long long>(std::floor(x.x / side));
      const auto cj = static_cast<long long>(std::floor(x.y / side));
      CHECK_MSG(corners_of(all, ci, cj), "eps " + std::to_string(eps));
    }
    const double radius = rho / side;
    const auto n = static_cast<long long>(std::ceil(radius)) + 1;
    for (long long i = -n; i < n; ++i) {
      for (long long j = -n; j < n; ++j) {
        if (std::hypot(nearest(i), nearest(j)) <= radius) {
          CHECK_MSG(corners_of(all, i, j), "eps " + std::to_string(eps) + ": a cell in the ball");
        }
      }
    }
  }
}

// The boundary starts are the corners of the cells that meet the circle of radius rho, so the
// circle, where the greedy may start, is covered the same way, with O(1 / eps) starts.
TEST("regions/the boundary starts cover the circle of radius rho, a subset growing as 1 / eps") {
  std::mt19937 rng(103);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::vector<std::pair<std::size_t, std::size_t>> counts;
  for (const double eps : kGridEps) {
    const double delta = 0.9;
    const double side = cell_side(delta, eps);
    const double rho = (1.0 + eps / 2.0) * delta;
    const auto all = cells(grid_points(delta, eps, false), side);
    const auto rim = cells(grid_points(delta, eps, true), side);
    CHECK(std::includes(all.begin(), all.end(), rim.begin(), rim.end()));
    for (int i = 0; i < 4000; ++i) {
      const Vec x = geo::Vec{std::cos(2 * kPi * u(rng)), std::sin(2 * kPi * u(rng))} * rho;
      if (std::abs(geo::norm(x) - rho) > 1e-12 * rho) {
        continue;
      }
      CHECK_MSG(corners_of(rim, static_cast<long long>(std::floor(x.x / side)),
                           static_cast<long long>(std::floor(x.y / side))),
                "eps " + std::to_string(eps));
    }
    const double radius = rho / side;
    const auto meets = [radius](long long a, long long b) {
      return std::hypot(nearest(a), nearest(b)) <= radius &&
             std::hypot(furthest(a), furthest(b)) > radius;
    };
    const auto n = static_cast<long long>(std::ceil(radius)) + 1;
    for (long long i = -n; i < n; ++i) {
      for (long long j = -n; j < n; ++j) {
        if (meets(i, j)) {
          CHECK_MSG(corners_of(rim, i, j), "eps " + std::to_string(eps) + ": a cell on the circle");
        }
      }
    }
    for (const auto& [i, j] : rim) {
      const bool of_one = meets(i - 1, j - 1) || meets(i - 1, j) || meets(i, j - 1) || meets(i, j);
      CHECK_MSG(of_one, "eps " + std::to_string(eps) + ": a start of no cell on the circle");
    }
    counts.emplace_back(rim.size(), all.size());
  }
  const double rim_growth = double(counts[4].first) / double(counts[3].first);
  const double all_growth = double(counts[4].second) / double(counts[3].second);
  CHECK_MSG(rim_growth > 1.6 && rim_growth < 2.4, "boundary starts x" + std::to_string(rim_growth));
  CHECK_MSG(all_growth > 3.2 && all_growth < 4.8, "all starts x" + std::to_string(all_growth));
}

// Cell by cell rather than by sampled points, and across eps rather than at a few: a cell that
// only just crosses the circle holds a sliver of it that samples miss, and whether such cells
// exist depends on eps. A cell that touches the circle at one point is rounding's to decide, in
// the code as here, and skipped: at every eps of the corpus rho is a whole number of cell
// diagonals, so the circle runs through the lattice point (2 / eps + 1)(1, 1).
TEST("regions/every cell meeting the ball or the circle of radius rho has its corners as starts") {
  for (const double eps : sweep(0.05)) {
    const double delta = 1.0;
    const double side = cell_side(delta, eps);
    const double radius = (1.0 + eps / 2.0) * delta / side;
    const auto all = cells(grid_points(delta, eps, false), side);
    const auto rim = cells(grid_points(delta, eps, true), side);
    const auto n = static_cast<long long>(std::ceil(radius)) + 1;
    for (long long i = -n; i < n; ++i) {
      for (long long j = -n; j < n; ++j) {
        const double near = std::hypot(nearest(i), nearest(j));
        const double far = std::hypot(furthest(i), furthest(j));
        if (std::abs(near - radius) <= 1e-9 * radius || std::abs(far - radius) <= 1e-9 * radius) {
          continue;
        }
        if (near <= radius) {
          CHECK_MSG(corners_of(all, i, j), "eps " + std::to_string(eps) + ": a cell in the ball");
        }
        if (near <= radius && far > radius) {
          CHECK_MSG(corners_of(rim, i, j), "eps " + std::to_string(eps) + ": a cell on the circle");
        }
      }
    }
  }
}

TEST("regions/conv(G_v) is the hull of the grid and holds the ball of radius rho") {
  for (const double eps : kGridEps) {
    const double delta = 3.1;
    const geo::ConvexPolygon hull = grid_hull(delta, eps);
    const std::vector<Vec> corners = geo::base(hull, {0, 0}).corners;
    const std::vector<Vec> points = grid_points(delta, eps, false);
    CHECK(geo::is_convex(corners));
    CHECK_MSG(cells(corners, cell_side(delta, eps)).size() == corners.size(),
              "eps " + std::to_string(eps) + ": a corner off the grid");
    for (const Vec c : corners) {
      CHECK(std::find_if(points.begin(), points.end(), [c](Vec q) {
              return geo::norm(q - c) <= 1e-12;
            }) != points.end());
    }
    for (const Vec q : points) {
      for (std::size_t k = 0; k < hull.normals.size(); ++k) {
        CHECK(geo::dot(hull.normals[k], q) <= hull.offsets[k] + 1e-12 * delta);
      }
    }
    for (const double b : hull.offsets) {
      CHECK_MSG(b >= (1.0 + eps / 2.0) * delta, "eps " + std::to_string(eps) + ": a side cuts rho");
    }
  }
}

// A cell meeting the ball of radius rho reaches at most a diagonal, eps delta / 2, past it, so
// conv(G_v) lies within (1 + eps) delta. Where 2 / eps is a whole number, every eps of the corpus,
// a corner is exactly there, and rounding puts it up to a few ulps out: unlike the m-gon the hull
// is not shrunk, and what keeps the output inside the bound is run()'s margin under delta.
TEST("regions/conv(G_v) reaches (1 + eps) delta, to a few ulps, at every eps of the corpus") {
  std::mt19937 rng(107);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (const double eps : kCorpusEps) {
    double furthest = 0.0;
    for (int i = 0; i < 10; ++i) {
      const double delta = std::pow(10.0, -6.0 + 12.0 * u(rng));
      for (const Vec c : geo::base(grid_hull(delta, eps), {0, 0}).corners) {
        CHECK_MSG(geo::norm(c) <= (1.0 + eps) * delta * (1.0 + 4.0 * kUlp), at(delta, eps));
        furthest = std::max(furthest, geo::norm(c) / ((1.0 + eps) * delta));
      }
    }
    CHECK_MSG(furthest >= 1.0 - 4.0 * kUlp, "eps " + std::to_string(eps) + ": short of the bound");
  }
}

TEST("regions/region and starts scale with delta") {
  for (const double eps : {0.5, 0.2, 0.05}) {
    const double k = 1024.0;
    const geo::ConvexPolygon small = mgon(0.3, eps);
    const geo::ConvexPolygon large = mgon(0.3 * k, eps);
    CHECK(small.normals == large.normals);
    for (std::size_t i = 0; i < small.offsets.size(); ++i) {
      CHECK_CLOSE(large.offsets[i], k * small.offsets[i], 1e-14 * large.offsets[i]);
    }
    for (const bool boundary : {false, true}) {
      const std::vector<Vec> a = grid_points(0.3, eps, boundary);
      const std::vector<Vec> b = grid_points(0.3 * k, eps, boundary);
      CHECK(a.size() == b.size());
      for (std::size_t i = 0; i < a.size(); ++i) {
        CHECK(b[i] == a[i] * k);
      }
    }
    const std::vector<Vec> ca = circle(0.3, eps);
    const std::vector<Vec> cb = circle(0.3 * k, eps);
    for (std::size_t i = 0; i < ca.size(); ++i) {
      CHECK(cb[i] == ca[i] * k);
    }
  }
}

}  // namespace
