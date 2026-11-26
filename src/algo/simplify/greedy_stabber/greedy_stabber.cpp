#include "simplify/greedy_stabber/greedy_stabber.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "geometry/geometry.hpp"

namespace ssk::simplify {
namespace {

using geometry::Vec;

// Moved off ceil(pi / acos) either way where acos rounds.
std::size_t sides(double eps) {
  const double ratio = (1.0 + eps / 2.0) / (1.0 + eps);
  const auto reaches = [ratio](std::size_t m) {
    return std::cos(std::numbers::pi / static_cast<double>(m)) >= ratio;
  };
  auto m = std::max<std::size_t>(
      3, static_cast<std::size_t>(std::ceil(std::numbers::pi / std::acos(ratio))));
  while (!reaches(m)) {
    ++m;
  }
  while (m > 3 && reaches(m - 1)) {
    --m;
  }
  return m;
}

double cell_side(double delta, double eps) { return eps * delta / (2.0 * std::sqrt(2.0)); }

// In cells around v, nearest v first.
std::vector<Vec> grid(double eps, bool boundary) {
  const double radius = (1.0 + eps / 2.0) * 2.0 * std::sqrt(2.0) / eps;
  const auto n = static_cast<long long>(std::ceil(radius));
  const auto near = [](long long k) { return std::clamp(0LL, k, k + 1); };
  const auto far = [](long long k) { return std::max(k * k, (k + 1) * (k + 1)); };
  std::vector<std::array<long long, 3>> corners;
  for (long long i = -n; i < n; ++i) {
    for (long long j = -n; j < n; ++j) {
      const auto near2 = static_cast<double>(near(i) * near(i) + near(j) * near(j));
      const auto far2 = static_cast<double>(far(i) + far(j));
      if (near2 <= radius * radius && (!boundary || far2 > radius * radius)) {
        for (const long long a : {i, i + 1}) {
          for (const long long b : {j, j + 1}) {
            corners.push_back({a * a + b * b, a, b});
          }
        }
      }
    }
  }
  std::sort(corners.begin(), corners.end());
  corners.erase(std::unique(corners.begin(), corners.end()), corners.end());
  std::vector<Vec> out;
  for (const auto& [distance, i, j] : corners) {
    out.push_back({static_cast<double>(i), static_cast<double>(j)});
  }
  return out;
}

std::vector<Vec> scaled(std::vector<Vec> points, double s) {
  for (Vec& q : points) {
    q = q * s;
  }
  return points;
}

}  // namespace

void check_params(const std::string& name, double delta, double eps) {
  if (delta <= 0.0) {
    throw std::runtime_error(name + ": delta must be positive");
  }
  if (eps <= 0.0 || eps >= 1.0) {
    throw std::runtime_error(name + ": eps must lie in (0, 1)");
  }
}

// A few ulps inside (1 + eps) * delta, as every test near the error bound errs inward: cos and
// sin round. Below eps 0.02, at 28 sides and more, a corner can still round past by under an
// ulp; run()'s margin, 2 ulps of delta at least, covers it. conv(G_v) is not shrunk at all.
geometry::ConvexPolygon mgon(double delta, double eps) {
  return geometry::regular_polygon(
      sides(eps), (1.0 + eps) * delta * (1.0 - 4.0 * std::numeric_limits<double>::epsilon()));
}

std::vector<Vec> circle(double delta, double eps) {
  const auto count = static_cast<std::size_t>(std::ceil(2.0 * std::numbers::pi / eps));
  std::vector<Vec> out;
  for (std::size_t s = 0; s < count; ++s) {
    const double angle =
        2.0 * std::numbers::pi * static_cast<double>(s) / static_cast<double>(count);
    out.push_back({delta * std::cos(angle), delta * std::sin(angle)});
  }
  return out;
}

// The hull is taken in cells, where grid points in a line are exactly in a line, and scaled after.
geometry::ConvexPolygon grid_hull(double delta, double eps) {
  return geometry::from_corners(
      scaled(geometry::full_convex_hull(grid(eps, false), false), cell_side(delta, eps)));
}

std::vector<Vec> grid_points(double delta, double eps, bool boundary) {
  return scaled(grid(eps, boundary), cell_side(delta, eps));
}

}  // namespace ssk::simplify
