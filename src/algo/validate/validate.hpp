// Checks of the invariants ConvexShadow's guarantee rests on, made while it runs.
//
// A failed check is counted, not thrown. A corpus run is hundreds of core-hours over a million
// curves, and the failures worth finding are a handful of rounding cases among billions of
// checks: stopping at the first one would hide how many there are and where they cluster. The
// counts travel with the result instead, so a run reports "N checks, K failed" per invariant.
//
// Every test here is written independently of the code it checks, with its own predicates,
// and accepts a relative slack of kSlack so that rounding alone does not fail it.
#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ssk::validate {

struct Vec {
  double x = 0.0;
  double y = 0.0;

  bool operator==(const Vec&) const = default;
};

using Polygon = std::vector<Vec>;

inline constexpr double kSlack = 1e-9;

enum class Check : std::size_t {
  // conv(G_v), once per run, with v at the origin
  GridHull,          // the polygon is convex, counter-clockwise and the hull of the grid points
  GridCoversBall,    // it contains the ball of radius (1 + eps/2) * delta: Lemma 2's witness fits
  GridInsideDisk,    // it lies in the disk of radius (1 + eps) * delta: the error bound
  // S[p], once per candidate per input vertex
  ShadowConvex,      // the new shadow is convex and counter-clockwise
  ShadowInsideGrid,  // it lies in conv(G_v) of the vertex just read
  ShadowBehind,      // each of its vertices is seen from p through the previous shadow
  // an emitted segment, once per output segment
  SegmentStabs,      // it meets the disk of every input vertex it replaces, in order
};

inline constexpr std::size_t kChecks = 7;

struct Count {
  std::size_t checked = 0;
  std::size_t failed = 0;
};

using Stats = std::array<Count, kChecks>;

[[nodiscard]] const char* name(Check check);
[[nodiscard]] std::size_t checked(const Stats& stats);
[[nodiscard]] std::size_t failed(const Stats& stats);
void add(Stats& into, const Stats& from);

class Validator {
 public:
  void grid(const Polygon& points, const Polygon& hull, double rho, double r);
  void shadow(const Polygon& next, const Polygon& grid, Vec p, const Polygon& previous);
  void segment(Vec a, Vec b, std::span<const Vec> stabbed, double r);

  [[nodiscard]] const Stats& stats() const { return stats_; }

 private:
  void record(Check check, bool ok);

  Stats stats_{};
};

}  // namespace ssk::validate
