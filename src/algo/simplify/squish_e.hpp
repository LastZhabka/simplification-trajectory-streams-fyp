// SQUISH-E(lambda, mu), ported from online/SQUISH-E/SQUISH_E.cpp in
// github.com/uestc-db/traj-compression.
//
// (Muckell et al., GeoInformatica 2014.) A priority queue of points keyed by an upper bound
// on the SED their removal would introduce. `cmp_ratio` is lambda: the queue's capacity
// grows by one point per `cmp_ratio` inputs, and the cheapest point is dropped whenever it
// fills. `sed_error` is mu: once the input ends, points keep being dropped while the
// cheapest priority is at most `sed_error`.
//
// The functions keep upstream's names and bodies. Upstream's globals are locals here, its
// `Point` is `GPSPoint`, and its lat is our y and its lon our x.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

class SquishE : public SubsetSimplifier<2> {
 public:
  explicit SquishE(const Params& params);

  [[nodiscard]] std::string name() const override { return "squish-e"; }
  [[nodiscard]] bool needs_time() const override { return true; }

 protected:
  [[nodiscard]] std::vector<std::size_t> run_indices(const Curve<2>& curve,
                                                     const Context& in) const override;

 private:
  int cmp_ratio_;
  double sed_error_;
};

}  // namespace ssk::simplify
