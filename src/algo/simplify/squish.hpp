// SQUISH, ported from SquishBatchSimplifier in github.com/caoweiquan322/dots.
//
// Plain SQUISH (Muckell et al., COM.Geo 2011): a fixed-size buffer of points, each keyed by
// an upper bound on the SED its removal would introduce. When the buffer fills, the
// cheapest point is dropped and its error is handed to both neighbours.
//
// The knob is a point budget, not an error bound -- upstream implements SQUISH, not
// SQUISH-E, so there is no mu.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

class Squish : public SubsetSimplifier<2> {
 public:
  explicit Squish(const Params& params);

  [[nodiscard]] std::string name() const override { return "squish"; }
  [[nodiscard]] bool needs_time() const override { return true; }

 protected:
  [[nodiscard]] std::vector<std::size_t> run_indices(const Curve<2>& curve,
                                                     const Context& in) const override;

 private:
  int buffer_size_;
};

}  // namespace ssk::simplify
