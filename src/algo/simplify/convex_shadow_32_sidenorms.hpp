// ConvexShadow32 with each shadow kept as what bounds it instead of as a polygon. Every region
// is the same regular m-gon moved to a vertex, so a shadow is bounded by those m side normals,
// each at an offset, and by two rays from its start p; a candidate is p, the m offsets, the two
// rays and the two tangents for the next step, and the shadow's corners are rebuilt from them
// when needed. A step keeps the offsets p is beyond, at the minimum with the new m-gon's.
// docs/convex-shadow-32.md §8.
#pragma once

#include <string>

#include "simplify/simplifier.hpp"
#include "validate/validate.hpp"

namespace ssk::simplify {

class ConvexShadow32SideNorms : public Simplifier<2> {
 public:
  explicit ConvexShadow32SideNorms(const Params& params, validate::Validator* validator = nullptr);

  [[nodiscard]] std::string name() const override { return "convex-shadow-32-sidenorms"; }

 protected:
  [[nodiscard]] Curve<2> run(const Curve<2>& curve, const Context& in) const override;

 private:
  double delta_;
  double eps_;
  validate::Validator* validator_;
};

}  // namespace ssk::simplify
