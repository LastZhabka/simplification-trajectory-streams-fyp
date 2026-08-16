// DOTS, ported from DotsSimplifier in github.com/caoweiquan322/dots.
//
// "Directed acyclic graph based Online Trajectory Simplification" (Cao & Li, JSS 2017). It
// keeps a frontier of candidate predecessors, joins each new point to one whose LSSD (the
// integral of squared SED along the shortcut) is under the threshold, and periodically
// minimises the accumulated error over the frontier and emits whatever prefix every
// surviving path agrees on.
//
// The streaming interface is kept, because being online is the whole point of the algorithm:
// feed(), read_index() and finish() mirror upstream's feedData/readOutputIndex/finish. The
// Simplifier::run override drives them exactly as upstream's batchDotsByIndex does.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "simplify/simplifier.hpp"

namespace ssk::simplify {

class Dots : public SubsetSimplifier<2> {
 public:
  explicit Dots(const Params& params);

  [[nodiscard]] std::string name() const override { return "dots"; }
  [[nodiscard]] bool needs_time() const override { return true; }

  // --- the online interface, as upstream ---
  void reset();
  void feed(double x, double y, double t);
  bool read_index(int& index);
  void finish();

 protected:
  [[nodiscard]] std::vector<std::size_t> run_indices(const Curve<2>& curve,
                                                     const Context& in) const override;

 private:
  void dag_search();
  [[nodiscard]] double lssd(int fst, int lst) const;
  [[nodiscard]] bool need_update_vk() const;
  void update_vk();
  void minimize_issed();
  void viterbi_decode();

  double lssd_th_;
  double lssd_upper_bound_;
  int max_vk_size_;

  std::vector<double> ptx_, pty_, ptt_;
  std::vector<double> x_sum_, y_sum_, t_sum_, x2_sum_, y2_sum_, t2_sum_, xt_sum_, yt_sum_;
  std::vector<int> v_k_, v_l_;
  std::vector<bool> terminated_;
  int num_terminated_ = 0;
  std::vector<std::vector<int>> path_k_;
  std::vector<double> issed_;
  std::vector<int> parents_;

  std::vector<int> simplified_index_;
  int input_count_ = 0;
  int output_count_ = 0;
  bool finished_ = false;
};

}  // namespace ssk::simplify
