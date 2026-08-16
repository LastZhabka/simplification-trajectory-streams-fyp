#include "simplify/dots.hpp"

#include <cmath>
#include <stdexcept>

namespace ssk::simplify {

Dots::Dots(const Params& params)
    : lssd_th_(params.get("lssd_threshold")),
      lssd_upper_bound_(lssd_th_ * params.get("k", 2.0)),
      max_vk_size_(static_cast<int>(params.get("max_vk_size", 1e6))) {}

void Dots::reset() {
  ptx_.clear();
  pty_.clear();
  ptt_.clear();
  x_sum_.clear();
  y_sum_.clear();
  t_sum_.clear();
  x2_sum_.clear();
  y2_sum_.clear();
  t2_sum_.clear();
  xt_sum_.clear();
  yt_sum_.clear();
  v_k_.clear();
  v_l_.clear();
  terminated_.clear();
  num_terminated_ = 0;
  path_k_.clear();
  issed_.clear();
  parents_.clear();
  simplified_index_.clear();
  input_count_ = 0;
  output_count_ = 0;
  finished_ = false;
}

void Dots::feed(double x, double y, double t) {
  if (finished_) {
    throw std::runtime_error("dots: cannot feed after finish(); call reset() first");
  }
  ptx_.push_back(x);
  pty_.push_back(y);
  ptt_.push_back(t);

  if (ptx_.size() == 1) {
    // Upstream's guard. The LSSD closed form below accumulates t*t and x*t, so a raw epoch
    // timestamp or a projected coordinate destroys the precision it needs.
    if (std::abs(t) > 3600.0 * 24.0 * 365.0) {
      throw std::runtime_error("dots: first timestamp too large; normalise the input first");
    }
    if (std::abs(x) > 2000.0 * 1000.0 || std::abs(y) > 2000.0 * 1000.0) {
      throw std::runtime_error("dots: first position too large; normalise the input first");
    }
    x_sum_.push_back(x);
    y_sum_.push_back(y);
    t_sum_.push_back(t);
    x2_sum_.push_back(x * x);
    y2_sum_.push_back(y * y);
    t2_sum_.push_back(t * t);
    xt_sum_.push_back(x * t);
    yt_sum_.push_back(y * t);

    v_k_.push_back(0);
    terminated_.push_back(false);
    num_terminated_ = 0;
    path_k_.push_back({0});

    input_count_ = 1;
    output_count_ = 0;
    simplified_index_.push_back(0);
  } else {
    const std::size_t n = x_sum_.size() - 1;
    x_sum_.push_back(x_sum_[n] + x);
    y_sum_.push_back(y_sum_[n] + y);
    t_sum_.push_back(t_sum_[n] + t);
    x2_sum_.push_back(x2_sum_[n] + x * x);
    y2_sum_.push_back(y2_sum_[n] + y * y);
    t2_sum_.push_back(t2_sum_[n] + t * t);
    xt_sum_.push_back(xt_sum_[n] + x * t);
    yt_sum_.push_back(yt_sum_[n] + y * t);
  }
  issed_.push_back(0.0);
  parents_.push_back(-1);
}

bool Dots::read_index(int& index) {
  if (!finished_ && output_count_ >= static_cast<int>(simplified_index_.size())) {
    dag_search();
  }
  if (output_count_ < static_cast<int>(simplified_index_.size())) {
    index = simplified_index_[static_cast<std::size_t>(output_count_)];
    ++output_count_;
    return true;
  }
  return false;
}

void Dots::finish() {
  if (!finished_) {
    finished_ = true;
    dag_search();
  }
}

double Dots::lssd(int fst, int lst) const {
  if (fst + 1 >= lst) {
    return 0.0;
  }
  const auto f = static_cast<std::size_t>(fst);
  const auto l = static_cast<std::size_t>(lst);
  const std::size_t p = l - 1;

  const double c1x = ptx_[f] * ptt_[l] - ptx_[l] * ptt_[f];
  const double c2x = c1x * c1x;
  const double c3x = ptt_[l] - ptt_[f];
  const double c4x = c3x * c3x;
  const double c5x = ptx_[l] - ptx_[f];
  const double c6x = c5x * c5x;

  const double c1y = pty_[f] * ptt_[l] - pty_[l] * ptt_[f];
  const double c2y = c1y * c1y;
  const double c3y = c3x;
  const double c4y = c3y * c3y;
  const double c5y = pty_[l] - pty_[f];
  const double c6y = c5y * c5y;

  const auto span = static_cast<double>(p - f);
  return span * c2x / c4x + c6x / c4x * (t2_sum_[p] - t2_sum_[f]) +
         (x2_sum_[p] - x2_sum_[f]) + 2 * c1x * c5x / c4x * (t_sum_[p] - t_sum_[f]) -
         2 * c1x / c3x * (x_sum_[p] - x_sum_[f]) - 2 * c5x / c3x * (xt_sum_[p] - xt_sum_[f]) +
         span * c2y / c4y + c6y / c4y * (t2_sum_[p] - t2_sum_[f]) +
         (y2_sum_[p] - y2_sum_[f]) + 2 * c1y * c5y / c4y * (t_sum_[p] - t_sum_[f]) -
         2 * c1y / c3y * (y_sum_[p] - y_sum_[f]) - 2 * c5y / c3y * (yt_sum_[p] - yt_sum_[f]);
}

bool Dots::need_update_vk() const {
  return static_cast<int>(v_k_.size()) == num_terminated_ ||
         static_cast<int>(v_l_.size()) >= max_vk_size_;
}

void Dots::update_vk() {
  std::vector<std::vector<int>> new_path;
  for (const int leaf : v_l_) {
    int index_k = -1;
    const int parent_pos = parents_[static_cast<std::size_t>(leaf)];
    for (std::size_t m = 0; m < v_k_.size(); ++m) {
      if (v_k_[m] == parent_pos) {
        index_k = static_cast<int>(m);
        break;
      }
    }
    std::vector<int> p = path_k_[static_cast<std::size_t>(index_k)];
    p.push_back(leaf);
    new_path.push_back(p);
  }
  path_k_ = new_path;

  v_k_ = v_l_;
  terminated_.assign(v_k_.size(), false);
  num_terminated_ = 0;
  v_l_.clear();
}

void Dots::minimize_issed() {
  for (const int i : v_l_) {
    const auto ui = static_cast<std::size_t>(i);
    double min_distance = issed_[ui];
    int min_parent = parents_[ui];
    for (const int j : v_k_) {
      const double local = lssd(j, i);
      const double distance = issed_[static_cast<std::size_t>(j)] + local;
      if (local < lssd_th_ && distance < min_distance) {
        min_distance = distance;
        min_parent = j;
      }
    }
    issed_[ui] = min_distance;
    parents_[ui] = min_parent;
  }
}

void Dots::viterbi_decode() {
  const auto input_layer = static_cast<int>(path_k_[0].size());
  auto start_layer = static_cast<int>(simplified_index_.size());
  while (start_layer < input_layer) {
    const int reference = path_k_[0][static_cast<std::size_t>(start_layer)];
    bool equal = true;
    for (const std::vector<int>& path : path_k_) {
      if (path[static_cast<std::size_t>(start_layer)] != reference) {
        equal = false;
        break;
      }
    }
    if (!equal) {
      break;
    }
    simplified_index_.push_back(reference);
    ++start_layer;
  }
}

void Dots::dag_search() {
  const auto num_points = static_cast<int>(ptx_.size());

  if (finished_) {
    while (input_count_ < num_points) {
      for (int i = input_count_; i < num_points; ++i) {
        const auto ui = static_cast<std::size_t>(i);
        if (parents_[ui] < 0) {
          for (std::size_t j = 0; j < v_k_.size(); ++j) {
            const int j_index = v_k_[j];
            if (terminated_[j]) {
              continue;
            }
            const double distance = lssd(j_index, i);
            if (distance < lssd_th_) {
              v_l_.push_back(i);
              issed_[ui] = issed_[static_cast<std::size_t>(j_index)] + distance;
              parents_[ui] = j_index;
              break;
            }
            if (distance > lssd_upper_bound_) {
              terminated_[j] = true;
              ++num_terminated_;
              if (need_update_vk()) {
                break;
              }
            }
          }
        }
        if (parents_[ui] >= 0 && input_count_ == i) {
          ++input_count_;
        }
        if (need_update_vk()) {
          break;
        }
      }
      minimize_issed();
      update_vk();
      // Upstream omits viterbi decoding here; the whole path is decoded below instead.
    }

    std::vector<int> temp;
    int idx = num_points - 1;
    int current_index = -1;
    if (!simplified_index_.empty()) {
      current_index = simplified_index_.back();
    }
    while (idx > current_index) {
      temp.push_back(idx);
      idx = parents_[static_cast<std::size_t>(idx)];
    }
    for (auto k = static_cast<int>(temp.size()) - 1; k >= 0; --k) {
      simplified_index_.push_back(temp[static_cast<std::size_t>(k)]);
    }
    return;
  }

  while (true) {
    bool vk_updated = false;
    if (input_count_ >= num_points) {
      break;
    }
    for (int i = input_count_; i < num_points; ++i) {
      const auto ui = static_cast<std::size_t>(i);
      if (parents_[ui] < 0) {
        for (std::size_t j = 0; j < v_k_.size(); ++j) {
          const int j_index = v_k_[j];
          if (terminated_[j]) {
            continue;
          }
          const double distance = lssd(j_index, i);
          if (distance < lssd_th_) {
            v_l_.push_back(i);
            issed_[ui] = issed_[static_cast<std::size_t>(j_index)] + distance;
            parents_[ui] = j_index;
            if (need_update_vk()) {
              minimize_issed();
              update_vk();
              viterbi_decode();
              vk_updated = true;
            }
            break;
          }
          if (distance > lssd_upper_bound_) {
            terminated_[j] = true;
            ++num_terminated_;
            if (need_update_vk()) {
              minimize_issed();
              update_vk();
              viterbi_decode();
              vk_updated = true;
              break;
            }
          }
        }
      }
      if (parents_[ui] >= 0 && input_count_ == i) {
        ++input_count_;
      }
      if (vk_updated) {
        break;
      }
    }
    if (!vk_updated) {
      break;
    }
  }
}

std::vector<std::size_t> Dots::run_indices(const Curve<2>& curve, const Context& in) const {
  // The guards in feed() reject raw epoch timestamps and projected coordinates outright.
  // Translating to the first point is what upstream's message asks the caller to do, and
  // LSSD is invariant under it.
  const double x0 = curve[0][0];
  const double y0 = curve[0][1];
  const double t0 = in.t[0];

  Dots run_state(*this);
  run_state.reset();

  std::vector<std::size_t> out;
  int idx = 0;
  for (std::size_t i = 0; i < curve.size(); ++i) {
    run_state.feed(curve[i][0] - x0, curve[i][1] - y0, in.t[i] - t0);
    if (run_state.read_index(idx)) {
      out.push_back(static_cast<std::size_t>(idx));
    }
  }
  run_state.finish();
  while (run_state.read_index(idx)) {
    out.push_back(static_cast<std::size_t>(idx));
  }
  return out;
}

}  // namespace ssk::simplify
