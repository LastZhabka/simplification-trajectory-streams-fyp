#include "simplify/squish_e.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <stdexcept>

namespace ssk::simplify {
namespace {

struct GPSPoint {
  int index;
  double lat;
  double lon;
  double time;
};

struct GPSPointWithSED {
  double priority;
  double pi;
  GPSPoint point;
};

using Queue = std::vector<GPSPointWithSED>;

GPSPointWithSED& at(Queue& Q, int index) { return Q[static_cast<std::size_t>(index)]; }

int size(const Queue& Q) { return static_cast<int>(Q.size()); }

double cacl_SED(GPSPoint s, GPSPoint m, GPSPoint e) {
  double numerator = m.time - s.time;
  double denominator = e.time - s.time;
  double time_ratio = denominator == 0 ? 1 : numerator / denominator;
  double lat = s.lat + (e.lat - s.lat) * time_ratio;
  double lon = s.lon + (e.lon - s.lon) * time_ratio;
  double lat_diff = lat - m.lat;
  double lon_diff = lon - m.lon;
  return std::sqrt(lat_diff * lat_diff + lon_diff * lon_diff);
}

void adjust_priority(Queue& Q, int pre_index, int Q_index, int succ_index) {
  if (Q_index == size(Q) - 1 || Q_index == 0) {
    return;
  }
  double p = at(Q, Q_index).pi +
             cacl_SED(at(Q, pre_index).point, at(Q, Q_index).point, at(Q, succ_index).point);
  at(Q, Q_index).priority = p;
}

void reduce(Queue& Q, int min_index, double min_p) {
  at(Q, min_index - 1).pi = std::max(min_p, at(Q, min_index - 1).pi);
  at(Q, min_index + 1).pi = std::max(min_p, at(Q, min_index + 1).pi);
  adjust_priority(Q, min_index - 2, min_index - 1, min_index + 1);
  adjust_priority(Q, min_index - 1, min_index + 1, min_index + 2);
  Q.erase(Q.begin() + min_index);
}

int find_mini_priority(Queue& Q) {
  auto to_remove = Q.end();
  int min_index = 1, k = 1;
  for (auto it = Q.begin() + 1; it != Q.end() - 1; ++it) {
    if (to_remove == Q.end() || it->priority < to_remove->priority) {
      to_remove = it;
      min_index = k;
    }
    k = k + 1;
  }
  return min_index;
}

}  // namespace

SquishE::SquishE(const Params& params)
    : cmp_ratio_(static_cast<int>(params.get("cmp_ratio"))),
      sed_error_(params.get("sed_error")) {
  if (cmp_ratio_ < 1) {
    throw std::runtime_error("squish-e: cmp_ratio must be at least 1");
  }
}

std::vector<std::size_t> SquishE::run_indices(const Curve<2>& curve, const Context& in) const {
  int capacity = 4;
  Queue Q;
  int i = 0, min_index;
  double min_p;
  while (i < static_cast<int>(curve.size())) {
    const auto ui = static_cast<std::size_t>(i);
    if (i / cmp_ratio_ >= capacity) {
      capacity = capacity + 1;
    }
    Q.push_back({DBL_MAX, 0, {i, curve[ui][1], curve[ui][0], in.t[ui]}});
    if (i > 0) {
      adjust_priority(Q, size(Q) - 3, size(Q) - 2, size(Q) - 1);
    }
    if (size(Q) == capacity) {
      min_index = find_mini_priority(Q);
      min_p = at(Q, min_index).priority;
      reduce(Q, min_index, min_p);
    }
    i = i + 1;
  }
  min_index = find_mini_priority(Q);
  min_p = at(Q, min_index).priority;
  while (min_p <= sed_error_) {
    reduce(Q, min_index, min_p);
    min_index = find_mini_priority(Q);
    min_p = at(Q, min_index).priority;
  }

  std::vector<std::size_t> out;
  out.reserve(Q.size());
  for (const GPSPointWithSED& kept : Q) {
    out.push_back(static_cast<std::size_t>(kept.point.index));
  }
  return out;
}

}  // namespace ssk::simplify
