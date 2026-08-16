#include "simplify/squish.hpp"

#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace ssk::simplify {
namespace {

// Upstream keeps priorities as integers scaled by this, so ties break on index.
constexpr double kSedAccuracy = 1000.0;
constexpr int kInf = 1 << 30;

struct Entry {
  double x = 0.0;
  double y = 0.0;
  double t = 0.0;
  int sed = 0;
};

using Buffer = std::set<std::pair<int, int>>;  // (priority, index), lowest first

int losing_sed(const std::map<int, Entry>& points, int i) {
  auto it = points.find(i);
  const Entry& curr = it->second;
  const Entry& prev = std::prev(it)->second;
  const Entry& next = std::next(it)->second;

  const double k = (curr.t - prev.t) / (next.t - prev.t);
  const double px = (1.0 - k) * prev.x + k * next.x;
  const double py = (1.0 - k) * prev.y + k * next.y;
  const double sed = std::sqrt((curr.x - px) * (curr.x - px) + (curr.y - py) * (curr.y - py));
  return static_cast<int>(std::lround(sed * kSedAccuracy));
}

void remove_first(Buffer& buffer, std::map<int, Entry>& points) {
  const std::pair<int, int> first = *buffer.begin();

  auto it = points.find(first.second);
  const int curr_index = it->first;
  const int curr_sed = it->second.sed;
  auto prev = std::prev(it);
  const int prev_index = prev->first;
  const int prev_sed = prev->second.sed;
  auto next = std::next(it);
  const int next_index = next->first;
  const int next_sed = next->second.sed;

  points.erase(curr_index);
  buffer.erase(buffer.begin());
  // Index 0 is pinned at infinity; re-keying it would let the start point be dropped.
  if (prev_index != 0) {
    points[prev_index].sed += curr_sed;
    buffer.erase({prev_sed, prev_index});
    buffer.insert({prev_sed + curr_sed, prev_index});
  }
  points[next_index].sed += curr_sed;
  if (buffer.erase({next_sed, next_index}) != 0) {
    buffer.insert({next_sed + curr_sed, next_index});
  }
}

}  // namespace

Squish::Squish(const Params& params)
    : buffer_size_(static_cast<int>(params.get("buffer_size"))) {
  if (buffer_size_ <= 4) {
    throw std::runtime_error("squish: buffer_size must be greater than 4");
  }
}

std::vector<std::size_t> Squish::run_indices(const Curve<2>& curve, const Context& in) const {
  std::map<int, Entry> points;
  Buffer buffer;
  const int count = static_cast<int>(curve.size());

  points[0] = Entry{curve[0][0], curve[0][1], in.t[0], kInf};
  buffer.insert({kInf, 0});
  points[1] = Entry{curve[1][0], curve[1][1], in.t[1], 0};

  for (int i = 2; i < count; ++i) {
    const int last_key = points.rbegin()->first;
    points[i] = Entry{curve[static_cast<std::size_t>(i)][0],
                      curve[static_cast<std::size_t>(i)][1],
                      in.t[static_cast<std::size_t>(i)], 0};

    const int sed = losing_sed(points, last_key);
    points[last_key].sed += sed;
    buffer.insert({points[last_key].sed, last_key});

    if (static_cast<int>(buffer.size()) >= buffer_size_) {
      remove_first(buffer, points);
    }
  }

  std::vector<std::size_t> out;
  out.reserve(points.size());
  for (const auto& [index, entry] : points) {
    out.push_back(static_cast<std::size_t>(index));
  }
  return out;
}

}  // namespace ssk::simplify
