#include "audit.hpp"

namespace ssk::pipelines {

bool is_subsequence(const io::Document& input, const io::Document& result) {
  const bool timed = !input.t.empty() && !result.t.empty();
  std::size_t at = 0;
  for (std::size_t i = 0; i < result.points.size(); ++i) {
    while (at < input.points.size()) {
      bool same = input.points[at].size() == result.points[i].size();
      for (std::size_t k = 0; same && k < result.points[i].size(); ++k) {
        same = input.points[at][k] == result.points[i][k];
      }
      if (same && timed) {
        same = input.t[at] == result.t[i];
      }
      if (same) {
        break;
      }
      ++at;
    }
    if (at == input.points.size()) {
      return false;
    }
    ++at;
  }
  return true;
}

}  // namespace ssk::pipelines
