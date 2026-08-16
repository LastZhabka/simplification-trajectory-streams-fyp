// The interface every trajectory-simplification algorithm implements.
//
// Dimension is a template parameter: an algorithm that only works in the plane derives from
// Simplifier<2>, one that works in any dimension stays a template. Hyper-parameters go
// through Params, so a driver can build any algorithm from a config file without knowing
// which one it is. Context carries the per-point inputs that are not geometry -- today that
// is the timestamp the SED-based algorithms need, and nothing else needs one.
//
// `simplify` returns points, and that is the only thing every algorithm can promise: a
// min-# or optimal-placement method may emit vertices that are nowhere in the input.
//
// The three ported here do all pick a subsequence, which is worth more than the points
// alone -- positions keep the link back to the source record and its timestamp, and make
// the correspondence an error measure needs explicit. Those derive from SubsetSimplifier
// and additionally offer `indices`. Code that wants points takes a Simplifier; code that
// wants provenance asks for a SubsetSimplifier and gets it checked at compile time.
#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "io/trajectory.hpp"

namespace ssk::simplify {

template <std::size_t D>
using Point = std::array<double, D>;

template <std::size_t D>
using Curve = std::vector<Point<D>>;

struct Context {
  std::vector<double> t;
};

class Params {
 public:
  Params() = default;
  Params(std::initializer_list<std::pair<const std::string, double>> init) : values_(init) {}

  void set(const std::string& key, double value) { values_[key] = value; }

  [[nodiscard]] double get(const std::string& key) const {
    const auto it = values_.find(key);
    if (it == values_.end()) {
      throw std::runtime_error("params: missing required parameter '" + key + "'");
    }
    return it->second;
  }

  [[nodiscard]] double get(const std::string& key, double fallback) const {
    const auto it = values_.find(key);
    return (it == values_.end()) ? fallback : it->second;
  }

 private:
  std::map<std::string, double> values_;
};

template <std::size_t D>
class Simplifier {
 public:
  virtual ~Simplifier() = default;

  [[nodiscard]] virtual std::string name() const = 0;
  [[nodiscard]] virtual bool needs_time() const { return false; }

  [[nodiscard]] Curve<D> simplify(const Curve<D>& curve, const Context& in = {}) const {
    check(curve, in);
    if (curve.size() < 3) {
      return curve;
    }
    return run(curve, in);
  }

 protected:
  [[nodiscard]] virtual Curve<D> run(const Curve<D>& curve, const Context& in) const = 0;

  void check(const Curve<D>& curve, const Context& in) const {
    if (needs_time() && in.t.size() != curve.size()) {
      throw std::runtime_error(name() + ": needs one timestamp per point, got " +
                               std::to_string(in.t.size()) + " for " +
                               std::to_string(curve.size()) + " points");
    }
  }
};

// An algorithm whose output is always a subsequence of its input.
template <std::size_t D>
class SubsetSimplifier : public Simplifier<D> {
 public:
  [[nodiscard]] std::vector<std::size_t> indices(const Curve<D>& curve,
                                                 const Context& in = {}) const {
    this->check(curve, in);
    // Nothing to drop between the two endpoints, and several algorithms index [0] and [1]
    // unguarded.
    if (curve.size() < 3) {
      std::vector<std::size_t> all(curve.size());
      for (std::size_t i = 0; i < all.size(); ++i) {
        all[i] = i;
      }
      return all;
    }
    return run_indices(curve, in);
  }

 protected:
  [[nodiscard]] Curve<D> run(const Curve<D>& curve, const Context& in) const final {
    Curve<D> out;
    for (const std::size_t i : run_indices(curve, in)) {
      out.push_back(curve[i]);
    }
    return out;
  }

  [[nodiscard]] virtual std::vector<std::size_t> run_indices(const Curve<D>& curve,
                                                             const Context& in) const = 0;
};

// --- bridges to the on-disk format -------------------------------------------------

template <std::size_t D>
Curve<D> curve_of(const io::Document& doc) {
  if (doc.dim != D) {
    throw std::runtime_error("simplify: document is " + std::to_string(doc.dim) +
                             "D but the algorithm is " + std::to_string(D) + "D");
  }
  Curve<D> out;
  out.reserve(doc.points.size());
  for (const io::Point& p : doc.points) {
    Point<D> q{};
    for (std::size_t k = 0; k < D; ++k) {
      q[k] = p[k];
    }
    out.push_back(q);
  }
  return out;
}

inline Context context_of(const io::Document& doc) { return Context{doc.t}; }

}  // namespace ssk::simplify
