// Checking a result document against the trajectory it claims to simplify.
//
// The sweep refuses to write a run whose indices are not increasing, so a corpus produced by
// the current driver should be clean. This exists to check one produced by an earlier driver,
// or by a hand-run of an algorithm, or after any change to a port -- verification that reads
// what is actually on disk rather than trusting what wrote it.
#pragma once

#include "io/trajectory.hpp"

namespace ssk::pipelines {

// Whether `result` is a subsequence of `input`, matching on position *and* timestamp.
//
// Greedy left-to-right matching is complete for subsequence testing, so a false answer means
// the document genuinely is not a subsequence. Matching on coordinates alone would report
// false failures wherever a trace repeats a fix -- GPS traces do, and stationary vehicles do.
[[nodiscard]] bool is_subsequence(const io::Document& input, const io::Document& result);

}  // namespace ssk::pipelines
