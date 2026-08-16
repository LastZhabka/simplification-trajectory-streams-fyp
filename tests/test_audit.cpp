// Checking a result against its input. The subtlety is duplicates: a trace that repeats a
// fix, or a stationary vehicle, makes coordinate-only matching report failures that are not
// there, which is why the check uses position and timestamp together.
#include "audit.hpp"
#include "test_main.hpp"

namespace {

ssk::io::Document doc(const char* json) {
  return ssk::io::read_trajectory(ssk::json::parse(json));
}

TEST("audit/a subsequence is recognised") {
  const auto input = doc(R"({"dim": 2, "t": [0,1,2,3,4],
                             "points": [[0,0],[1,1],[2,2],[3,3],[4,4]]})");
  CHECK(ssk::pipelines::is_subsequence(
      input, doc(R"({"dim": 2, "t": [0,4], "points": [[0,0],[4,4]]})")));
  CHECK(ssk::pipelines::is_subsequence(
      input, doc(R"({"dim": 2, "t": [0,2,4], "points": [[0,0],[2,2],[4,4]]})")));
  CHECK(ssk::pipelines::is_subsequence(input, input));
}

TEST("audit/points in the wrong order are not a subsequence") {
  const auto input = doc(R"({"dim": 2, "t": [0,1,2,3,4],
                             "points": [[0,0],[1,1],[2,2],[3,3],[4,4]]})");
  CHECK_MSG(!ssk::pipelines::is_subsequence(
                input, doc(R"({"dim": 2, "t": [0,3,2,4], "points": [[0,0],[3,3],[2,2],[4,4]]})")),
            "this is exactly what DOTS' decode produces on some inputs");
}

TEST("audit/a point that is not in the input is not a subsequence") {
  const auto input = doc(R"({"dim": 2, "t": [0,1,2], "points": [[0,0],[1,1],[2,2]]})");
  CHECK(!ssk::pipelines::is_subsequence(
      input, doc(R"({"dim": 2, "t": [0,1], "points": [[0,0],[9,9]]})")));
}

TEST("audit/a repeated fix does not cause a false failure") {
  // The same coordinate twice, a second apart. Matching on position alone would take the
  // first occurrence for the second output point and then fail; with the clock it matches.
  const auto input = doc(R"({"dim": 2, "t": [0,1,2,3],
                             "points": [[5,5],[7,7],[7,7],[9,9]]})");
  CHECK_MSG(ssk::pipelines::is_subsequence(
                input, doc(R"({"dim": 2, "t": [0,2,3], "points": [[5,5],[7,7],[9,9]]})")),
            "the second occurrence of a repeated point must still match");
}

TEST("audit/a document with no clock is matched on position alone") {
  const auto input = doc(R"({"dim": 2, "points": [[0,0],[1,1],[2,2]]})");
  CHECK(ssk::pipelines::is_subsequence(input, doc(R"({"dim": 2, "points": [[0,0],[2,2]]})")));
  CHECK(!ssk::pipelines::is_subsequence(input, doc(R"({"dim": 2, "points": [[2,2],[0,0]]})")));
}

}  // namespace
