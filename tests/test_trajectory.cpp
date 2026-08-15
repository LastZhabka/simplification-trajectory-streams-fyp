// read_trajectory: the exported documents carry an explicit "dim", but a hand-written one
// need not, and a document whose points disagree with either must not load silently.
#include "io/trajectory.hpp"
#include "test_main.hpp"

#include <stdexcept>
#include <string>

namespace {

ssk::io::Document parse(const std::string& text) {
  return ssk::io::read_trajectory(ssk::json::parse(text));
}

bool raises(const std::string& text) {
  try {
    parse(text);
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

TEST("io/trajectory reads dim, name and points") {
  const auto doc = parse(R"({"dim": 2, "name": "mopsi/7/1400", )"
                         R"("points": [[116.318417, 39.984702], [116.31845, 39.984683]]})");
  CHECK(doc.dim == 2);
  CHECK(doc.name == "mopsi/7/1400");
  CHECK(doc.points.size() == 2);
  CHECK_CLOSE(doc.points[0][0], 116.318417, 1e-9);
  CHECK_CLOSE(doc.points[1][1], 39.984683, 1e-9);
}

TEST("io/trajectory reads 3D points") {
  const auto doc = parse(R"({"dim": 3, "points": [[1, 2, 3], [4, 5, 6]]})");
  CHECK(doc.dim == 3);
  CHECK_CLOSE(doc.points[1][2], 6.0, 1e-12);
}

TEST("io/trajectory infers the dimension when dim is absent") {
  const auto two = parse(R"({"points": [[1, 2], [3, 4]]})");
  CHECK(two.dim == 2);
  const auto three = parse(R"({"points": [[1, 2, 3], [4, 5, 6]]})");
  CHECK(three.dim == 3);
  CHECK_CLOSE(three.points[1][2], 6.0, 1e-12);
}

TEST("io/trajectory name is optional") {
  const auto doc = parse(R"({"points": [[1, 2]]})");
  CHECK(doc.name.empty());
}

TEST("io/trajectory accepts a document with no points") {
  const auto doc = parse(R"({"dim": 2, "points": []})");
  CHECK(doc.dim == 2);
  CHECK(doc.points.empty());
}

TEST("io/trajectory rejects a point that disagrees with dim") {
  CHECK_MSG(raises(R"({"dim": 2, "points": [[1, 2], [3, 4, 5]]})"),
            "a point wider than \"dim\" must raise");
}

TEST("io/trajectory inferred dimension is still enforced on later points") {
  CHECK_MSG(raises(R"({"points": [[1, 2], [3, 4, 5]]})"),
            "a later point of a different width must raise");
}

TEST("io/trajectory rejects a missing or non-array points field") {
  CHECK_MSG(raises(R"({"dim": 2})"), "no \"points\" field must raise");
  CHECK_MSG(raises(R"({"points": "1,2"})"), "a non-array \"points\" must raise");
}

TEST("io/trajectory rejects a non-numeric coordinate") {
  CHECK_MSG(raises(R"({"points": [[1, "2"]]})"), "a string coordinate must raise");
  CHECK_MSG(raises(R"({"points": [1, 2]})"), "a point that is not an array must raise");
}

TEST("io/trajectory round-trips through to_json") {
  const auto in = parse(R"({"dim": 2, "name": "spiral", "points": [[0, 0], [1.5, -2.25]]})");
  const auto out = ssk::io::read_trajectory(ssk::json::parse(ssk::io::to_json(in).dump(2)));
  CHECK(out.dim == in.dim);
  CHECK(out.name == in.name);
  CHECK(out.points.size() == in.points.size());
  CHECK_CLOSE(out.points[1][0], 1.5, 1e-12);
  CHECK_CLOSE(out.points[1][1], -2.25, 1e-12);
}

}  // namespace
