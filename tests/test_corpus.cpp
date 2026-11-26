// What every driver shares: which files a shard runs, and a document write that leaves a whole
// file or none. On real directories under the system's temporary one.
#include "corpus.hpp"
#include "test_main.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

// A fresh directory, removed when it goes out of scope.
struct Scratch {
  fs::path dir;
  Scratch() {
    std::random_device seed;
    dir = fs::temp_directory_path() / ("ssk-test-" + std::to_string(seed()));
    fs::create_directories(dir);
  }
  ~Scratch() { fs::remove_all(dir); }
  Scratch(const Scratch&) = delete;
  Scratch& operator=(const Scratch&) = delete;
};

void touch(const fs::path& path) { std::ofstream(path) << "{}\n"; }

bool throws(int shard, int shards, const std::string& dir) {
  try {
    (void)ssk::pipelines::trajectory_files(dir, 0, shard, shards);
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

TEST("corpus/the shards hold every trajectory once, in name order, and nothing else") {
  const Scratch scratch;
  std::vector<std::string> names;
  for (int i = 0; i < 23; ++i) {
    names.push_back("mopsi-" + std::to_string(1000 + i) + ".json");
    touch(scratch.dir / names.back());
  }
  touch(scratch.dir / "index.json");
  touch(scratch.dir / "notes.txt");
  fs::create_directories(scratch.dir / "sub.json.d");
  const std::string dir = scratch.dir.string();

  for (const int shards : {1, 2, 3, 8, 23, 40}) {
    std::multiset<std::string> seen;
    for (int shard = 0; shard < shards; ++shard) {
      const std::vector<fs::path> mine = ssk::pipelines::trajectory_files(dir, 0, shard, shards);
      CHECK(std::is_sorted(mine.begin(), mine.end()));
      for (const fs::path& p : mine) {
        seen.insert(p.filename().string());
      }
    }
    CHECK_MSG(seen == std::multiset<std::string>(names.begin(), names.end()),
              std::to_string(shards) + " shards");
  }
  const std::vector<fs::path> first = ssk::pipelines::trajectory_files(dir, 5, 0, 1);
  CHECK(first.size() == 5 && first.front().filename() == names.front());
  CHECK(ssk::pipelines::trajectory_files(dir, 5, 1, 2).size() == 2);
}

// --shards 0 stepped by zero for ever, and --shard 8 --shards 8 ran shard 0's files again.
TEST("corpus/a shard outside 0 to shards - 1, or no shards, is refused") {
  const Scratch scratch;
  touch(scratch.dir / "a.json");
  const std::string dir = scratch.dir.string();
  CHECK(throws(0, 0, dir));
  CHECK(throws(0, -2, dir));
  CHECK(throws(-1, 4, dir));
  CHECK(throws(4, 4, dir));
  CHECK(!throws(3, 4, dir));
}

TEST("corpus/a document is written whole under its name, with no part file left") {
  const Scratch scratch;
  const fs::path path = scratch.dir / "deep" / "er" / "doc.json";
  const ssk::json::Value doc(ssk::json::Object{{"dim", ssk::json::Value(2)}});
  ssk::pipelines::write_document(path, doc);
  CHECK(ssk::json::parse_file(path.string()).find("dim")->num() == 2.0);
  CHECK(!fs::exists(path.string() + ".part"));
  ssk::pipelines::write_document(path, ssk::json::Value(ssk::json::Object{}));
  CHECK(ssk::json::parse_file(path.string()).obj().empty());
}

// A write that fails, here because the part file cannot be opened, used to go unnoticed and the
// rename put whatever was there in place: a disk filling up would leave a truncated document
// that --skip-existing then trusts.
TEST("corpus/a document that cannot be written throws and is not put in place") {
  const Scratch scratch;
  const fs::path path = scratch.dir / "doc.json";
  fs::create_directories(path.string() + ".part");
  bool threw = false;
  try {
    ssk::pipelines::write_document(path, ssk::json::Value(ssk::json::Object{}));
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  CHECK(!fs::exists(path));
}

TEST("corpus/the dataset is the directory's name, with or without a trailing separator") {
  CHECK(ssk::pipelines::dataset_of("data/trajectories/mopsi") == "mopsi");
  CHECK(ssk::pipelines::dataset_of("data/trajectories/mopsi/") == "mopsi");
  CHECK(ssk::pipelines::dataset_of("data/trajectories/./mot-dancetrack//") == "mot-dancetrack");
}

}  // namespace
