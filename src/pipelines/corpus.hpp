#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "io/json.hpp"

namespace ssk::pipelines {

inline const std::vector<std::string> kBaselines{"douglas-peucker-n", "squish", "squish-e",
                                                 "dots"};

// Compression rates 2^first to 2^last: at each m in between, a trajectory of N points is asked
// for ceil(N / 2^m) vertices, and its documents go under m<m>/.
struct Rates {
  int first = 2;
  int last = 6;
};

// `FIRST..LAST`, as --rates takes it; nothing unless 1 <= FIRST <= LAST.
std::optional<Rates> parse_rates(const std::string& text);

std::string dataset_of(const std::string& dir);

// Sorted, so every shard of a run derives the same list and takes positions
// shard, shard + shards, ... of it.
std::vector<std::filesystem::path> trajectory_files(const std::string& dir, int limit = 0,
                                                    int shard = 0, int shards = 1);

// --skip-existing reads a file's existence as "done", so the file appears only once it is
// whole: a killed run leaves a .part behind, never a truncated document.
void write_document(const std::filesystem::path& path, const json::Value& doc);

}  // namespace ssk::pipelines
