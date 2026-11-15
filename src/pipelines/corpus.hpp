#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "io/json.hpp"

namespace ssk::pipelines {

inline const std::vector<std::string> kBaselines{"douglas-peucker-n", "squish", "squish-e",
                                                 "dots"};

std::string dataset_of(const std::string& dir);

// Sorted, so every shard of a run derives the same list and takes positions
// shard, shard + shards, ... of it.
std::vector<std::filesystem::path> trajectory_files(const std::string& dir, int limit = 0,
                                                    int shard = 0, int shards = 1);

// --skip-existing reads a file's existence as "done", so the file appears only once it is
// whole: a killed run leaves a .part behind, never a truncated document.
void write_document(const std::filesystem::path& path, const json::Value& doc);

}  // namespace ssk::pipelines
