#include "corpus.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>

namespace fs = std::filesystem;

namespace ssk::pipelines {

std::optional<Rates> parse_rates(const std::string& text) {
  Rates rates;
  char rest = 0;
  if (std::sscanf(text.c_str(), "%d..%d%c", &rates.first, &rates.last, &rest) != 2 ||
      rates.first < 1 || rates.first > rates.last) {
    return std::nullopt;
  }
  return rates;
}

std::string dataset_of(const std::string& dir) {
  // A trailing separator -- which tab completion adds -- makes filename() empty.
  const fs::path path = fs::path(dir).lexically_normal();
  return path.has_filename() ? path.filename().string() : path.parent_path().filename().string();
}

std::vector<fs::path> trajectory_files(const std::string& dir, int limit, int shard, int shards) {
  std::vector<fs::path> files;
  for (const auto& entry : fs::directory_iterator(dir)) {
    if (entry.path().extension() == ".json" && entry.path().filename() != "index.json") {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  if (limit > 0 && files.size() > static_cast<std::size_t>(limit)) {
    files.resize(static_cast<std::size_t>(limit));
  }
  std::vector<fs::path> mine;
  for (auto i = static_cast<std::size_t>(shard); i < files.size();
       i += static_cast<std::size_t>(shards)) {
    mine.push_back(files[i]);
  }
  return mine;
}

void write_document(const fs::path& path, const json::Value& doc) {
  fs::create_directories(path.parent_path());
  fs::path part = path;
  part += ".part";
  {
    std::ofstream out(part);
    doc.write(out);
    out << '\n';
  }
  fs::rename(part, path);
}

}  // namespace ssk::pipelines
