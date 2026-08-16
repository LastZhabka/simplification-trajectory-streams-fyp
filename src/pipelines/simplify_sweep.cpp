// Simplify every trajectory in a dataset directory with all three baselines, at a fixed
// compression rate of 1/2^m for m = 1..rates.
//
//   ssk_simplify --in data/trajectories/mopsi --out data/simplified-trajectories
//
// writes data/simplified-trajectories/<algorithm>/mopsi/m<rate>/<name>.json, one document per
// trajectory per algorithm per rate. Each is a trajectory document -- points and their
// timestamps -- with the result fields on top, so it reads and draws like any other. SQUISH
// and DOTS need timestamps, so a document without a clock is simplified by DPn alone.
//
// One dataset can be split across processes with `--shard I --shards N`, which is worth
// doing on GeoLife: DOTS' cost grows faster than linearly, and a handful of 60k-point tracks
// there dominate the run. Shards stride through the file list rather than take blocks, so
// those long tracks spread across workers instead of landing on one.
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "io/trajectory.hpp"
#include "sweep.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string out = "data/simplified-trajectories";
  int rates = 6;
  int limit = 0;   // 0 means every trajectory
  int shard = 0;   // this worker's index, for splitting one dataset across processes
  int shards = 1;
};

Options parse_args(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const bool has_value = i + 1 < argc;
    if (arg == "--in" && has_value) {
      opt.in = argv[++i];
    } else if (arg == "--out" && has_value) {
      opt.out = argv[++i];
    } else if (arg == "--rates" && has_value) {
      opt.rates = std::atoi(argv[++i]);
    } else if (arg == "--limit" && has_value) {
      opt.limit = std::atoi(argv[++i]);
    } else if (arg == "--shard" && has_value) {
      opt.shard = std::atoi(argv[++i]);
    } else if (arg == "--shards" && has_value) {
      opt.shards = std::atoi(argv[++i]);
    } else {
      std::fprintf(stderr, "usage: ssk_simplify --in DIR [--out DIR] [--rates N] [--limit N] [--shard I --shards N]\n");
      std::exit(2);
    }
  }
  if (opt.in.empty()) {
    std::fprintf(stderr, "usage: ssk_simplify --in DIR [--out DIR] [--rates N] [--limit N] [--shard I --shards N]\n");
    std::exit(2);
  }
  return opt;
}

void write_document(const fs::path& path, const json::Value& doc) {
  fs::create_directories(path.parent_path());
  std::ofstream out(path);
  doc.write(out);
  out << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  const Options opt = parse_args(argc, argv);
  // A trailing separator -- which tab completion adds -- makes filename() empty.
  const fs::path in_dir = fs::path(opt.in).lexically_normal();
  const std::string dataset = in_dir.has_filename()
                                  ? in_dir.filename().string()
                                  : in_dir.parent_path().filename().string();

  std::vector<fs::path> files;
  for (const auto& entry : fs::directory_iterator(in_dir)) {
    if (entry.path().extension() == ".json" && entry.path().filename() != "index.json") {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  if (opt.limit > 0 && files.size() > static_cast<std::size_t>(opt.limit)) {
    files.resize(static_cast<std::size_t>(opt.limit));
  }
  if (opt.shards > 1) {
    std::vector<fs::path> mine;
    for (std::size_t i = static_cast<std::size_t>(opt.shard); i < files.size();
         i += static_cast<std::size_t>(opt.shards)) {
      mine.push_back(files[i]);
    }
    files = std::move(mine);
  }

  const auto started = std::chrono::steady_clock::now();
  std::size_t points = 0, written = 0, no_clock = 0;
  for (const fs::path& file : files) {
    const io::Document doc = io::read_trajectory_file(file.string());
    const auto curve = simplify::curve_of<2>(doc);
    const simplify::Context in = simplify::context_of(doc);
    points += curve.size();

    std::vector<pipelines::Sweep> sweeps{pipelines::sweep_dpn(curve, opt.rates)};
    if (in.t.size() == curve.size()) {
      sweeps.push_back(pipelines::sweep_squish(curve, in, opt.rates));
      sweeps.push_back(pipelines::sweep_dots(curve, in, opt.rates));
    } else {
      ++no_clock;
    }

    const std::string source = dataset + "/" + file.filename().string();
    for (const pipelines::Sweep& sweep : sweeps) {
      for (const pipelines::Run& run : sweep.runs) {
        const std::string rate = "m" + std::to_string(run.m);
        write_document(fs::path(opt.out) / sweep.algorithm / dataset / rate / file.filename(),
                       pipelines::run_to_json(doc, source, sweep.algorithm, run));
        ++written;
      }
    }
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s: %zu trajectories, %zu points -> %zu documents in %.1f s\n", dataset.c_str(),
              files.size(), points, written, seconds);
  if (no_clock > 0) {
    std::printf("  %zu without timestamps, simplified by DPn only\n", no_clock);
  }
  std::printf("  out: %s/<algorithm>/%s/\n", opt.out.c_str(), dataset.c_str());
  return 0;
}
