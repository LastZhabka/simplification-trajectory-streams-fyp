// Measure every simplification against the trajectory it came from.
//
//   ssk_frechet --in data/trajectories/mopsi
//               --simplified data/simplified-trajectories --out data/frechet-distances
//
// writes data/frechet-distances/<algorithm>/mopsi/m<rate>/<name>.json, mirroring the
// simplified tree one for one: one document per simplification, carrying its Frechet distance
// plus the algorithm, mode, params and stats that produced it.
//
// Cost is quadratic in trajectory length at a fixed rate -- the free-space diagram is n by m
// cells -- so the long GeoLife tracks dominate, and `--shard I --shards N` is here for the
// same reason it is in the sweep. A whole-corpus pass runs for hours, so `--skip-existing`
// makes an interrupted one resumable: a measurement is a pure function of the two documents
// it reads, so one already written is one that need not be computed again.
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "measure.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string simplified = "data/simplified-trajectories";
  std::string out = "data/frechet-distances";
  double relative_tol = 1e-9;
  int rates = 6;
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: ssk_frechet --in DIR [--simplified DIR] [--out DIR] [--tol REL] "
               "[--rates N] [--limit N] [--shard I --shards N] [--skip-existing]\n");
  std::exit(2);
}

Options parse_args(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const bool has_value = i + 1 < argc;
    if (arg == "--in" && has_value) {
      opt.in = argv[++i];
    } else if (arg == "--simplified" && has_value) {
      opt.simplified = argv[++i];
    } else if (arg == "--out" && has_value) {
      opt.out = argv[++i];
    } else if (arg == "--tol" && has_value) {
      opt.relative_tol = std::atof(argv[++i]);
    } else if (arg == "--rates" && has_value) {
      opt.rates = std::atoi(argv[++i]);
    } else if (arg == "--limit" && has_value) {
      opt.limit = std::atoi(argv[++i]);
    } else if (arg == "--shard" && has_value) {
      opt.shard = std::atoi(argv[++i]);
    } else if (arg == "--shards" && has_value) {
      opt.shards = std::atoi(argv[++i]);
    } else if (arg == "--skip-existing") {
      opt.skip_existing = true;
    } else {
      usage();
    }
  }
  if (opt.in.empty()) {
    usage();
  }
  return opt;
}

frechet::Curve<2> curve_of(const io::Document& doc) {
  frechet::Curve<2> c;
  c.reserve(doc.points.size());
  for (const io::Point& p : doc.points) {
    c.push_back({p[0], p[1]});
  }
  return c;
}

}  // namespace

int main(int argc, char** argv) {
  const Options opt = parse_args(argc, argv);
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

  const std::vector<std::string> algorithms{"douglas-peucker-n", "squish", "dots"};
  const frechet::Frechet<2> computer;
  const auto started = std::chrono::steady_clock::now();
  std::size_t written = 0, cells = 0, skipped = 0;
  double worst = 0.0;

  for (const fs::path& file : files) {
    const io::Document input = io::read_trajectory_file(file.string());
    const frechet::Curve<2> curve = curve_of(input);
    const std::string name = file.filename().string();

    for (const std::string& algorithm : algorithms) {
      for (int m = 1; m <= opt.rates; ++m) {
        const std::string rate = "m" + std::to_string(m);
        const fs::path from = fs::path(opt.simplified) / algorithm / dataset / rate / name;
        if (!fs::exists(from)) {
          continue;
        }
        const fs::path to = fs::path(opt.out) / algorithm / dataset / rate / name;
        // Resuming an interrupted run: the measurement is a pure function of the two
        // documents, so one already on disk is one that need not be computed again.
        if (opt.skip_existing && fs::exists(to)) {
          ++skipped;
          continue;
        }
        const json::Value doc = json::parse_file(from.string());
        const frechet::Curve<2> simplified = curve_of(io::read_trajectory(doc));

        const pipelines::Measurement measured =
            pipelines::measure(computer, curve, simplified, opt.relative_tol);
        cells += curve.size() * simplified.size();
        worst = std::max(worst, measured.distance);

        fs::create_directories(to.parent_path());
        std::ofstream stream(to);
        pipelines::to_json(input, doc, dataset, dataset + "/" + name,
                           algorithm + "/" + dataset + "/" + rate + "/" + name, measured)
            .write(stream);
        stream << '\n';
        ++written;
      }
    }
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s: %zu trajectories -> %zu measurements in %.1f s (%.2f G cells, worst %.6g)\n",
              dataset.c_str(), files.size(), written, seconds,
              static_cast<double>(cells) / 1e9, worst);
  if (skipped > 0) {
    std::printf("  %zu already measured, skipped\n", skipped);
  }
  std::printf("  out: %s/<algorithm>/%s/\n", opt.out.c_str(), dataset.c_str());
  return 0;
}
