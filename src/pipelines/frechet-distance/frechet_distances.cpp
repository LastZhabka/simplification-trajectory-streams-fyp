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
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "corpus.hpp"
#include "measure.hpp"
#include "simplify/simplifier.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string simplified = "data/simplified-trajectories";
  std::string out = "data/frechet-distances";
  double relative_tol = 1e-9;
  pipelines::Rates rates;
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: ssk_frechet --in DIR [--simplified DIR] [--out DIR] [--tol REL] "
               "[--rates A..B] [--limit N] [--shard I --shards N] [--skip-existing]\n");
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
      const auto rates = pipelines::parse_rates(argv[++i]);
      if (!rates) {
        usage();
      }
      opt.rates = *rates;
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

}  // namespace

int main(int argc, char** argv) {
  const Options opt = parse_args(argc, argv);
  const std::string dataset = pipelines::dataset_of(opt.in);
  const std::vector<fs::path> files =
      pipelines::trajectory_files(opt.in, opt.limit, opt.shard, opt.shards);

  const frechet::Frechet<2> computer;
  const auto started = std::chrono::steady_clock::now();
  std::size_t written = 0, cells = 0, skipped = 0;
  double worst = 0.0;

  for (const fs::path& file : files) {
    const io::Document input = io::read_trajectory_file(file.string());
    const frechet::Curve<2> curve = simplify::curve_of<2>(input);
    const std::string name = file.filename().string();

    for (const std::string& algorithm : pipelines::kBaselines) {
      for (int m = opt.rates.first; m <= opt.rates.last; ++m) {
        const std::string rate = "m" + std::to_string(m);
        const fs::path from = fs::path(opt.simplified) / algorithm / dataset / rate / name;
        if (!fs::exists(from)) {
          continue;
        }
        const fs::path to = fs::path(opt.out) / algorithm / dataset / rate / name;
        if (opt.skip_existing && fs::exists(to)) {
          ++skipped;
          continue;
        }
        const json::Value doc = json::parse_file(from.string());
        const frechet::Curve<2> simplified = simplify::curve_of<2>(io::read_trajectory(doc));

        const pipelines::Measurement measured =
            pipelines::measure(computer, curve, simplified, opt.relative_tol);
        cells += curve.size() * simplified.size();
        worst = std::max(worst, measured.distance);

        pipelines::write_document(
            to, pipelines::to_json(input, doc, dataset, dataset + "/" + name,
                                   algorithm + "/" + dataset + "/" + rate + "/" + name, measured));
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
