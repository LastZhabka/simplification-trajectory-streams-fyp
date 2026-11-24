// Every baseline at every error target of every trajectory in a dataset directory: the fewest
// vertices each was found to reach within T = f * diagonal, f in pipelines::kFractions.
//
//   ssk_target_search --in data/trajectories/mopsi [--out data/target-simplified]
//                     [--limit N] [--shard I --shards N] [--skip-existing] [--no-dots]
//
// writes, for each trajectory,
//
//   data/target-simplified/<algorithm>/mopsi/f<fraction>/<name>.json   the simplification
//   data/target-simplified/summary/mopsi/<name>.json                    every result, written last
//
// SQUISH, SQUISH-E and DOTS need timestamps, so a document without a clock gets DPn alone.
// --no-dots is for the trajectories over 30 000 points, on which DOTS' search takes minutes to
// hours; docs/running.md step 8.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "corpus.hpp"
#include "io/trajectory.hpp"
#include "measure.hpp"
#include "target.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string out = "data/target-simplified";
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
  bool dots = true;
};

void usage() {
  std::fprintf(stderr, "usage: ssk_target_search --in DIR [--out DIR] [--limit N] "
                       "[--shard I --shards N] [--skip-existing] [--no-dots]\n");
  std::exit(2);
}

Options parse_args(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const bool has_value = i + 1 < argc;
    if (arg == "--in" && has_value) {
      opt.in = argv[++i];
    } else if (arg == "--out" && has_value) {
      opt.out = argv[++i];
    } else if (arg == "--limit" && has_value) {
      opt.limit = std::atoi(argv[++i]);
    } else if (arg == "--shard" && has_value) {
      opt.shard = std::atoi(argv[++i]);
    } else if (arg == "--shards" && has_value) {
      opt.shards = std::atoi(argv[++i]);
    } else if (arg == "--skip-existing") {
      opt.skip_existing = true;
    } else if (arg == "--no-dots") {
      opt.dots = false;
    } else {
      usage();
    }
  }
  if (opt.in.empty()) {
    usage();
  }
  return opt;
}

bool increasing(const std::vector<std::size_t>& indices) {
  for (std::size_t i = 1; i < indices.size(); ++i) {
    if (indices[i - 1] >= indices[i]) {
      return false;
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const Options opt = parse_args(argc, argv);
  const std::string dataset = pipelines::dataset_of(opt.in);
  const std::vector<fs::path> files =
      pipelines::trajectory_files(opt.in, opt.limit, opt.shard, opt.shards);
  const frechet::Frechet<2> computer;

  const auto started = std::chrono::steady_clock::now();
  std::size_t done = 0, skipped = 0, no_extent = 0, written = 0, unordered = 0;
  long runs = 0;
  for (const fs::path& file : files) {
    const std::string name = file.filename().string();
    const fs::path summary = fs::path(opt.out) / "summary" / dataset / name;
    if (opt.skip_existing && fs::exists(summary)) {
      ++skipped;
      continue;
    }
    const io::Document doc = io::read_trajectory_file(file.string());
    const auto curve = simplify::curve_of<2>(doc);
    const simplify::Context in = simplify::context_of(doc);
    const std::string source = dataset + "/" + name;
    const std::vector<pipelines::Target> targets =
        curve.size() < 3 ? std::vector<pipelines::Target>{}
                         : pipelines::targets(pipelines::diagonal(curve));
    no_extent += targets.empty();

    std::vector<pipelines::Found> found;
    if (!targets.empty()) {
      const bool clock = in.t.size() == curve.size();
      for (const std::string& algorithm : pipelines::kTargetBaselines) {
        if ((algorithm != "douglas-peucker-n" && !clock) || (algorithm == "dots" && !opt.dots)) {
          continue;
        }
        for (pipelines::Best& best :
             pipelines::search(algorithm, computer, curve, in, targets)) {
          runs += best.algorithm_runs;
          if (!increasing(best.indices)) {
            ++unordered;
            std::fprintf(stderr, "dropped %s %s %s: indices not increasing\n", algorithm.c_str(),
                         source.c_str(), best.target.level.c_str());
            continue;
          }
          simplify::Curve<2> kept;
          for (const std::size_t i : best.indices) {
            kept.push_back(curve[i]);
          }
          const pipelines::Measurement measured = pipelines::measure(computer, curve, kept, 1e-9);
          pipelines::write_document(
              fs::path(opt.out) / algorithm / dataset / best.target.level / name,
              pipelines::best_to_json(doc, dataset, source, algorithm, best, measured));
          ++written;
          found.push_back({algorithm, std::move(best), measured});
        }
      }
    }
    pipelines::write_document(
        summary, pipelines::search_summary_to_json(doc, dataset, source, targets, found));
    ++done;
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s: %zu trajectories -> %zu documents, %ld algorithm runs, in %.1f s\n",
              dataset.c_str(), done, written, runs, seconds);
  if (no_extent > 0) {
    std::printf("  %zu with no extent or under 3 points, no targets\n", no_extent);
  }
  if (skipped > 0) {
    std::printf("  %zu already done, skipped\n", skipped);
  }
  if (unordered > 0) {
    std::printf("  %zu results dropped: indices not increasing, so not a subsequence\n",
                unordered);
  }
  std::printf("  out: %s/{<algorithm>,summary}/%s/\n", opt.out.c_str(), dataset.c_str());
  return 0;
}
