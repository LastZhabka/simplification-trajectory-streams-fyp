// Ours at every error target ssk_target_search set, at delta = T / (1 + eps), beside the baselines.
//
//   ssk_target_match --in data/trajectories/mopsi --eps 0.1 [--ours NAME]
//                    [--targets data/target-simplified] [--out data/target-matched]
//                    [--limit N] [--shard I --shards N] [--skip-existing]
//
// writes, for each trajectory,
//
//   data/target-matched/<ours>/eps0.1/results/mopsi/f<fraction>/<name>.json   ours at a target
//   data/target-matched/<ours>/eps0.1/summary/mopsi/<name>.json               beside the baselines
//
// The summary is written last, so with --skip-existing an interrupted run resumes.
//
// --unscaled-delta runs ours at delta = T: within (1 + eps) T, so with more error than the
// baselines, and with its size bound against kappa(T). It is step 10 of docs/running.md, and
// writes under data/target-matched-unscaled unless --out says else.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "corpus.hpp"
#include "match.hpp"
#include "target.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string ours = pipelines::kOurs.front();
  std::string targets = "data/target-simplified";
  std::string out;
  double eps = 0.0;
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
  bool unscaled = false;
};

void usage() {
  std::fprintf(stderr, "usage: ssk_target_match --in DIR --eps E [--ours NAME] [--targets DIR] "
                       "[--out DIR] [--limit N] [--shard I --shards N] [--skip-existing] "
                       "[--unscaled-delta]\n");
  std::exit(2);
}

Options parse_args(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const bool has_value = i + 1 < argc;
    if (arg == "--in" && has_value) {
      opt.in = argv[++i];
    } else if (arg == "--eps" && has_value) {
      opt.eps = std::atof(argv[++i]);
    } else if (arg == "--ours" && has_value) {
      opt.ours = argv[++i];
    } else if (arg == "--targets" && has_value) {
      opt.targets = argv[++i];
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
    } else if (arg == "--unscaled-delta") {
      opt.unscaled = true;
    } else {
      usage();
    }
  }
  if (opt.out.empty()) {
    opt.out = opt.unscaled ? "data/target-matched-unscaled" : "data/target-matched";
  }
  if (opt.in.empty() || opt.eps <= 0.0 ||
      std::find(pipelines::kOurs.begin(), pipelines::kOurs.end(), opt.ours) ==
          pipelines::kOurs.end()) {
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
  std::ostringstream eps_dir;
  eps_dir << "eps" << opt.eps;
  const fs::path root = fs::path(opt.out) / opt.ours / eps_dir.str();
  const frechet::Frechet<2> computer;

  const auto started = std::chrono::steady_clock::now();
  std::size_t done = 0, skipped = 0, missing = 0, runs = 0, over = 0, over_target = 0;
  for (const fs::path& file : files) {
    const std::string name = file.filename().string();
    const fs::path summary = root / "summary" / dataset / name;
    if (opt.skip_existing && fs::exists(summary)) {
      ++skipped;
      continue;
    }
    const fs::path searched = fs::path(opt.targets) / "summary" / dataset / name;
    if (!fs::exists(searched)) {
      ++missing;
      continue;
    }
    const json::Value search_summary = json::parse_file(searched.string());
    const io::Document doc = io::read_trajectory_file(file.string());
    const auto curve = simplify::curve_of<2>(doc);
    const std::string source = dataset + "/" + name;

    std::vector<pipelines::OursRun> matched;
    for (const pipelines::Target& target : pipelines::targets_of(search_summary)) {
      const pipelines::OursRun& run =
          matched.emplace_back(pipelines::run_ours(computer, opt.ours, curve, target, opt.eps,
                                                   opt.unscaled));
      pipelines::write_document(root / "results" / dataset / target.level / name,
                                pipelines::ours_to_json(doc, dataset, source, opt.ours, opt.eps,
                                                        opt.unscaled, run));
      ++runs;
      over += !run.within;
      over_target += !run.within_target;
    }
    pipelines::write_document(
        summary, pipelines::match_summary_to_json(doc, dataset, source, opt.ours, opt.eps,
                                                  opt.unscaled, search_summary, matched));
    ++done;
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s %s eps=%g: %zu trajectories -> %zu runs in %.1f s\n", opt.ours.c_str(),
              dataset.c_str(), opt.eps, done, runs, seconds);
  if (over > 0) {
    std::printf("  %zu over (1+eps)*delta -- a broken guarantee, look at these\n", over);
  }
  if (over_target > 0) {
    std::printf("  %zu over the target%s\n", over_target,
                opt.unscaled ? ", which delta = T allows" : " -- a broken guarantee");
  }
  if (missing > 0) {
    std::printf("  %zu with no ssk_target_search summary, skipped\n", missing);
  }
  if (skipped > 0) {
    std::printf("  %zu already done, skipped\n", skipped);
  }
  std::printf("  out: %s/%s/%s/{results,summary}/%s/\n", opt.out.c_str(), opt.ours.c_str(),
              eps_dir.str().c_str(), dataset.c_str());
  return 0;
}
