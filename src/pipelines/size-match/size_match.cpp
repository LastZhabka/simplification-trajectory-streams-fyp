// Ours at the size of every baseline run that ssk_frechet measured -- step 9 of
// docs/running.md.
//
//   ssk_size_match --in data/trajectories/mopsi --eps 0.1 [--ours greedy-stabbing-shared]
//                  [--distances data/frechet-distances] [--out data/size-matched]
//
// writes, for each trajectory, in step 7's layout,
//
//   data/size-matched/<ours>/eps0.1/<algorithm>/mopsi/m<rate>/<name>.json   one per baseline run
//   data/size-matched/<ours>/eps0.1/summary/mopsi/<name>.json               the numbers of all
//
// The summary is written last, so with --skip-existing an interrupted run resumes.
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
#include "size.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string ours = pipelines::kOurs.front();
  std::string distances = "data/frechet-distances";
  std::string out = "data/size-matched";
  double eps = 0.0;
  pipelines::Rates rates;
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: ssk_size_match --in DIR --eps E [--ours NAME] [--distances DIR] "
               "[--out DIR] [--rates A..B] [--limit N] [--shard I --shards N] "
               "[--skip-existing]\n");
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
    } else if (arg == "--distances" && has_value) {
      opt.distances = argv[++i];
    } else if (arg == "--out" && has_value) {
      opt.out = argv[++i];
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
  std::size_t done = 0, matched = 0, runs = 0, over = 0, skipped = 0;

  for (const fs::path& file : files) {
    const std::string name = file.filename().string();
    const fs::path summary = root / "summary" / dataset / name;
    if (opt.skip_existing && fs::exists(summary)) {
      ++skipped;
      continue;
    }
    const io::Document input = io::read_trajectory_file(file.string());
    const simplify::Curve<2> curve = simplify::curve_of<2>(input);
    const std::string source = dataset + "/" + name;

    std::vector<pipelines::Fitted> fitted;
    std::vector<fs::path> where;
    std::vector<std::size_t> budgets;
    for (const std::string& algorithm : pipelines::kBaselines) {
      for (int m = opt.rates.first; m <= opt.rates.last; ++m) {
        const fs::path at = fs::path(algorithm) / dataset / ("m" + std::to_string(m)) / name;
        const fs::path from = fs::path(opt.distances) / at;
        if (!fs::exists(from)) {
          continue;
        }
        const json::Value& baseline =
            fitted.emplace_back(pipelines::Fitted{json::parse_file(from.string()), {}}).baseline;
        budgets.push_back(pipelines::budget_of(baseline));
        where.push_back(at);
      }
    }
    const std::vector<pipelines::Fit> fits =
        pipelines::fit(computer, opt.ours, curve, opt.eps, budgets);
    for (std::size_t i = 0; i < fitted.size(); ++i) {
      fitted[i].fit = fits[i];
      pipelines::write_document(root / where[i],
                                pipelines::fit_to_json(input, dataset, source, opt.ours, opt.eps,
                                                       fits[i], fitted[i].baseline));
      runs += static_cast<std::size_t>(fits[i].runs);
      over += !fits[i].within_bound;
    }
    pipelines::write_document(summary, pipelines::fit_summary_to_json(input, dataset, source,
                                                                      opt.ours, opt.eps, fitted));
    matched += fitted.size();
    ++done;
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s %s eps=%g: %zu trajectories -> %zu fits, %zu runs of ours, in %.1f s\n",
              opt.ours.c_str(), dataset.c_str(), opt.eps, done, matched, runs, seconds);
  if (over > 0) {
    std::printf("  %zu over (1+eps)*delta -- a broken guarantee, look at these\n", over);
  }
  if (skipped > 0) {
    std::printf("  %zu already done, skipped\n", skipped);
  }
  std::printf("  out: %s/%s/%s/{<algorithm>,summary}/%s/\n", opt.out.c_str(), opt.ours.c_str(),
              eps_dir.str().c_str(), dataset.c_str());
  return 0;
}
