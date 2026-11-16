// Run ours at every baseline operating point that ssk_frechet measured.
//
//   ssk_match --in data/trajectories/mopsi --eps 0.1 [--ours convex-shadow]
//             [--distances data/frechet-distances] [--out data/error-matched]
//
// writes, for each trajectory,
//
//   data/error-matched/<ours>/eps0.1/<algorithm>/mopsi/m<rate>/<name>.json   one per baseline run
//   data/error-matched/<ours>/eps0.1/summary/mopsi/<name>.json               the numbers of all
//
// where <ours> is the implementation run, one of pipelines::kOurs, and <algorithm> the
// baseline matched against. --eps is required: it is ours alone, so each value is its own run,
// and a result that does not say which eps produced it cannot be reproduced. Likewise each
// implementation is its own run, under its own root. The summary is written last, so with
// --skip-existing a trajectory that has one is finished, and an interrupted run resumes where
// it stopped.
//
// --validate checks the algorithm's invariants as it runs (src/algo/validate) and puts the
// counts in every result, in the summary and in the closing lines. It is off by default
// because it costs time.
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

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string ours = pipelines::kOurs.front();
  std::string distances = "data/frechet-distances";
  std::string out = "data/error-matched";
  double eps = 0.0;
  double relative_tol = 1e-9;
  int rates = 6;
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
  bool validate = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: ssk_match --in DIR --eps E [--ours NAME] [--distances DIR] [--out DIR] "
               "[--tol REL] [--rates N] [--limit N] [--shard I --shards N] [--skip-existing] "
               "[--validate]\n"
               "  ours:");
  for (const std::string& ours : pipelines::kOurs) {
    std::fprintf(stderr, " %s", ours.c_str());
  }
  std::fprintf(stderr, "\n");
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
    } else if (arg == "--validate") {
      opt.validate = true;
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
  std::size_t matched = 0, points = 0, over = 0, skipped = 0;
  validate::Stats validation{};

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

    std::vector<pipelines::Match> matches;
    for (const std::string& algorithm : pipelines::kBaselines) {
      for (int m = 1; m <= opt.rates; ++m) {
        const std::string rate = "m" + std::to_string(m);
        const fs::path from = fs::path(opt.distances) / algorithm / dataset / rate / name;
        if (!fs::exists(from)) {
          continue;
        }
        const pipelines::Match& match = matches.emplace_back(pipelines::match(
            computer, opt.ours, curve, json::parse_file(from.string()), opt.eps,
            opt.relative_tol, opt.validate));
        pipelines::write_document(root / algorithm / dataset / rate / name,
                       pipelines::match_to_json(input, dataset, source, opt.eps, match));

        points += curve.size();
        over += !match.within_bound;
        if (match.validation) {
          validate::add(validation, *match.validation);
        }
      }
    }
    pipelines::write_document(
        summary,
        pipelines::summary_to_json(input, dataset, source, opt.ours, opt.eps, matches));
    matched += matches.size();
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s %s eps=%g: %zu trajectories -> %zu matches over %zu points in %.1f s\n",
              opt.ours.c_str(), dataset.c_str(), opt.eps, files.size() - skipped, matched, points,
              seconds);
  if (over > 0) {
    std::printf("  %zu over (1+eps)*delta -- a broken guarantee, look at these\n", over);
  }
  if (skipped > 0) {
    std::printf("  %zu already done, skipped\n", skipped);
  }
  if (opt.validate) {
    std::printf("  validation: %zu checks, %zu failed\n", validate::checked(validation),
                validate::failed(validation));
    for (std::size_t i = 0; i < validate::kChecks; ++i) {
      std::printf("    %-18s %zu of %zu failed\n",
                  validate::name(static_cast<validate::Check>(i)), validation[i].failed,
                  validation[i].checked);
    }
  }
  std::printf("  out: %s/%s/%s/{<algorithm>,summary}/%s/\n", opt.out.c_str(), opt.ours.c_str(),
              eps_dir.str().c_str(), dataset.c_str());
  return 0;
}
