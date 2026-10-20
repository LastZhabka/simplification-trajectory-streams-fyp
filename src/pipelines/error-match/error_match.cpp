// Run ours at every baseline operating point that ssk_frechet measured.
//
//   ssk_match --in data/trajectories/mopsi --eps 0.1
//             [--distances data/frechet-distances] [--out data/error-matched]
//
// writes, for each trajectory,
//
//   data/error-matched/eps0.1/<algorithm>/mopsi/m<rate>/<name>.json   ours, one per baseline run
//   data/error-matched/eps0.1/summary/mopsi/<name>.json               the numbers of all of them
//
// where <algorithm> is the baseline matched against. --eps is required: it is ours alone, so
// each value is its own run, and a result that does not say which eps produced it cannot be
// reproduced. The summary is written last, so with --skip-existing a trajectory that has one
// is finished, and an interrupted run resumes where it stopped.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "match.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string distances = "data/frechet-distances";
  std::string out = "data/error-matched";
  double eps = 0.0;
  double relative_tol = 1e-9;
  int rates = 6;
  int limit = 0;
  int shard = 0;
  int shards = 1;
  bool skip_existing = false;
};

void usage() {
  std::fprintf(stderr,
               "usage: ssk_match --in DIR --eps E [--distances DIR] [--out DIR] [--tol REL] "
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
    } else if (arg == "--eps" && has_value) {
      opt.eps = std::atof(argv[++i]);
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
    } else {
      usage();
    }
  }
  if (opt.in.empty() || opt.eps <= 0.0) {
    usage();
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

  std::ostringstream eps_dir;
  eps_dir << "eps" << opt.eps;
  const fs::path root = fs::path(opt.out) / eps_dir.str();
  const std::vector<std::string> algorithms{"douglas-peucker-n", "squish", "dots"};
  const frechet::Frechet<2> computer;
  const auto started = std::chrono::steady_clock::now();
  std::size_t matched = 0, points = 0, over = 0, skipped = 0;

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
    for (const std::string& algorithm : algorithms) {
      for (int m = 1; m <= opt.rates; ++m) {
        const std::string rate = "m" + std::to_string(m);
        const fs::path from = fs::path(opt.distances) / algorithm / dataset / rate / name;
        if (!fs::exists(from)) {
          continue;
        }
        const pipelines::Match& match = matches.emplace_back(pipelines::match(
            computer, curve, json::parse_file(from.string()), opt.eps, opt.relative_tol));
        write_document(root / algorithm / dataset / rate / name,
                       pipelines::match_to_json(input, dataset, source, opt.eps, match));

        points += curve.size();
        over += !match.within_bound;
      }
    }
    write_document(summary,
                   pipelines::summary_to_json(input, dataset, source, opt.eps, matches));
    matched += matches.size();
  }

  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
  std::printf("%s eps=%g: %zu trajectories -> %zu matches over %zu points in %.1f s\n",
              dataset.c_str(), opt.eps, files.size() - skipped, matched, points, seconds);
  if (over > 0) {
    std::printf("  %zu over (1+eps)*delta -- a broken guarantee, look at these\n", over);
  }
  if (skipped > 0) {
    std::printf("  %zu already done, skipped\n", skipped);
  }
  std::printf("  out: %s/%s/{<algorithm>,summary}/%s/\n", opt.out.c_str(), eps_dir.str().c_str(),
              dataset.c_str());
  return 0;
}
