// Check every result document in a corpus against the trajectory it came from.
//
//   ssk_audit --in data/trajectories/mopsi --results data/simplified-trajectories
//
// Prints one line per offender and a summary. `--paths` prints only the offending paths, so
// removing them is a pipe:
//
//   ssk_audit --in ... --results ... --paths | xargs rm          # POSIX
//   ssk_audit --in ... --results ... --paths | Remove-Item       # PowerShell
//
// The tool never deletes anything itself: verification and mutation are worth keeping apart.
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "audit.hpp"

namespace fs = std::filesystem;
using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string results = "data/simplified-trajectories";
  int rates = 6;
  bool paths_only = false;
};

Options parse_args(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const bool has_value = i + 1 < argc;
    if (arg == "--in" && has_value) {
      opt.in = argv[++i];
    } else if (arg == "--results" && has_value) {
      opt.results = argv[++i];
    } else if (arg == "--rates" && has_value) {
      opt.rates = std::atoi(argv[++i]);
    } else if (arg == "--paths") {
      opt.paths_only = true;
    } else {
      std::fprintf(stderr,
                   "usage: ssk_audit --in DIR [--results DIR] [--rates N] [--paths]\n");
      std::exit(2);
    }
  }
  if (opt.in.empty()) {
    std::fprintf(stderr, "usage: ssk_audit --in DIR [--results DIR] [--rates N] [--paths]\n");
    std::exit(2);
  }
  return opt;
}

}  // namespace

int main(int argc, char** argv) {
  const Options opt = parse_args(argc, argv);
  const fs::path in_dir = fs::path(opt.in).lexically_normal();
  const std::string dataset = in_dir.has_filename()
                                  ? in_dir.filename().string()
                                  : in_dir.parent_path().filename().string();
  const std::vector<std::string> algorithms{"douglas-peucker-n", "squish", "dots"};

  std::size_t checked = 0, bad = 0;
  for (const auto& entry : fs::directory_iterator(in_dir)) {
    if (entry.path().extension() != ".json" || entry.path().filename() == "index.json") {
      continue;
    }
    const io::Document input = io::read_trajectory_file(entry.path().string());
    const std::string name = entry.path().filename().string();

    for (const std::string& algorithm : algorithms) {
      for (int m = 1; m <= opt.rates; ++m) {
        const fs::path p =
            fs::path(opt.results) / algorithm / dataset / ("m" + std::to_string(m)) / name;
        if (!fs::exists(p)) {
          continue;
        }
        ++checked;
        if (pipelines::is_subsequence(input, io::read_trajectory_file(p.string()))) {
          continue;
        }
        ++bad;
        if (opt.paths_only) {
          std::printf("%s\n", p.string().c_str());
        } else {
          std::printf("BAD %s %s %s m%d\n", algorithm.c_str(), dataset.c_str(), name.c_str(),
                      m);
        }
        std::fflush(stdout);
      }
    }
  }

  if (!opt.paths_only) {
    std::printf("%s: %zu documents checked, %zu not a subsequence of their input\n",
                dataset.c_str(), checked, bad);
  }
  return bad == 0 ? 0 : 1;
}
