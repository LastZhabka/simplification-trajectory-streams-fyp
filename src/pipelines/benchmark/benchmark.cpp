// Time one implementation of ours on one trajectory, and print the timing as one JSON line.
//
//   ssk_bench --in data/trajectories/mopsi/mopsi-000001.json --ours greedy-stabbing-shared
//             --eps 0.1 --delta 0.0004
//
// Only the simplification is timed, not reading the file. scripts/benchmark.py runs it once per
// cell of docs/benchmark.md and kills it at the time limit.
#include <cstdio>
#include <cstdlib>
#include <string>

#include "bench.hpp"
#include "match.hpp"

using namespace ssk;

namespace {

struct Options {
  std::string in;
  std::string ours;
  double eps = 0.0;
  double delta = 0.0;
};

[[noreturn]] void usage() {
  std::fprintf(stderr, "usage: ssk_bench --in FILE --ours NAME --eps E --delta D\n  ours:");
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
    } else if (arg == "--ours" && has_value) {
      opt.ours = argv[++i];
    } else if (arg == "--eps" && has_value) {
      opt.eps = std::atof(argv[++i]);
    } else if (arg == "--delta" && has_value) {
      opt.delta = std::atof(argv[++i]);
    } else {
      usage();
    }
  }
  if (opt.in.empty() || opt.ours.empty() || opt.eps <= 0.0 || opt.delta <= 0.0) {
    usage();
  }
  return opt;
}

}  // namespace

int main(int argc, char** argv) {
  const Options opt = parse_args(argc, argv);
  const io::Document input = io::read_trajectory_file(opt.in);
  const simplify::Curve<2> curve = simplify::curve_of<2>(input);
  const auto ours = pipelines::make_ours(opt.ours, {{"delta", opt.delta}, {"eps", opt.eps}});
  const pipelines::Timing timing = pipelines::time_runs(*ours, curve);
  std::printf("%s\n",
              pipelines::timing_to_json(opt.ours, opt.eps, opt.delta, curve.size(), timing)
                  .dump(0)
                  .c_str());
  return 0;
}
