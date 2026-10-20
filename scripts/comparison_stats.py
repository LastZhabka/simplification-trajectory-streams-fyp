"""The statistics docs/comparison.md asks for, gathered from what ssk_match wrote.

    python scripts/comparison_stats.py [--matched data/error-matched] [--out FILE] [--jobs N]

Reads every <matched>/eps<E>/summary/<dataset>/<name>.json -- one per trajectory per eps, each
baseline operating point beside our run at its error -- and prints, per dataset:

  coverage    trajectories in each (rate, algorithm) cell, and our matches beside them
  fixed rate  each baseline's measured error at a budget of N/2^m: the headline table
  matched     ours at each baseline's error, per eps: size and error against the baseline's,
              and Pareto dominance in both directions

then the same matched table pooled over every dataset, and the paper's Theorem 3 checked: error
within (1+eps)*delta, and size within 2*kappa - 2 wherever a curve already on disk bounds kappa.
All of it, plus the rate-distortion curves, is also written as JSON to --out.
"""

import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from statistics import mean, median

ALGORITHMS = ("douglas-peucker-n", "squish", "dots")
SHORT = {"douglas-peucker-n": "DPn", "squish": "SQUISH", "dots": "DOTS"}
NOTES = (
    "DOTS reaches each budget by a search on its threshold: its rows describe the "
    "simplification, not an online algorithm.",
    "Ours places vertices that are not input points; every baseline keeps a subsequence.",
    "Errors are in each dataset's own coordinate units; only ratios pool across datasets.",
    "Everything is scored under the Frechet distance; SQUISH and DOTS minimise SED, which is "
    "not scored here.",
)


def read(path):
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    return [{"dataset": doc["dataset"], "source": doc["source"], "n": doc["input_size"],
             "eps": doc["eps"], **pair} for pair in doc["pairs"]]


def load(root, jobs):
    # One small file per trajectory per eps, so reading is the cost; spread it over processes.
    paths = sorted(glob.glob(os.path.join(root, "eps*", "summary", "*", "*.json")))
    with ProcessPoolExecutor(jobs) as pool:
        return [pair for pairs in pool.map(read, paths, chunksize=256) for pair in pairs]


def group(rows, *keys):
    out = defaultdict(list)
    for row in rows:
        out[tuple(row[k] for k in keys)].append(row)
    return out


def fixed_rate(rows):
    return {
        "trajectories": len(rows),
        "error_median": median(r["baseline"]["distance"] for r in rows),
        "error_mean": mean(r["baseline"]["distance"] for r in rows),
        "compression_median": median(r["n"] / r["baseline"]["output_size"] for r in rows),
        "exact_budget": mean(r["baseline"]["output_size"] == r["baseline"]["target"] for r in rows),
        "runs_mean": mean(r["baseline"]["algorithm_runs"] for r in rows),
    }


def matched(rows):
    ko = [r["ours"]["output_size"] for r in rows]
    kb = [r["baseline"]["output_size"] for r in rows]
    do = [r["ours"]["distance"] for r in rows]
    db = [r["baseline"]["distance"] for r in rows]
    pairs = list(zip(ko, kb, do, db))
    return {
        "pairs": len(rows),
        "size_ratio_median": median(o / b for o, b in zip(ko, kb)),
        "error_ratio_median": median(o / b for o, b in zip(do, db)),
        "smaller": mean(o < b for o, b in zip(ko, kb)),
        "equal": mean(o == b for o, b in zip(ko, kb)),
        "larger": mean(o > b for o, b in zip(ko, kb)),
        "ours_dominates": mean(a < b and c <= d for a, b, c, d in pairs),
        "baseline_dominates": mean(b < a and d <= c for a, b, c, d in pairs),
    }


def curve(rows, size, error):
    return [{"m": m, "pairs": len(cell),
             "compression": median(r["n"] / size(r) for r in cell),
             "error": median(error(r) for r in cell)}
            for (m,), cell in sorted(group(rows, "m").items())]


def bounds(pairs, base):
    # Every distance on disk is the top of a bisection, so its curve is certainly within it:
    # any such curve bounds kappa(tau, delta) from above, and so does the input itself.
    known = defaultdict(list)
    for b in base:
        known[b["source"]].append((b["baseline"]["distance"], b["baseline"]["output_size"]))
    for p in pairs:
        known[p["source"]].append((p["ours"]["distance"], p["ours"]["output_size"]))

    out = {}
    for (eps,), rows in sorted(group(pairs, "eps").items()):
        error_over = size_over = tighter = 0
        worst_error = worst_size = 0.0
        for r in rows:
            ours = r["ours"]
            bound = (1 + eps) * ours["delta"]
            error_over += not ours["within_bound"]
            # The bisection's answer is at most a tolerance above the truth, so this is how far
            # over the bound a run certainly is; at or below 1 means not at all.
            worst_error = max(worst_error, (ours["distance"] - ours["tolerance"]) / bound)

            kappa = min([r["n"]] + [k for d, k in known[r["source"]] if d <= ours["delta"]])
            tighter += kappa < r["n"]
            limit = 2 * max(kappa, 2) - 2
            size_over += ours["output_size"] > limit
            worst_size = max(worst_size, ours["output_size"] / limit)
        out[eps] = {"pairs": len(rows), "error_over": error_over, "worst_error": worst_error,
                    "size_over": size_over, "worst_size": worst_size, "kappa_bounded": tighter}
    return out


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
    return "\n".join(lines)


def pct(x):
    return f"{100 * x:.1f}%"


def matched_rows(stats):
    return [[SHORT[a], m, s["pairs"], f"{s['size_ratio_median']:.3f}", pct(s["smaller"]),
             pct(s["equal"]), pct(s["larger"]), f"{s['error_ratio_median']:.3f}",
             pct(s["ours_dominates"]), pct(s["baseline_dominates"])]
            for (a, m), s in stats]


MATCHED_HEADER = ["baseline", "m", "pairs", "size ours/theirs", "smaller", "equal", "larger",
                  "error ours/theirs", "ours dominates", "theirs dominates"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--matched", default=os.path.join("data", "error-matched"))
    parser.add_argument("--out", help="JSON output, default <matched>/comparison.json")
    parser.add_argument("--jobs", type=int, default=os.cpu_count(),
                        help="processes reading summaries, default one per core")
    args = parser.parse_args()

    pairs = load(args.matched, args.jobs)
    if not pairs:
        sys.exit(f"no summaries under {args.matched}/eps*/summary/ -- run ssk_match first")
    # Each baseline operating point appears once per eps; it is one operating point.
    base = list({(p["source"], p["algorithm"], p["m"]): p for p in pairs}.values())
    eps_values = sorted({p["eps"] for p in pairs})
    report = {"datasets": {}, "pooled": {}, "bounds": bounds(pairs, base), "notes": NOTES}

    for (dataset,), rows in sorted(group(base, "dataset").items()):
        ours = [p for p in pairs if p["dataset"] == dataset]
        rates = sorted({r["m"] for r in rows})
        cells = group(rows, "m", "algorithm")
        fixed = {(m, a): fixed_rate(c) for (m, a), c in cells.items()}
        by_eps = {eps: sorted(((a, m), matched(c)) for (a, m), c in
                              group([p for p in ours if p["eps"] == eps], "algorithm", "m").items())
                  for eps in eps_values}
        report["datasets"][dataset] = {
            "fixed_rate": {f"{SHORT[a]} m{m}": s for (m, a), s in fixed.items()},
            "matched": {f"eps={eps}": {f"{SHORT[a]} m{m}": s for (a, m), s in stats}
                        for eps, stats in by_eps.items()},
            "curves": {
                "baselines": {SHORT[a]: curve(c, lambda r: r["baseline"]["output_size"],
                                              lambda r: r["baseline"]["distance"])
                              for (a,), c in group(rows, "algorithm").items()},
                "ours": {f"eps={eps}": {SHORT[a]: curve(c, lambda r: r["ours"]["output_size"],
                                                        lambda r: r["ours"]["distance"])
                                        for (a,), c in group([p for p in ours if p["eps"] == eps],
                                                             "algorithm").items()}
                         for eps in eps_values},
            },
        }

        print(f"\n## {dataset}\n\n### Coverage: trajectories per cell\n")
        matches = group(ours, "m", "eps")
        print(table(["m"] + [SHORT[a] for a in ALGORITHMS] + [f"ours eps={e}" for e in eps_values],
                    [[m] + [len(cells.get((m, a), [])) for a in ALGORITHMS]
                     + [len(matches.get((m, e), [])) for e in eps_values] for m in rates]))

        print("\n### Fixed rate: median measured error (dataset units)\n")
        print(table(["m", "rate"] + [SHORT[a] for a in ALGORITHMS] + ["DOTS exact", "DOTS runs"],
                    [[m, 2 ** m]
                     + [f"{fixed[m, a]['error_median']:.4g}" if (m, a) in fixed else "-"
                        for a in ALGORITHMS]
                     + ([pct(fixed[m, "dots"]["exact_budget"]),
                         f"{fixed[m, 'dots']['runs_mean']:.1f}"]
                        if (m, "dots") in fixed else ["-", "-"])
                     for m in rates]))

        for eps, stats in by_eps.items():
            print(f"\n### Ours at each baseline's error, eps={eps}\n")
            print(table(MATCHED_HEADER, matched_rows(stats)))

    for eps in eps_values:
        stats = sorted(((a, m), matched(c)) for (a, m), c in
                       group([p for p in pairs if p["eps"] == eps], "algorithm", "m").items())
        report["pooled"][f"eps={eps}"] = {f"{SHORT[a]} m{m}": s for (a, m), s in stats}
        print(f"\n## All datasets, eps={eps} (ratios only)\n")
        print(table(MATCHED_HEADER, matched_rows(stats)))

    print("\n## Theorem 3, checked\n")
    print(table(["eps", "pairs", "over (1+eps)delta", "worst certain d/bound", "kappa bounded by a curve",
                 "over 2kappa-2", "worst size/bound"],
                [[eps, b["pairs"], b["error_over"], f"{b['worst_error']:.9f}", b["kappa_bounded"],
                  b["size_over"], f"{b['worst_size']:.3f}"]
                 for eps, b in report["bounds"].items()]))
    print("\nkappa is only bounded from above here, so passing the size check does not prove "
          "the bound; failing it does break it.")

    print("\n## Notes\n")
    for note in NOTES:
        print(f"- {note}")

    out = args.out or os.path.join(args.matched, "comparison.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"\nwritten: {out}")


if __name__ == "__main__":
    main()
