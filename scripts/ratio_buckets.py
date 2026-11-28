"""How ours compares with each baseline, pair by pair, as histograms of ratios -- step 10 of
docs/running.md.

    python scripts/ratio_buckets.py [--error-matched data/error-matched]
                                    [--error-matched-unscaled data/error-matched-unscaled]
                                    [--size-matched data/size-matched]
                                    [--out data/ratio-buckets.json] [--jobs N]

Three comparisons, each read from the per-trajectory summaries its step wrote, and each keyed by
m, the baseline at N / 2^m vertices:

  fixed-rate             step 7: ours at delta = d*/(1+eps), within the baseline's measured
                         error d*.
  fixed-rate-unscaled    step 8: as step 7 with ours at delta = d*, within (1+eps) d*.
  fixed-size             step 9: ours at the delta a bisection found it to fit in the
                         baseline's own size. Here the distance ratio is the result, and the
                         size ratio at most 100.

The comparison is taken from the summary's `delta_rule`, so a tree in the wrong place is
labelled by what it is. For every pair, two ratios, as percentages of the baseline's:

  size          ours' vertices / the baseline's vertices * 100. Under 100: ours is smaller.
  distance      ours' measured Frechet distance / the baseline's * 100: how much of the error
                the baseline spent ours spent. Where the baseline is within its tolerance,
                at 0 within measurement, and ours within the two tolerances, it is 100.

each counted into the buckets [0, 10), [10, 20), ..., [190, 200), [200, inf), separately for
every (algorithm, level, eps) and implementation of ours, pooled over the datasets and per
dataset. Each group also counts the pairs over ours' error bound (1+eps)*delta, over the
baseline's error (steps 8 and 9 allow it), and over Theorem 3's size bound 2*kappa - 2, kappa bounded
by the input and by every baseline in the same summary within delta. Prints the pooled tables
as Markdown; writes all of it to --out.
"""

import argparse
import glob
import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from statistics import median

EDGES = list(range(0, 201, 10))
LABELS = [f"{a}-{a + 10}" for a in EDGES[:-1]] + ["200+"]
WHAT = {"fixed-rate": "step 7: baseline at N/2^m, ours at delta = d*/(1+eps)",
        "fixed-rate-unscaled": "step 8: as step 7, ours at delta = d*",
        "fixed-size": "step 9: baseline at N/2^m, ours bisected on delta to its size"}
COMPARISON = {"d*/(1+eps)": "fixed-rate", "d*": "fixed-rate-unscaled", "size": "fixed-size"}
SHORT = {"douglas-peucker-n": "DPn", "squish": "SQUISH", "squish-e": "SQUISH-E", "dots": "DOTS"}
ALGORITHMS = tuple(SHORT)


def bucket(percent):
    return min(int(percent // 10), len(LABELS) - 1)


def ratio(ours, baseline):
    if baseline > 0:
        return 100.0 * ours / baseline
    return 100.0 if ours == 0 else float("inf")


def distance_ratio(ours, base):
    # A baseline within its tolerance is at 0 within measurement, and a ratio of two numbers under
    # the tolerance is noise: ours within the two tolerances is the same error.
    if base["distance"] <= base["tolerance"] and \
            ours["distance"] <= ours["tolerance"] + base["tolerance"]:
        return 100.0
    return ratio(ours["distance"], base["distance"])


def read(path):
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    comparison = COMPARISON[doc["delta_rule"]]
    rows = []
    for p in doc["pairs"]:
        ours, base = p["ours"], p["baseline"]
        # Every distance on disk is the top of a bisection, so that curve is within it.
        kappa = min([doc["input_size"]] + [q["baseline"]["output_size"] for q in doc["pairs"]
                                           if q["baseline"]["distance"] <= ours["delta"]])
        rows.append((comparison, doc["ours"], p["algorithm"], f"m{p['m']}", doc["eps"],
                     doc["dataset"],
                     ratio(ours["output_size"], base["output_size"]),
                     distance_ratio(ours, base),
                     not ours["within_bound"],
                     ours["distance"] - ours["tolerance"] > base["distance"],
                     ours["output_size"] > 2 * max(kappa, 2) - 2))
    return rows


def summarise(values, counts):
    buckets = [0] * len(LABELS)
    for v in values:
        buckets[bucket(v)] += 1
    n = len(values)
    finite = sorted(v for v in values if v != float("inf"))
    return {"n": n, "buckets": dict(zip(LABELS, buckets)),
            "percent": {b: 100.0 * c / n for b, c in zip(LABELS, buckets)},
            "median": median(values), "mean": sum(finite) / len(finite) if finite else None,
            "under_100": sum(v < 100 for v in values) / n,
            "equal_100": sum(v == 100 for v in values) / n,
            **counts}


def level_order(level):
    return float(level[1:])


def table(title, groups, metric, comparison):
    keys = sorted((k for k in groups if k[0] == comparison),
                  key=lambda k: (k[1], ALGORITHMS.index(k[2]), level_order(k[3]), -k[4]))
    if not keys:
        return []
    lines = [f"\n### {title}\n",
             "| ours | algorithm | m | eps | n | median | < 100% | "
             + " | ".join(LABELS) + " |",
             "|---" * (len(LABELS) + 7) + "|"]
    for k in keys:
        s = groups[k][metric]
        lines.append(f"| {k[1]} | {SHORT[k[2]]} | {k[3]} | {k[4]:g} | {s['n']} | "
                     f"{s['median']:.1f} | {100 * s['under_100']:.1f} | "
                     + " | ".join(f"{s['percent'][b]:.1f}" for b in LABELS) + " |")
    return lines


def checks(groups):
    lines = ["\n### Checks, per comparison and eps\n",
             "| comparison | eps | pairs | over (1+eps)delta | over the baseline's error | "
             "over 2kappa-2 |", "|---|---|---|---|---|---|"]
    totals = defaultdict(lambda: [0, 0, 0, 0])
    for k, g in groups.items():
        s = g["size"]
        t = totals[(k[0], k[4])]
        t[0] += s["n"]
        t[1] += s["over_bound"]
        t[2] += s["over_error"]
        t[3] += s["over_size"]
    for (comparison, eps), (n, bound, error, size) in sorted(totals.items(),
                                                            key=lambda kv: (kv[0][0], -kv[0][1])):
        lines.append(f"| {comparison} | {eps:g} | {n} | {bound} | {error} | {size} |")
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--error-matched", default="data/error-matched")
    ap.add_argument("--error-matched-unscaled", default="data/error-matched-unscaled")
    ap.add_argument("--size-matched", default="data/size-matched")
    ap.add_argument("--out", default="data/ratio-buckets.json")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()

    paths = []
    for root in (args.error_matched, args.error_matched_unscaled, args.size_matched):
        paths += glob.glob(os.path.join(root, "*", "eps*", "summary", "*", "*.json"))
    if not paths:
        raise SystemExit("no summaries found: run steps 7 to 9 first")

    values = defaultdict(lambda: {"size": [], "distance": [],
                                  "counts": {"over_bound": 0, "over_error": 0, "over_size": 0}})
    with ProcessPoolExecutor(args.jobs) as pool:
        for rows in pool.map(read, paths, chunksize=256):
            for (comparison, ours, algorithm, level, eps, dataset, size, dist, over_bound,
                 over_error, over_size) in rows:
                for key in ((comparison, ours, algorithm, level, eps, "all"),
                            (comparison, ours, algorithm, level, eps, dataset)):
                    v = values[key]
                    v["size"].append(size)
                    v["distance"].append(dist)
                    v["counts"]["over_bound"] += over_bound
                    v["counts"]["over_error"] += over_error
                    v["counts"]["over_size"] += over_size

    groups = {k: {"size": summarise(v["size"], v["counts"]),
                  "distance": summarise(v["distance"], v["counts"])} for k, v in values.items()}
    pooled = {k[:5]: g for k, g in groups.items() if k[5] == "all"}

    lines = ["Ratios of ours to the baseline, in percent, pair by pair; each bucket column is "
             "the share of pairs in it. `< 100%` is the share where ours is smaller (size) or "
             "closer (distance). Pooled over the datasets."]
    lines += checks(pooled)
    for comparison in WHAT:
        lines += table(f"Size, {comparison} ({WHAT[comparison]})", pooled, "size", comparison)
        lines += table(f"Distance, {comparison} ({WHAT[comparison]})", pooled, "distance",
                       comparison)
    print("\n".join(lines))

    doc = {"buckets": LABELS, "groups": [
        {"comparison": k[0], "ours": k[1], "algorithm": k[2], "level": k[3], "eps": k[4],
         "dataset": k[5], **g} for k, g in sorted(groups.items(), key=lambda kv: str(kv[0]))]}
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out + ".part", "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=1)
    os.replace(args.out + ".part", args.out)
    print(f"\nwrote {args.out}: {len(groups)} groups")


if __name__ == "__main__":
    main()
