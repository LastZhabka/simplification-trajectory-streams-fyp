"""The benchmark of docs/benchmark.md: every implementation of ours, timed on a fixed selection
of trajectories at a fixed selection of deltas, for each eps.

    python scripts/benchmark.py [--eps 0.5,0.25,0.2,0.1,0.05] [--time-limit 300] [--jobs 8]
                                [--per-dataset 10] [--out data/benchmark/benchmark.json]

Each cell -- one trajectory, delta, eps and implementation -- is one ssk_bench process, killed
at --time-limit seconds and then recorded as time_limit_exceeded. Writes, beside --out:

  benchmark.json         every cell, the selection and the settings
  times/<ours>.json      every time of one implementation, grouped by trajectory
  summary/<ours>.json    its distributions per eps, overall and per dataset: seconds,
                         ns per input vertex, output size
  pairwise.json          for each pair a, b: the distribution of b's time over a's

then prints the tables of docs/benchmark.md.
"""

import argparse
import json
import os
import platform
import random
import shutil
import subprocess
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from statistics import median

OURS = ("greedy-stabbing-shared", "greedy-stabbing-naive",
        "greedy-stabbing-naive-paper-boundary", "greedy-stabbing-naive-paper")
SHORT = {"greedy-stabbing-shared": "shared", "greedy-stabbing-naive": "naive",
         "greedy-stabbing-naive-paper-boundary": "paper-boundary",
         "greedy-stabbing-naive-paper": "paper"}
SHARED = OURS[0]
# Below this the shared stabber takes microseconds, which measures the timer.
MIN_POINTS = 500
# The long tail, beside the sample: the largest trajectory, and the one at the 99th percentile.
TAIL = (("geolife", "largest"), ("geolife", "p99"), ("mopsi", "largest"), ("mopsi", "p99"),
        ("ngsim-us-101", "largest"), ("ngsim-lankershim", "largest"),
        ("mot-mot20", "largest"), ("mot-dancetrack", "largest"))
# delta is where Douglas-Peucker-n, given N / 2^m vertices, measured its Frechet distance.
RATES = (2, 4, 6)
BASELINE = "douglas-peucker-n"
EXE = "build/src/pipelines/benchmark/ssk_bench"


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def deltas(distances, dataset, file):
    return {f"m{m}": read_json(f"{distances}/{BASELINE}/{dataset}/m{m}/{file}")
            ["frechet"]["distance"] for m in RATES}


def entry(trajectories, distances, dataset, t, rule):
    return {"dataset": dataset, "rule": rule, "file": t["file"], "name": t["name"],
            "points": t["points"], "path": f"{trajectories}/{dataset}/{t['file']}",
            "deltas": deltas(distances, dataset, t["file"])}


def select(trajectories, distances, per_dataset):
    chosen = {}
    for dataset in sorted(os.listdir(trajectories)):
        if not os.path.isfile(f"{trajectories}/{dataset}/index.json"):
            continue
        index = sorted(read_json(f"{trajectories}/{dataset}/index.json")["trajectories"],
                       key=lambda t: (t["points"], t["file"]))
        for d, rule in TAIL:
            if d == dataset:
                t = index[-1] if rule == "largest" else index[int(0.99 * (len(index) - 1))]
                chosen[(dataset, t["file"])] = entry(trajectories, distances, dataset, t, rule)
        # Seeded by the dataset's name, so adding a dataset leaves the others' samples alone.
        eligible = sorted((t for t in index if t["points"] >= MIN_POINTS),
                          key=lambda t: t["file"])
        random.Random(dataset).shuffle(eligible)
        taken = 0
        for t in eligible:
            if taken == per_dataset:
                break
            if (dataset, t["file"]) in chosen:
                continue
            e = entry(trajectories, distances, dataset, t, "sample")
            if min(e["deltas"].values()) > 0.0:
                chosen[(dataset, t["file"])] = e
                taken += 1
    return list(chosen.values())


def run(exe, limit, cell):
    args = [exe, "--in", cell["path"], "--ours", cell["ours"], "--eps", repr(cell["eps"]),
            "--delta", repr(cell["delta"])]
    started = time.monotonic()
    try:
        done = subprocess.run(args, capture_output=True, text=True, timeout=limit)
    except subprocess.TimeoutExpired:
        return {**cell, "status": "time_limit_exceeded", "wall": time.monotonic() - started}
    if done.returncode != 0:
        return {**cell, "status": "error", "stderr": done.stderr.strip()}
    return {**cell, "status": "ok", **json.loads(done.stdout)}


def groups(cells):
    by = defaultdict(dict)
    for c in cells:
        by[(c["dataset"], c["file"], c["rate"], c["eps"])][c["ours"]] = c
    return by


def speedups(gs, other, limit):
    ratios, bounded = [], 0
    for g in gs:
        s, o = g[SHARED], g[other]
        if s["status"] != "ok":
            continue
        if o["status"] == "ok":
            ratios.append(o["best"] / s["best"])
        elif o["status"] == "time_limit_exceeded":
            ratios.append(limit / s["best"])
            bounded += 1
    return ratios, bounded


def speedup_text(ratios, bounded):
    if not ratios:
        return "-"
    sign = ">= " if bounded * 2 >= len(ratios) else ""
    return f"{sign}{median(ratios):.1f}x (min {min(ratios):.1f}x)" + \
        (f", {bounded} bounds" if bounded else "")


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|---" * len(header) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def report(cells, limit, eps_list):
    by = groups(cells)
    per_eps = {eps: [g for k, g in by.items() if k[3] == eps] for eps in eps_list}
    lines = []

    lines.append(f"## Time summed per eps, over the cells that finished, of "
                 f"{len(per_eps[eps_list[0]])}; and how many hit the {limit:g} s limit\n")
    rows = []
    for eps in eps_list:
        row = [f"{eps:g}"]
        for o in OURS:
            mine = [g[o] for g in per_eps[eps]]
            ok = [c["best"] for c in mine if c["status"] == "ok"]
            tle = sum(c["status"] == "time_limit_exceeded" for c in mine)
            err = sum(c["status"] == "error" for c in mine)
            row.append(f"{sum(ok):.1f} s over {len(ok)}" + (f", {tle} TLE" if tle else "") +
                       (f", {err} error" if err else ""))
        rows.append(row)
    lines += table(["eps"] + [SHORT[o] for o in OURS], rows)

    lines.append("\n## Speedup of shared: median over cells, and its smallest\n")
    lines.append("Where the other implementation hit the limit, the speedup is at least "
                 "limit / shared and is counted at that lower bound; a median marked >= is one "
                 "where those are at least half the cells.\n")
    lines += table(["eps"] + [f"vs {SHORT[o]}" for o in OURS[1:]],
                   [[f"{eps:g}"] + [speedup_text(*speedups(per_eps[eps], o, limit))
                                    for o in OURS[1:]] for eps in eps_list])

    lines.append("\n## Pairwise: how many times faster the row is than the column, the median "
                 "over cells of column time / row time; cells where one hit the limit count at "
                 "the limit\n")
    for eps in eps_list:
        lines.append(f"eps {eps:g}\n")
        rows = []
        for a in OURS:
            row = [SHORT[a]]
            for b in OURS:
                if a == b:
                    row.append("-")
                    continue
                s = ratio_summary([(g[a], g[b]) for g in per_eps[eps]], limit)["with_bounds"]
                row.append(f"{s['median']:.2f}x" if s else "-")
            rows.append(row)
        lines += table(["faster than"] + [SHORT[o] for o in OURS], rows)
        lines.append("")

    lines.append("\n## Time per input vertex, ns: summed time over summed N of the cells "
                 "that finished\n")
    rows = []
    for eps in eps_list:
        row = [f"{eps:g}"]
        for o in OURS:
            ok = [g[o] for g in per_eps[eps] if g[o]["status"] == "ok"]
            ns = 1e9 * sum(c["best"] for c in ok) / max(1, sum(c["input_size"] for c in ok))
            row.append(f"{ns:,.0f}" + ("" if len(ok) == len(per_eps[eps]) else
                                       f" ({len(ok)} cells)"))
        rows.append(row)
    lines += table(["eps"] + [SHORT[o] for o in OURS], rows)

    lines.append("\n## Vertices: shared's extra over paper-boundary, on the cells where both "
                 "finished\n")
    rows = []
    for eps in eps_list:
        both = [g for g in per_eps[eps]
                if g[SHARED]["status"] == "ok" and g[OURS[2]]["status"] == "ok"]
        vs = sum(g[SHARED]["output_size"] for g in both)
        vb = sum(g[OURS[2]]["output_size"] for g in both)
        rows.append([f"{eps:g}", f"{len(both)}", f"{vs:,}", f"{vb:,}",
                     f"+{100 * (vs / vb - 1):.1f}%" if vb else "-"])
    lines += table(["eps", "cells", "shared", "paper-boundary", "shared's extra"], rows)

    datasets = sorted({k[0] for k in by})
    for o in OURS[1:]:
        lines.append(f"\n## Per dataset: median speedup of shared over {SHORT[o]}\n")
        rows = []
        for ds in datasets:
            row = [ds, f"{len({k[1] for k in by if k[0] == ds})}"]
            for eps in eps_list:
                gs = [g for k, g in by.items() if k[0] == ds and k[3] == eps]
                row.append(speedup_text(*speedups(gs, o, limit)))
            rows.append(row)
        lines += table(["dataset", "trajectories"] + [f"eps {e:g}" for e in eps_list], rows)

    lines.append("\n## The long tail: every cell of the largest and p99 trajectories, best-of "
                 "seconds, and vertices of shared / paper-boundary / paper\n")
    rows = []
    for k in sorted((k for k, g in by.items() if g[SHARED]["rule"] != "sample"),
                    key=lambda k: (k[0], k[1], k[2], -k[3])):
        g = by[k]
        times = [f"{g[o]['best']:.3f}" if g[o]["status"] == "ok" else
                 ("TLE" if g[o]["status"] == "time_limit_exceeded" else "error") for o in OURS]
        sizes = " / ".join(str(g[o].get("output_size", "-")) for o in (SHARED, *OURS[2:]))
        rows.append([f"{k[0]}/{k[1]}", f"{g[SHARED]['points']}", k[2], f"{k[3]:g}", *times,
                     sizes])
    lines += table(["trajectory", "N", "delta", "eps", *[SHORT[o] for o in OURS], "vertices"],
                   rows)

    agree = [g for g in by.values() if all(g[o]["status"] == "ok" for o in OURS[:2])]
    same = sum(g[OURS[0]]["digest"] == g[OURS[1]]["digest"] for g in agree)
    lines.append(f"\nshared and naive output identical to the bit in {same} of {len(agree)} "
                 "cells where both finished.")
    return "\n".join(lines)


def write_json(path, doc):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".part", "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=1)
    os.replace(path + ".part", path)


def dist(xs):
    if not xs:
        return None
    xs = sorted(xs)
    at = lambda q: xs[min(len(xs) - 1, int(q * len(xs)))]
    return {"count": len(xs), "min": xs[0], "p10": at(0.10), "p25": at(0.25),
            "median": median(xs), "p75": at(0.75), "p90": at(0.90), "p99": at(0.99),
            "max": xs[-1], "mean": sum(xs) / len(xs)}


def statuses(cs):
    return {s: sum(c["status"] == s for c in cs) for s in ("ok", "time_limit_exceeded", "error")}


def algorithm_summary(cs):
    ok = [c for c in cs if c["status"] == "ok"]
    return {"cells": len(cs), **statuses(cs),
            "seconds": dist([c["best"] for c in ok]),
            "ns_per_vertex": dist([1e9 * c["best"] / c["input_size"] for c in ok]),
            "output_size": dist([c["output_size"] for c in ok])}


def ratio_summary(pairs, limit):
    """`pairs` are (a, b) cells of one (trajectory, delta, eps). The ratio is b's time over a's:
    how many times faster a is. Where one hit the limit, the limit stands in for its time, so
    the ratio is a bound; `with_bounds` counts those, `finished` does not."""
    exact, bounded = [], []
    a_only = b_only = neither = 0
    for a, b in pairs:
        if a["status"] == "ok" and b["status"] == "ok":
            exact.append(b["best"] / a["best"])
        elif a["status"] == "ok" and b["status"] == "time_limit_exceeded":
            bounded.append(limit / a["best"])
            a_only += 1
        elif a["status"] == "time_limit_exceeded" and b["status"] == "ok":
            bounded.append(b["best"] / limit)
            b_only += 1
        else:
            neither += 1
    return {"finished": dist(exact), "with_bounds": dist(exact + bounded),
            "only_a_finished": a_only, "only_b_finished": b_only, "neither_finished": neither}


def per_vertex(cell, limit):
    if cell["status"] == "ok":
        return {"ns_per_vertex": 1e9 * cell["best"] / cell["points"]}
    if cell["status"] == "time_limit_exceeded":
        return {"ns_per_vertex_at_least": 1e9 * limit / cell["points"]}
    return {}


def write_derived(root, results, selection, settings):
    by = groups(results)
    eps_list, limit = settings["eps"], settings["time_limit"]
    datasets = sorted({t["dataset"] for t in selection})
    written = []

    for o in OURS:
        trajectories = []
        for t in selection:
            cells = sorted((g[o] for k, g in by.items()
                            if k[0] == t["dataset"] and k[1] == t["file"]),
                           key=lambda c: (c["rate"], -c["eps"]))
            trajectories.append({
                "dataset": t["dataset"], "file": t["file"], "name": t["name"],
                "rule": t["rule"], "points": t["points"],
                "cells": [{"rate": c["rate"], "delta": c["delta"], "eps": c["eps"],
                           "status": c["status"], "best": c.get("best"),
                           "seconds": c.get("seconds"), "output_size": c.get("output_size"),
                           **per_vertex(c, limit)}
                          for c in cells]})
        path = f"{root}/times/{o}.json"
        write_json(path, {"algorithm": o, "settings": settings, "trajectories": trajectories})
        written.append(path)

        mine = lambda eps, ds=None: [g[o] for k, g in by.items()
                                     if k[3] == eps and (ds is None or k[0] == ds)]
        path = f"{root}/summary/{o}.json"
        write_json(path, {
            "algorithm": o, "settings": settings,
            "per_eps": {f"{eps:g}": algorithm_summary(mine(eps)) for eps in eps_list},
            "per_dataset": {ds: {f"{eps:g}": algorithm_summary(mine(eps, ds))
                                 for eps in eps_list} for ds in datasets}})
        written.append(path)

    comparisons = []
    for i, a in enumerate(OURS):
        for b in OURS[i + 1:]:
            pairs = lambda eps, ds=None: [(g[a], g[b]) for k, g in by.items()
                                          if k[3] == eps and (ds is None or k[0] == ds)]
            comparisons.append({
                "a": a, "b": b, "ratio": f"time of {b} / time of {a}",
                "per_eps": {f"{eps:g}": ratio_summary(pairs(eps), limit) for eps in eps_list},
                "per_dataset": {ds: {f"{eps:g}": ratio_summary(pairs(eps, ds), limit)
                                     for eps in eps_list} for ds in datasets}})
    path = f"{root}/pairwise.json"
    write_json(path, {"settings": settings, "comparisons": comparisons})
    written.append(path)
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--eps", default="0.5,0.25,0.2,0.1,0.05")
    ap.add_argument("--time-limit", type=float, default=300.0)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--per-dataset", type=int, default=10)
    ap.add_argument("--trajectories", default="data/trajectories")
    ap.add_argument("--distances", default="data/frechet-distances")
    ap.add_argument("--out", default="data/benchmark/benchmark.json")
    args = ap.parse_args()

    exe = shutil.which(EXE)
    if exe is None:
        raise SystemExit(f"no {EXE}: build first (docs/running.md step 2)")
    eps_list = [float(e) for e in args.eps.split(",")]
    selection = select(args.trajectories, args.distances, args.per_dataset)
    print(f"{len(selection)} trajectories, {sum(t['points'] for t in selection):,} points",
          flush=True)

    cells = [{"dataset": t["dataset"], "file": t["file"], "rule": t["rule"],
              "points": t["points"], "path": t["path"], "rate": rate, "delta": delta,
              "eps": eps, "ours": ours}
             for t in selection for rate, delta in t["deltas"].items()
             for eps in eps_list for ours in OURS]
    # Shuffled, so no implementation always runs beside the slow ones or late in a warm machine.
    random.Random(0).shuffle(cells)

    started = time.monotonic()
    results = []
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run, exe, args.time_limit, c) for c in cells]
        for i, f in enumerate(as_completed(futures), 1):
            results.append(f.result())
            if i % 1000 == 0 or i == len(cells):
                print(f"  {i}/{len(cells)} cells, {time.monotonic() - started:.0f} s",
                      flush=True)

    settings = {"eps": eps_list, "time_limit": args.time_limit, "jobs": args.jobs,
                "per_dataset": args.per_dataset, "min_points": MIN_POINTS,
                "rates": list(RATES), "baseline": BASELINE,
                "machine": {"processor": platform.processor(), "cpus": os.cpu_count(),
                            "system": platform.platform()}}
    write_json(args.out, {"settings": settings, "selection": selection,
                          "wall_seconds": time.monotonic() - started, "cells": results})
    print(f"wrote {args.out}")
    for path in write_derived(os.path.dirname(args.out), results, selection, settings):
        print(f"wrote {path}")
    print()
    print(report(results, args.time_limit, eps_list))


if __name__ == "__main__":
    main()
