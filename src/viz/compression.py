"""Distribution of our compression rate at each baseline's error, per baseline and m.

    python src/viz/compression.py --out "data/error-matched/[2026-10-29] Results - Compression rate distribution"

Reads every <matched>/eps<E>/summary/<dataset>/<name>.json that ssk_match wrote. For each eps
it draws one PNG: a grid of histograms, a row per baseline and a column per m, of n/|ours|
over every trajectory, with the baseline's median n/|theirs| marked. Then it writes <out>.md
with the same distributions as tables, the gain over the baseline, and coverage per dataset.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import mean, median

import matplotlib

matplotlib.use("Agg")  # write files without needing a display

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402
import numpy as np  # noqa: E402

ALGORITHMS = ("douglas-peucker-n", "squish", "dots")
SHORT = {"douglas-peucker-n": "DPn", "squish": "SQUISH", "dots": "DOTS"}
RATES = range(1, 7)
# A bound under this many measurement tolerances (1e-9 x diagonal each) is below 1e-5 x the
# diagonal: the stored verdicts there are unreliable, see the 2026-10-28 incident note.
TINY = 1e4

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
BARS = "#2a78d6"


def read(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    return [{"dataset": doc["dataset"], "source": doc["source"], "eps": doc["eps"],
             "n": doc["input_size"],
             "algorithm": p["algorithm"], "m": p["m"],
             "ours": p["ours"]["output_size"], "theirs": p["baseline"]["output_size"],
             "error": p["baseline"]["distance"],
             "tiny": p["baseline"]["distance"] < TINY * p["baseline"]["tolerance"]}
            for p in doc["pairs"]]


def load(root: str) -> list[dict]:
    paths = sorted(glob.glob(os.path.join(root, "eps*", "summary", "*", "*.json")))
    with ProcessPoolExecutor() as pool:
        return [row for rows in pool.map(read, paths, chunksize=256) for row in rows]


def group(rows: list[dict], *keys: str) -> dict[tuple, list[dict]]:
    out = defaultdict(list)
    for row in rows:
        out[tuple(row[k] for k in keys)].append(row)
    return out


def ours_rate(r: dict) -> float:
    return r["n"] / r["ours"]


def theirs_rate(r: dict) -> float:
    return r["n"] / r["theirs"]


def quantiles(values: list[float], qs: tuple[float, ...]) -> list[float]:
    return list(np.quantile(values, qs))


def nice_ticks(lo: float, hi: float) -> list[float]:
    decades = range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1)
    for steps in ((1,), (1, 2, 5), (1, 1.5, 2, 3, 4, 6, 8)):
        ticks = [s * 10.0 ** k for k in decades for s in steps if lo <= s * 10.0 ** k <= hi]
        if len(ticks) >= 2:
            return ticks if len(ticks) <= 3 else ticks[::(len(ticks) + 2) // 3]
    return [lo, hi]


def histograms(rows: list[dict], eps: float, out: Path) -> None:
    # Sized for the report form's 5 in cell, so one file serves the note and the report.
    plt.rcParams.update({"font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6})
    cells = group(rows, "algorithm", "m")
    fig, axes = plt.subplots(len(ALGORITHMS), len(RATES), figsize=(5.05, 3.9), facecolor=SURFACE)
    for row, a in enumerate(ALGORITHMS):
        for col, m in enumerate(RATES):
            ax = axes[row][col]
            ax.set_facecolor(SURFACE)
            cell = cells[a, m]
            rates = [ours_rate(r) for r in cell]
            theirs = median(theirs_rate(r) for r in cell)
            lo, hi = quantiles(rates, (0.005, 0.995))
            lo, hi = min(lo, theirs / 1.05), max(hi, theirs * 1.05)
            weights = np.full(len(rates), 100.0 / len(rates))
            heights, _, _ = ax.hist(np.clip(rates, lo, hi), bins=np.geomspace(lo, hi, 25),
                                    weights=weights, color=BARS, edgecolor=SURFACE,
                                    linewidth=0.5, zorder=2)
            ax.set_ylim(0, 1.45 * heights.max())
            ax.axvline(theirs, ymax=0.74, color=INK, linestyle=(0, (3, 2)), linewidth=0.8, zorder=3)
            ax.text(0.97, 0.96, f"ours {median(rates):.3g}×\ntheirs {theirs:.3g}×",
                    transform=ax.transAxes, ha="right", va="top", fontsize=5,
                    color=INK_SECONDARY)
            ax.set_xscale("log")
            ax.set_xlim(lo, hi)
            ax.set_xticks(nice_ticks(lo, hi))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(
                lambda x, _: f"{x / 1000:g}k" if x >= 1000 else f"{x:g}"))
            ax.xaxis.set_minor_locator(ticker.NullLocator())
            ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
            ax.tick_params(colors=INK_SECONDARY, labelsize=5, length=0, pad=2)
            for spine in ("top", "right", "left"):
                ax.spines[spine].set_visible(False)
            ax.spines["bottom"].set_color(GRID)
            if row == 0:
                ax.set_title(f"m = {m}  (N/{2 ** m})", color=INK)
            if col == 0:
                ax.set_ylabel(f"vs {SHORT[a]}\n% of pairs", color=INK)
    fig.supxlabel("ours: n / output size, log scale; dashed line: the baseline's median "
                  "n / output size", fontsize=6, color=INK_SECONDARY)
    fig.suptitle(f"Our compression rate at each baseline's measured Fréchet error, "
                 f"eps = {eps}", fontsize=7, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=0.6, w_pad=0.4)
    fig.savefig(out, dpi=300, facecolor=SURFACE)
    plt.close(fig)


def table(header: list[str], rows: list[list]) -> str:
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join("---" if i < 2 else "---:" for i in range(len(header))) + "|"]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
    return "\n".join(lines)


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def distribution_rows(rows: list[dict]) -> list[list]:
    out = []
    for (a, m), cell in sorted(group(rows, "algorithm", "m").items(),
                               key=lambda kv: (ALGORITHMS.index(kv[0][0]), kv[0][1])):
        rates = [ours_rate(r) for r in cell]
        p10, p25, p50, p75, p90 = quantiles(rates, (0.1, 0.25, 0.5, 0.75, 0.9))
        gain = [r["theirs"] / r["ours"] for r in cell]
        out.append([SHORT[a], m, f"{len(cell):,}", f"{median(theirs_rate(r) for r in cell):.2f}",
                    f"{p10:.2f}", f"{p25:.2f}", f"**{p50:.2f}**", f"{p75:.2f}", f"{p90:.2f}",
                    f"{mean(rates):.2f}", f"**{median(gain):.3f}**",
                    pct(mean(r["ours"] < r["theirs"] for r in cell)),
                    pct(mean(r["ours"] == r["theirs"] for r in cell)),
                    pct(mean(r["ours"] > r["theirs"] for r in cell))])
    return out


DISTRIBUTION_HEADER = ["baseline", "m", "pairs", "theirs median", "ours p10", "p25", "median",
                       "p75", "p90", "mean", "gain median", "ours smaller", "equal",
                       "ours larger"]


def by_dataset_rows(rows: list[dict], datasets: list[str]) -> list[list]:
    cells = group(rows, "algorithm", "m", "dataset")
    out = []
    for a in ALGORITHMS:
        for m in RATES:
            line = [SHORT[a], m]
            for d in datasets:
                cell = cells.get((a, m, d))
                line.append(f"{median(ours_rate(r) for r in cell):.2f} "
                            f"({median(r['theirs'] / r['ours'] for r in cell):.2f})"
                            if cell else "-")
            out.append(line)
    return out


def report(rows: list[dict], out: Path, trajectories: dict[str, int]) -> str:
    eps_values = sorted({r["eps"] for r in rows})
    datasets = sorted({r["dataset"] for r in rows})
    lines = [f"# {out.name}", "",
             "How many vertices ours needs at each baseline's own measured Fréchet error.",
             "Ours runs at delta = d*/(1+eps), so its guarantee is exactly the baseline's error d*,"
             " and the compression rate is n / output size.", "",
             "- **theirs median**: the baseline's n / output size. It is set by the budget "
             "N/2^m, so it is close to 2^m.",
             "- **ours p10 ... p90, mean**: our n / output size over every trajectory in the "
             "cell.",
             "- **gain**: theirs / ours in vertices, the factor by which ours compresses more. "
             "Above 1 means fewer vertices than the baseline at the same error.",
             "- **ours smaller / equal / larger**: the share of pairs where ours has fewer, "
             "the same or more vertices than the baseline.", "",
             "## Coverage", ""]
    per = group(rows, "dataset", "eps")
    lines.append(table(["dataset", "trajectories on disk"]
                       + [f"eps={e}: trajectories / pairs" for e in eps_values],
                       [[d, f"{trajectories.get(d, 0):,}"]
                        + [f"{len({r['source'] for r in per[d, e]}):,} / {len(per[d, e]):,}"
                           if (d, e) in per else "-" for e in eps_values]
                        for d in datasets]))
    lines += ["", "A pair is one trajectory at one baseline and one m. Datasets not listed "
              "have no results on this machine yet.", ""]

    lines += ["## Headline: every m pooled", ""]
    headline = []
    for eps in eps_values:
        for a in ALGORITHMS:
            cell = [r for r in rows if r["eps"] == eps and r["algorithm"] == a]
            headline.append([eps, SHORT[a], f"{len(cell):,}",
                             f"**{median(r['theirs'] / r['ours'] for r in cell):.3f}**",
                             pct(mean(r["ours"] < r["theirs"] for r in cell)),
                             pct(mean(r["ours"] == r["theirs"] for r in cell)),
                             pct(mean(r["ours"] > r["theirs"] for r in cell))])
    lines += [table(["eps", "baseline", "pairs", "gain median", "ours smaller", "equal",
                     "ours larger"], headline), "",
              "Each histogram has its own axis, log scale, from the cell's 0.5th to 99.5th "
              "percentile; the few values beyond are counted in the end bars. Heights are "
              "the share of that cell's pairs.", ""]

    for eps in eps_values:
        at = [r for r in rows if r["eps"] == eps]
        png = f"{out.name} eps{eps}.png"
        lines += [f"## eps = {eps}", "", f"![eps = {eps}](<{png}>)", "",
                  "### All datasets pooled", "", table(DISTRIBUTION_HEADER, distribution_rows(at)),
                  "", "### Per dataset: ours' median n / size (median gain)", "",
                  table(["baseline", "m"] + datasets, by_dataset_rows(at, datasets)), ""]

    lines += ["## Why the gains differ: the baselines' own errors", "",
              "Median over trajectories of the baseline's measured error over DPn's at the same "
              "budget. It does not depend on eps. Ours is matched to that error, so a baseline "
              "that is far less accurate than DPn at a budget leaves ours room for far fewer "
              "vertices.", ""]
    first = [r for r in rows if r["eps"] == eps_values[0]]
    dpn = {(r["source"], r["m"]): r["error"] for r in first if r["algorithm"] == ALGORITHMS[0]}
    cells = group([r for r in first if r["algorithm"] != ALGORITHMS[0]
                   and dpn.get((r["source"], r["m"]), 0) > 0], "algorithm", "m", "dataset")
    lines += [table(["baseline", "m"] + datasets,
                    [[SHORT[a], m] + [f"{median(r['error'] / dpn[r['source'], m] for r in cells[a, m, d]):.2f}"
                                      if (a, m, d) in cells else "-" for d in datasets]
                     for a in ALGORITHMS[1:] for m in RATES]), ""]

    lines += ["## Bounds below 1e-5 × diagonal", "",
              "The Fréchet check is unreliable there (archive incident note, 2026-10-28). It does "
              "not change the sizes above, which are counted, not measured. They are shown pooled "
              "without those pairs as a check that they do not drive the medians.", ""]
    summary = []
    for eps in eps_values:
        for a in ALGORITHMS:
            cell = [r for r in rows if r["eps"] == eps and r["algorithm"] == a]
            rest = [r for r in cell if not r["tiny"]]
            summary.append([eps, SHORT[a], f"{len(cell):,}", pct(1 - len(rest) / len(cell)),
                            f"{median(r['theirs'] / r['ours'] for r in cell):.3f}",
                            f"{median(r['theirs'] / r['ours'] for r in rest):.3f}"])
    lines += [table(["eps", "baseline", "pairs, all m", "tiny bound", "gain median, all",
                     "gain median, without tiny"], summary), ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--matched", default=os.path.join("data", "error-matched"))
    parser.add_argument("--trajectories", default=os.path.join("data", "trajectories"))
    parser.add_argument("--out", required=True, help="path prefix for the .md and the PNGs")
    args = parser.parse_args()

    rows = load(args.matched)
    out = Path(args.out)
    for eps in sorted({r["eps"] for r in rows}):
        histograms([r for r in rows if r["eps"] == eps], eps, Path(f"{out} eps{eps}.png"))
    trajectories = {d.name: len(list(d.glob("*.json")))
                    for d in Path(args.trajectories).iterdir() if d.is_dir()}
    with open(f"{out}.md", "w", encoding="utf-8") as f:
        f.write(report(rows, out, trajectories))
    print(f"written: {out}.md and {len({r['eps'] for r in rows})} PNGs")


if __name__ == "__main__":
    main()
