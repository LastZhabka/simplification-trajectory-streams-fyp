"""Histograms of the benchmark's time per input vertex, from what scripts/benchmark.py wrote.

    python scripts/plot_benchmark.py [--times data/benchmark/times] [--out <times>/../plots]

Needs matplotlib, as scripts/plot_results.py does (the venv has it). Writes

  <out>/histograms/<ours>_eps<E>.png   one per implementation and eps
  <out>/01_overview.png                all of them, a row per implementation, a column per eps
  <out>/02_by_eps.png                  per implementation, its eps on one axis
  <out>/03_eps_lines_log.png           the median against eps, a line per implementation, log y
  <out>/04_eps_lines_linear.png        the same, linear y
  <out>/05_eps_band.png                as 03, shaded from p10 to p90
  <out>/06_writeup_eps.pdf             05 for the write-up: white, no title, labels on the lines
  <out>/07_writeup_distributions.pdf   02 for the write-up: one axis, a row per implementation

06 and 07 also go out as .png. They are drawn at the size they are printed, two in a row of a
6.5 in text block, so their 8 pt is the page's.

A cell is one trajectory at one delta, and its time per vertex is its best over N. A cell that
hit the time limit has only a lower bound, which is no value to bin, so its panel counts it in
a note instead.
"""

import argparse
import json
import math
import os
from statistics import median

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullLocator  # noqa: E402

from plot_results import INK, INK_2, MUTED, SURFACE, style  # noqa: E402

OURS = ("greedy-stabbing-shared", "greedy-stabbing-naive",
        "greedy-stabbing-naive-paper-boundary", "greedy-stabbing-naive-paper",
        "greedy-stabbing-no-decomp", "greedy-stabbing-no-decomp-paper-boundary",
        "greedy-stabbing-no-decomp-paper")
# The write-up's names: Algorithms 1 and 2, Algorithm 1 on the paper's grid, and each of the three
# with the paper's own step, which keeps S[p] as a polygon: Direct-Grid is the paper's algorithm.
SHORT = {"greedy-stabbing-shared": "Algorithm 2 (shared)",
         "greedy-stabbing-naive": "Algorithm 1 (naive)",
         "greedy-stabbing-naive-paper-boundary": "Paper-Boundary",
         "greedy-stabbing-naive-paper": "Paper-Grid",
         "greedy-stabbing-no-decomp": "Direct-m-gon",
         "greedy-stabbing-no-decomp-paper-boundary": "Direct-Boundary",
         "greedy-stabbing-no-decomp-paper": "Direct-Grid"}
# A direct stabber is drawn as its decomposed twin is, dashed.
TWIN = {"greedy-stabbing-no-decomp": "greedy-stabbing-naive",
        "greedy-stabbing-no-decomp-paper-boundary": "greedy-stabbing-naive-paper-boundary",
        "greedy-stabbing-no-decomp-paper": "greedy-stabbing-naive-paper"}
BAR = "#2a78d6"
SERIES = {"greedy-stabbing-shared": "#2a78d6", "greedy-stabbing-naive": "#eb6834",
          "greedy-stabbing-naive-paper-boundary": "#1baf7a",
          "greedy-stabbing-naive-paper": "#eda100"}
SERIES.update({o: SERIES[t] for o, t in TWIN.items()})
DASH = (0, (3, 1.6))
# Blue 250 to 650, light to dark as eps shrinks and the work grows; 250 is the lightest an
# ordinal step may be on this surface.
RAMP = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")

# The write-up's figures: on white, in the text's Computer Modern.
WRITEUP = {"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
           "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
           "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
           "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK_2,
           "xtick.color": INK_2, "ytick.color": INK_2, "xtick.labelcolor": INK,
           "ytick.labelcolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
           "savefig.facecolor": "white", "axes.grid": False, "pdf.fonttype": 42}
WRITEUP_SIZE = (3.15, 2.4)
# The same bottom and top in both, so their x axes line up side by side.
WRITEUP_FRAME = {"bottom": 0.17, "top": 0.98, "right": 0.97}
NAME = {"greedy-stabbing-shared": "Algorithm 2", "greedy-stabbing-naive": "Algorithm 1",
        "greedy-stabbing-naive-paper-boundary": "Paper-Boundary",
        "greedy-stabbing-naive-paper": "Paper-Grid",
        "greedy-stabbing-no-decomp": r"Direct-$m$-gon",
        "greedy-stabbing-no-decomp-paper-boundary": "Direct-Boundary",
        "greedy-stabbing-no-decomp-paper": "Direct-Grid"}
# Told apart in grey print too; a direct stabber has its twin's marker, hollow.
MARKER = {"greedy-stabbing-shared": "o", "greedy-stabbing-naive": "s",
          "greedy-stabbing-naive-paper-boundary": "^", "greedy-stabbing-naive-paper": "D"}
MARKER.update({o: MARKER[t] for o, t in TWIN.items()})
# The ridgeline's rows, bottom up, each direct stabber over its twin.
ROWS = ("greedy-stabbing-shared", "greedy-stabbing-naive", "greedy-stabbing-no-decomp",
        "greedy-stabbing-naive-paper-boundary", "greedy-stabbing-no-decomp-paper-boundary",
        "greedy-stabbing-naive-paper", "greedy-stabbing-no-decomp-paper")


def ns(x):
    if x >= 1e6:
        return f"{x / 1e6:.3g}M"
    if x >= 1e3:
        return f"{x / 1e3:.3g}k"
    return f"{x:.3g}"


def quantile(xs, q):
    return xs[min(len(xs) - 1, int(q * len(xs)))]


def read(times):
    data, settings, trajectories = {}, None, 0
    for o in OURS:
        with open(f"{times}/{o}.json", encoding="utf-8") as f:
            doc = json.load(f)
        settings, trajectories = doc["settings"], len(doc["trajectories"])
        for t in doc["trajectories"]:
            for c in t["cells"]:
                cell = data.setdefault((o, c["eps"]), {"values": [], "bounds": []})
                if "ns_per_vertex" in c:
                    cell["values"].append(c["ns_per_vertex"])
                elif "ns_per_vertex_at_least" in c:
                    cell["bounds"].append(c["ns_per_vertex_at_least"])
    for cell in data.values():
        cell["values"].sort()
    return data, settings, trajectories


def log_edges(lo, hi, bins):
    a, b = math.log10(lo), math.log10(hi)
    return [10 ** (a + (b - a) * i / bins) for i in range(bins + 1)]


def log_ticks(ax, lo, hi):
    marks = []
    for steps in ((1, 2, 5), (1, 2, 3, 5, 7), range(1, 10)):
        marks = [m * 10 ** k for k in range(int(math.floor(math.log10(lo))),
                                             int(math.ceil(math.log10(hi))) + 1)
                 for m in steps if lo <= m * 10 ** k <= hi]
        if len(marks) >= 3:
            break
    ax.xaxis.set_major_locator(FixedLocator(marks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([ns(m) for m in marks])


def note(cell, limit):
    v = cell["values"]
    text = (f"median {ns(median(v))} ns\np10-p90 {ns(quantile(v, 0.1))}-{ns(quantile(v, 0.9))}"
            f"\nn = {len(v):,} cells")
    if cell["bounds"]:
        text += (f"\n+{len(cell['bounds'])} over the {limit:g} s limit,\n"
                 f"at least {ns(min(cell['bounds']))}-{ns(max(cell['bounds']))}, not binned")
    return text


def draw(ax, cell, limit, bins, compact=False):
    v = cell["values"]
    lo, hi = v[0] / 1.05, v[-1] * 1.05
    ax.hist(v, bins=log_edges(lo, hi, bins), color=BAR, edgecolor=SURFACE, linewidth=0.6)
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    log_ticks(ax, lo, hi)
    ax.axvline(median(v), color=INK_2, linewidth=1, linestyle=(0, (3, 2)))
    ax.grid(axis="x", visible=False)
    ax.set_axisbelow(True)
    # On the side away from the median, over the low tail rather than the peak.
    right = math.log(median(v) / lo) < 0.5 * math.log(hi / lo)
    ax.text(0.98 if right else 0.02, 0.95, note(cell, limit), transform=ax.transAxes,
            ha="right" if right else "left", va="top", fontsize=7 if compact else 9,
            color=INK_2)


def histogram(data, o, eps, limit, footer, path):
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=130)
    draw(ax, data[(o, eps)], limit, 40)
    ax.set_title(f"{SHORT[o]} at ε = {eps:g}: time per input vertex", loc="left", fontsize=11,
                 color=INK)
    ax.set_xlabel("ns per input vertex (log scale); dashed: median")
    ax.set_ylabel("cells")
    fig.text(0.012, 0.005, footer, fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path)
    plt.close(fig)


def overview(data, eps_list, limit, footer, path):
    fig, axes = plt.subplots(len(OURS), len(eps_list), figsize=(3.0 * len(eps_list), 2.3 * len(OURS)),
                             dpi=130, squeeze=False)
    for i, o in enumerate(OURS):
        for j, eps in enumerate(eps_list):
            ax = axes[i][j]
            draw(ax, data[(o, eps)], limit, 30, compact=True)
            ax.tick_params(labelsize=7)
            if i == 0:
                ax.set_title(f"ε = {eps:g}", fontsize=10, color=INK)
            if j == 0:
                ax.set_ylabel(f"{SHORT[o]}\ncells", fontsize=9)
    fig.suptitle("Time per input vertex, ns (log scale), a row per implementation", x=0.01,
                 ha="left", fontsize=12, color=INK)
    fig.text(0.01, 0.005, footer, fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(path)
    plt.close(fig)


def by_eps(data, eps_list, limit, footer, path):
    fig, axes = plt.subplots(len(OURS), 1, figsize=(8.0, 2.4 * len(OURS)), dpi=130)
    for ax, o in zip(axes, OURS):
        cells = [data[(o, eps)] for eps in eps_list]
        lo = min(c["values"][0] for c in cells) / 1.05
        hi = max(c["values"][-1] for c in cells) * 1.05
        edges = log_edges(lo, hi, int(36 * math.log10(hi / lo)))
        for eps, cell, color in zip(eps_list, cells, RAMP):
            counts, _, _ = ax.hist(cell["values"], bins=edges, histtype="step", linewidth=2,
                                   color=color, label=f"ε = {eps:g}")
            peak = max(range(len(counts)), key=lambda k: counts[k])
            ax.text(math.sqrt(edges[peak] * edges[peak + 1]), counts[peak], f"{eps:g}",
                    ha="center", va="bottom", fontsize=8, color=INK_2)
        ax.set_xscale("log")
        ax.set_xlim(lo, hi)
        log_ticks(ax, lo, hi)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
        ax.grid(axis="x", visible=False)
        ax.set_ylabel("cells")
        bounds = sum(len(c["bounds"]) for c in cells)
        ax.set_title(SHORT[o] + (f"  (+{bounds} over the {limit:g} s limit, not binned)"
                                 if bounds else ""), loc="left", fontsize=10, color=INK)
    axes[0].legend(loc="upper right", fontsize=8, ncol=len(eps_list))
    axes[-1].set_xlabel("ns per input vertex (log scale); each curve labelled with its ε")
    fig.suptitle("Time per input vertex as ε shrinks", x=0.01, ha="left", fontsize=12,
                 color=INK)
    fig.text(0.01, 0.005, footer, fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    fig.savefig(path)
    plt.close(fig)


def slope(xs, ys):
    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    return (sum((a - mx) * (b - my) for a, b in zip(lx, ly)) /
            sum((a - mx) ** 2 for a in lx))


def against_eps(data, eps_list, limit, footer, path, log_y, band):
    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=130)
    xs = sorted(eps_list)
    for o in OURS:
        cells = [data[(o, eps)]["values"] for eps in xs]
        mid = [median(v) for v in cells]
        color = SERIES[o]
        if band:
            ax.fill_between(xs, [quantile(v, 0.1) for v in cells],
                            [quantile(v, 0.9) for v in cells], color=color, alpha=0.2,
                            linewidth=0)
        ax.plot(xs, mid, color=color, linewidth=2, marker="o", markersize=7,
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                label=f"{SHORT[o]} (k = {slope(xs, mid):.2f})")
        if log_y:
            ax.annotate(SHORT[o], (xs[0], mid[0]), xytext=(8, 0), textcoords="offset points",
                        va="center", fontsize=9, color=INK_2)
    # From coarse to fine: eps falls to the right, as the work grows.
    span = xs[-1] - xs[0]
    ax.set_xlim(xs[-1] + 0.04 * span, xs[0] - (0.2 if log_y else 0.04) * span)
    ax.xaxis.set_major_locator(FixedLocator(xs))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{e:g}" for e in xs])
    if log_y:
        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(NullLocator())
    else:
        ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: ns(y)))
    ax.set_xlabel("ε")
    ax.set_ylabel("ns per input vertex" + (" (log scale)" if log_y else ""))
    ax.legend(loc="upper left", fontsize=9, title=r"k fitted over the five ε: time $\propto$ ε$^k$",
              title_fontsize=8)
    if not log_y:
        top = {o: max(median(data[(o, eps)]["values"]) for eps in xs) for o in OURS}
        flat = [o for o in OURS if top[o] < 0.02 * ax.get_ylim()[1]]
        ax.text(0.02, 0.6, f"On this scale {', '.join(SHORT[o] for o in flat)}\n"
                           f"are near the axis: at most {ns(max(top[o] for o in flat))} ns.",
                transform=ax.transAxes, ha="left", va="top", fontsize=9, color=INK_2)
    what = "median and p10-p90 band" if band else "median"
    ax.set_title(f"Time per input vertex against ε, {what}", loc="left", fontsize=11,
                 color=INK)
    bounds = sum(len(data[(o, eps)]["bounds"]) for o in OURS for eps in xs)
    fig.text(0.012, 0.005, footer + (f" {bounds} over the {limit:g} s limit left out."
                                     if bounds else ""), fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path)
    plt.close(fig)


def save_writeup(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def writeup_eps(data, eps_list, path):
    fig, ax = plt.subplots(figsize=WRITEUP_SIZE)
    fig.subplots_adjust(left=0.16, **WRITEUP_FRAME)
    xs = sorted(eps_list)
    # A direct stabber dashed over its twin, hollow and without a band, as its spread is its
    # twin's; only the four decomposed ones are named, and a key tells dashed from solid.
    at = {}
    for o in OURS:
        cells = [data[(o, eps)]["values"] for eps in xs]
        mid = [median(v) for v in cells]
        at[o] = math.sqrt(mid[-1] * mid[-2])
        if o in TWIN:
            ax.plot(xs, mid, color=SERIES[o], linewidth=1.0, linestyle=DASH, marker=MARKER[o],
                    markersize=3.2, markerfacecolor="white", markeredgecolor=SERIES[o],
                    markeredgewidth=0.8)
            continue
        ax.fill_between(xs, [quantile(v, 0.1) for v in cells], [quantile(v, 0.9) for v in cells],
                        color=SERIES[o], alpha=0.2, linewidth=0)
        ax.plot(xs, mid, color=SERIES[o], linewidth=1.2, marker=MARKER[o], markersize=3.5,
                markeredgecolor="white", markeredgewidth=0.6)
    # Each name halfway between the two coarsest eps: under its line, in its band, when nothing
    # is within half a decade below, else over the line and its dashed twin, never between them.
    twin_of = {t: o for o, t in TWIN.items()}
    for o in (o for o in OURS if o not in TWIN):
        pair = [at[o]] + ([at[twin_of[o]]] if o in twin_of else [])
        under = [y for y in at.values() if y < min(pair)]
        below = not under or math.log10(min(pair) / max(under)) >= 0.5
        ax.annotate(NAME[o], ((xs[-1] + xs[-2]) / 2, min(pair) if below else max(pair)),
                    xytext=(0, -5 if below else 4), textcoords="offset points", ha="center",
                    va="top" if below else "bottom", fontsize=7)
    key = [Line2D([], [], color=INK_2, linewidth=1.2, marker="o", markersize=3.5,
                  markeredgecolor="white", markeredgewidth=0.6, label="decomposed"),
           Line2D([], [], color=INK_2, linewidth=1.0, linestyle=DASH, marker="o",
                  markersize=3.2, markerfacecolor="white", markeredgecolor=INK_2,
                  markeredgewidth=0.8, label="direct")]
    ax.legend(handles=key, loc="upper left", frameon=False, fontsize=7, handlelength=2.6)
    span = xs[-1] - xs[0]
    ax.set_xlim(xs[-1] + 0.04 * span, xs[0] - 0.04 * span)
    ax.xaxis.set_major_locator(FixedLocator(xs))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels([f"{e:g}" for e in xs])
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.grid(axis="y", which="major", color="#e6e6e6", linewidth=0.5)
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("time per input vertex (ns)")
    save_writeup(fig, path)


def density(values, grid):
    # Gaussian kernel density over log10 of the values, Silverman's bandwidth.
    v = np.log10(values)
    q75, q25 = np.percentile(v, [75, 25])
    h = 0.9 * min(v.std(), (q75 - q25) / 1.34) * len(v) ** -0.2
    z = (grid[:, None] - v[None, :]) / h
    return np.exp(-0.5 * z * z).sum(axis=1) / (len(v) * h * math.sqrt(2 * math.pi))


def shades(color, n):
    # Tints to shades of one hue, light for coarse eps, dark for fine.
    c = np.array(to_rgb(color))
    return [c + (1 - c) * -s if s < 0 else c * (1 - s) for s in np.linspace(-0.55, 0.5, n)]


def writeup_distributions(data, eps_list, path):
    fig, ax = plt.subplots(figsize=WRITEUP_SIZE)
    fig.subplots_adjust(left=0.03, **WRITEUP_FRAME)
    lo = min(data[(o, e)]["values"][0] for o in ROWS for e in eps_list) / 1.6
    hi = max(data[(o, e)]["values"][-1] for o in ROWS for e in eps_list) * 1.6
    grid = np.linspace(math.log10(lo), math.log10(hi), 1200)
    curves = {(o, e): density(data[(o, e)]["values"], grid) for o in ROWS for e in eps_list}
    x = 10 ** grid
    height = 0.82
    for row, o in enumerate(ROWS):
        ax.axhline(row, color="#d9d9d9", linewidth=0.5, zorder=0)
        # Each row to the same height, under the next row's baseline, so no row hides another.
        scale = height / max(curves[(o, e)].max() for e in eps_list)
        for eps, color in zip(eps_list, shades(SERIES[o], len(eps_list))):
            y = scale * curves[(o, eps)]
            # Without the flat tails, which would draw over the baseline the whole width.
            seen = np.nonzero(y > 0.004)[0]
            part = slice(seen[0], seen[-1] + 1)
            ax.fill_between(x[part], row, row + y[part], color=color,
                            alpha=0.08 if o in TWIN else 0.15, linewidth=0)
            ax.plot(x[part], row + y[part], color=color, linewidth=0.8 if o in TWIN else 0.9,
                    linestyle=DASH if o in TWIN else "solid")
        values = [data[(o, e)]["values"] for e in eps_list]
        # The lower rows have room on the right, the upper ones on the left.
        if row < len(ROWS) / 2:
            at, ha = max(quantile(v, 0.99) for v in values) * 1.4, "left"
        else:
            at, ha = min(quantile(v, 0.01) for v in values) / 1.4, "right"
        ax.text(at, row + 0.08, NAME[o], ha=ha, va="bottom", fontsize=7)
    top = ROWS[-1]
    scale = height / max(curves[(top, e)].max() for e in eps_list)
    for eps, text in ((eps_list[0], rf"$\varepsilon = {eps_list[0]:g}$"),
                      (eps_list[-1], f"{eps_list[-1]:g}")):
        peak = int(curves[(top, eps)].argmax())
        ax.text(x[peak], len(ROWS) - 1 + scale * curves[(top, eps)][peak] + 0.05, text,
                ha="center", va="bottom", fontsize=7)
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.set_ylim(-0.05, len(ROWS) - 1 + height + 0.3)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("time per input vertex (ns)")
    save_writeup(fig, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--times", default="data/benchmark/times")
    ap.add_argument("--out")
    args = ap.parse_args()
    out = args.out or os.path.join(os.path.dirname(os.path.normpath(args.times)), "plots")
    data, settings, trajectories = read(args.times)
    eps_list = sorted({e for _, e in data}, reverse=True)
    limit = settings["time_limit"]
    cells = len(data[(OURS[0], eps_list[0])]["values"])
    footer = (f"{trajectories:,} trajectories × {len(settings['rates'])} δ = {cells:,} cells "
              f"per implementation and ε. A cell: one trajectory at one δ, its best time "
              f"over N.")

    style()
    os.makedirs(f"{out}/histograms", exist_ok=True)
    for o in OURS:
        for eps in eps_list:
            histogram(data, o, eps, limit, footer, f"{out}/histograms/{o}_eps{eps:g}.png")
    overview(data, eps_list, limit, footer, f"{out}/01_overview.png")
    by_eps(data, eps_list, limit, footer, f"{out}/02_by_eps.png")
    against_eps(data, eps_list, limit, footer, f"{out}/03_eps_lines_log.png", True, False)
    against_eps(data, eps_list, limit, footer, f"{out}/04_eps_lines_linear.png", False, False)
    against_eps(data, eps_list, limit, footer, f"{out}/05_eps_band.png", True, True)
    with plt.rc_context(WRITEUP):
        writeup_eps(data, eps_list, f"{out}/06_writeup_eps")
        writeup_distributions(data, eps_list, f"{out}/07_writeup_distributions")
    print(f"wrote {len(OURS) * len(eps_list)} histograms and 01-07 under {out}")


if __name__ == "__main__":
    main()
