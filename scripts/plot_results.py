"""Figures for step 10 of docs/running.md, from what scripts/ratio_buckets.py wrote.

    python scripts/plot_results.py [--buckets data/ratio-buckets.json] [--out data/plots]
                                   [--eps 0.1] [--jobs N]

The one script that needs a package: matplotlib (the venv has it). Writes

  <out>/histograms/<comparison>/<metric>/<ours>/<algorithm>_<level>_eps<E>.png
        one histogram per row of ratio_buckets.py's tables, pooled over the datasets
  <out>/summary/
        01_win_tie_loss.png        share of pairs where ours is smaller, the same, larger
        02_median_size_ratio.png   median size ratio by m, one line per baseline
        03_eps_effect.png          the same against eps, at m = 4
        04_histograms_overview.png the size histograms at m = 4
        05_median_distance_ratio.png  as 02, for the measured Frechet distance
        06_delta_rule.png          step 7 against step 8: ours at delta = d*/(1+eps) beside
                                   delta = d*, median size ratio by m
        07_win_tie_loss_unscaled.png  as 01, for step 8
        08_median_distance_ratio_unscaled.png  as 05, for step 8
        09_fixed_size_distance.png as 05, for step 9: ours at each baseline's size

  <out>/writeup/
        01_fixed-rate.pdf ... 03_fixed-size.pdf, and .png
                                   the write-up's: per step, a row per baseline of the size
                                   histogram, or for step 9 the distance histogram, the rates
                                   pooled, in plot_benchmark.py's style
        04_all-experiments.pdf, and .png
                                   the three in the text's width: a row per step, a column per
                                   baseline, each column on one scale

01-05 are step 7, where ours is held to the baseline's error; 06-08 appear once step 8 has
run, 09 once step 9 has. The summary and write-up figures are at --eps unless they vary it.
"""

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = {"douglas-peucker-n": "#2a78d6", "squish": "#eb6834", "squish-e": "#1baf7a",
          "dots": "#eda100"}
SHORT = {"douglas-peucker-n": "DPn", "squish": "SQUISH", "squish-e": "SQUISH-E", "dots": "DOTS"}
SMALLER, SAME, LARGER = "#2a78d6", "#c3c2b7", "#e34948"
COMPARISONS = {
    "fixed-rate": ("Fixed compression rate (step 7)", "baseline at N/2^m, ours at delta = d*/(1+eps)"),
    "fixed-rate-unscaled": ("Fixed rate, delta = d* (step 8)", "as step 7, ours at delta = d*"),
    "fixed-size": ("Fixed size (step 9)", "baseline at N/2^m, ours bisected on delta to its size"),
}
MAIN = ("fixed-rate",)
UNSCALED = ("fixed-rate-unscaled",)
RULE = "#898781"
REPRESENTATIVE = {"fixed-rate": "m4"}
METRICS = {"size": "ours' vertices / baseline's", "distance": "ours' Frechet distance / baseline's"}

# The write-up's figures: on white, in the text's Computer Modern, drawn 3.15 in wide to sit two
# in a row of a 6.5 in text block. plot_benchmark.py draws its own in the same style.
WRITEUP = {"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
           "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
           "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
           "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK_2,
           "xtick.color": INK_2, "ytick.color": INK_2, "xtick.labelcolor": INK,
           "ytick.labelcolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
           "savefig.facecolor": "white", "axes.grid": False, "pdf.fonttype": 42}
WRITEUP_SIZE = (3.15, 2.4)
# Top down from the baseline nearest ours to the furthest.
WRITEUP_ORDER = ("douglas-peucker-n", "dots", "squish-e", "squish")


def style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
        "font.size": 10, "text.color": INK, "axes.labelcolor": INK_2,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    })


def level_value(level):
    return float(level[1:])


def bucket_colors(labels):
    return [SMALLER if int(b.split("-")[0].rstrip("+")) < 100 else
            SAME if b.startswith("100-") else LARGER for b in labels]


def draw_histogram(ax, group, labels, metric, compact=False):
    s = group[metric]
    values = [s["percent"][b] for b in labels]
    ax.bar(range(len(labels)), values, width=0.62, color=bucket_colors(labels), linewidth=0)
    ax.axvline(9.5, color=INK_2, linewidth=1)
    ax.set_xticks(range(0, len(labels), 2 if not compact else 4))
    ax.set_xticklabels([labels[i].split("-")[0] if not labels[i].endswith("+") else "200+"
                        for i in range(0, len(labels), 2 if not compact else 4)])
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.grid(axis="x", visible=False)
    ax.set_axisbelow(True)
    ax.text(0.98, 0.95, f"median {s['median']:.0f}%\nsmaller {100 * s['under_100']:.0f}%\n"
                        f"n = {s['n']:,}", transform=ax.transAxes, ha="right", va="top",
            fontsize=8 if compact else 9, color=INK_2)


def histogram(args):
    group, labels, metric, path = args
    style()
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=130)
    draw_histogram(ax, group, labels, metric)
    title, sub = COMPARISONS[group["comparison"]]
    ax.set_title(f"{SHORT[group['algorithm']]} at {group['level']}, eps {group['eps']:g}: "
                 f"{metric} ratio", loc="left", fontsize=11, color=INK)
    fig.text(0.012, 0.005, f"{title}; {sub}. Bars left of the line: ours smaller "
                           f"({'fewer vertices' if metric == 'size' else 'closer'}).",
             fontsize=8, color=MUTED)
    ax.set_xlabel(f"{METRICS[metric]}, % (bucket start)")
    ax.set_ylabel("share of pairs, %")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path)
    plt.close(fig)
    return path


def pooled(groups):
    return {(g["comparison"], g["algorithm"], g["level"], g["eps"]): g
            for g in groups if g["dataset"] == "all" and g["ours"] == "greedy-stabbing-shared"}


def levels_of(table, comparison):
    return sorted({k[2] for k in table if k[0] == comparison}, key=level_value)


def win_tie_loss(table, eps, path, comparisons=MAIN, title=None):
    fig, axes = plt.subplots(len(comparisons), 4, figsize=(13, 0.6 + 2.5 * len(comparisons)),
                             dpi=130, sharex=True, squeeze=False)
    for row, comparison in enumerate(comparisons):
        levels = levels_of(table, comparison)
        for col, algorithm in enumerate(SERIES):
            ax = axes[row][col]
            for y, level in enumerate(levels):
                s = table[(comparison, algorithm, level, eps)]["size"]
                smaller, same = 100 * s["under_100"], 100 * s["equal_100"]
                larger = 100 - smaller - same
                left = 0
                for width, color in ((smaller, SMALLER), (same, SAME), (larger, LARGER)):
                    ax.barh(y, width - 0.6 if width > 1 else width, left=left, height=0.55,
                            color=color, linewidth=0)
                    left += width
                if smaller >= 12:
                    ax.text(1.5, y, f"{smaller:.0f}%", va="center", fontsize=8, color="#ffffff")
            ax.set_yticks(range(len(levels)))
            ax.set_yticklabels(levels)
            ax.invert_yaxis()
            ax.set_xlim(0, 100)
            ax.grid(axis="y", visible=False)
            if row == 0:
                ax.set_title(SHORT[algorithm], loc="left", fontsize=11)
            if col == 0:
                ax.set_ylabel(COMPARISONS[comparison][0], fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("share of pairs, %")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (SMALLER, SAME, LARGER)]
    fig.legend(handles, ["ours smaller", "same size", "ours larger"], loc="upper right",
               ncol=3, fontsize=9)
    fig.suptitle(title or f"How often ours needs fewer vertices than each baseline, eps {eps:g}",
                 x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def lines_by_level(table, eps, metric, path, title, comparisons=MAIN):
    fig, axes = plt.subplots(1, len(comparisons), figsize=(max(6 * len(comparisons), 9), 4.6),
                             dpi=130, squeeze=False)
    axes = axes[0]
    for ax, comparison in zip(axes, comparisons):
        levels = levels_of(table, comparison)
        series = {a: [table[(comparison, a, lv, eps)][metric]["median"] for lv in levels]
                  for a in SERIES}
        top = max(115.0, 1.08 * max(max(ys) for ys in series.values()))
        ends = sorted(ys[-1] for ys in series.values())
        # Converging end labels would stack on each other; the legend names the lines instead.
        label_ends = min(b - a for a, b in zip(ends, ends[1:])) > 0.04 * top
        for algorithm, ys in series.items():
            ax.plot(range(len(levels)), ys, color=SERIES[algorithm], linewidth=2, marker="o",
                    markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.5,
                    label=SHORT[algorithm], solid_capstyle="round")
            if label_ends:
                ax.annotate(SHORT[algorithm], (len(levels) - 1, ys[-1]), xytext=(8, 0),
                            textcoords="offset points", va="center", fontsize=9, color=INK_2)
        ax.axhline(100, color=INK_2, linewidth=1)
        above_clear = all(not 100 < ys[0] < 100 + 0.12 * top for ys in series.values())
        ax.text(-0.25, 100 + (0.035 if above_clear else -0.035) * top,
                f"same {metric} as the baseline", fontsize=8, color=MUTED,
                va="bottom" if above_clear else "top")
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels)
        ax.set_xlim(-0.3, len(levels) - 0.3 + 0.7)
        ax.set_ylim(0, top)
        ax.set_title(COMPARISONS[comparison][0], loc="left", fontsize=11)
        ax.set_xlabel("compression rate 2^m")
    axes[0].set_ylabel(f"median {METRICS[metric]}, %")
    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="upper right", ncol=4, fontsize=9)
    fig.suptitle(title.format(eps=eps), x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def eps_effect(table, path):
    fig, axes = plt.subplots(1, len(MAIN), figsize=(max(6 * len(MAIN), 9), 4.4), dpi=130,
                             sharey=True, squeeze=False)
    axes = axes[0]
    for ax, comparison in zip(axes, MAIN):
        level = REPRESENTATIVE[comparison]
        eps_list = sorted({k[3] for k in table if k[0] == comparison}, reverse=True)
        for algorithm, color in SERIES.items():
            ys = [table[(comparison, algorithm, level, e)]["size"]["median"] for e in eps_list]
            ax.plot(range(len(eps_list)), ys, color=color, linewidth=2, marker="o",
                    markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.5,
                    label=SHORT[algorithm])
            ax.annotate(SHORT[algorithm], (len(eps_list) - 1, ys[-1]), xytext=(8, 0),
                        textcoords="offset points", va="center", fontsize=9, color=INK_2)
        ax.axhline(100, color=INK_2, linewidth=1)
        ax.set_xticks(range(len(eps_list)))
        ax.set_xticklabels([f"{e:g}" for e in eps_list])
        ax.set_xlim(-0.3, len(eps_list) - 0.3 + 0.7)
        ax.set_ylim(0, 115)
        ax.set_title(f"{COMPARISONS[comparison][0]}, at {level}", loc="left", fontsize=11)
        ax.set_xlabel("eps (finer to the right)")
    axes[0].set_ylabel("median ours' vertices / baseline's, %")
    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="upper right", ncol=4, fontsize=9)
    fig.suptitle("A smaller eps helps a little: the median size ratio against eps",
                 x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def overview(table, labels, eps, path):
    fig, axes = plt.subplots(len(MAIN), 4, figsize=(14, 0.6 + 2.8 * len(MAIN)), dpi=130,
                             sharey=True, squeeze=False)
    for row, comparison in enumerate(MAIN):
        level = REPRESENTATIVE[comparison]
        for col, algorithm in enumerate(SERIES):
            ax = axes[row][col]
            draw_histogram(ax, table[(comparison, algorithm, level, eps)], labels, "size",
                           compact=True)
            if row == 0:
                ax.set_title(SHORT[algorithm], loc="left", fontsize=11)
            if col == 0:
                ax.set_ylabel(f"{COMPARISONS[comparison][0]}\nat {level}: share, %", fontsize=9)
            if row == len(MAIN) - 1:
                ax.set_xlabel("ours / baseline vertices, %")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (SMALLER, SAME, LARGER)]
    fig.legend(handles, ["ours smaller", "100-110%", "ours larger"], loc="upper right", ncol=3,
               fontsize=9)
    fig.suptitle(f"Size ratio, pair by pair, eps {eps:g}", x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def delta_rule(table, eps, path):
    fig, axes = plt.subplots(len(MAIN), 4, figsize=(14, 0.6 + 2.8 * len(MAIN)), dpi=130,
                             sharey=True, squeeze=False)
    for row, (scaled, unscaled) in enumerate(zip(MAIN, UNSCALED)):
        levels = levels_of(table, scaled)
        for col, algorithm in enumerate(SERIES):
            ax = axes[row][col]
            for comparison, color, name in ((scaled, RULE, "delta = error / (1+eps)"),
                                            (unscaled, SERIES[algorithm], "delta = error")):
                ys = [table[(comparison, algorithm, lv, eps)]["size"]["median"] for lv in levels]
                ax.plot(range(len(levels)), ys, color=color, linewidth=2, marker="o",
                        markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.5, label=name)
            ax.axhline(100, color=INK_2, linewidth=1)
            ax.set_xticks(range(len(levels)))
            ax.set_xticklabels(levels, fontsize=8)
            ax.set_ylim(0, 115)
            if row == 0:
                ax.set_title(SHORT[algorithm], loc="left", fontsize=11)
            if col == 0:
                ax.set_ylabel(f"{COMPARISONS[scaled][0].split(' (')[0]}\nmedian size ratio, %",
                              fontsize=9)
    handles = [plt.Line2D([], [], color=RULE, linewidth=2, marker="o"),
               plt.Line2D([], [], color=INK_2, linewidth=2, marker="o")]
    fig.legend(handles, ["delta = error / (1+eps): ours within the baseline's error",
                         "delta = error: ours within (1+eps) x it (line in the baseline's colour)"],
               loc="upper right", ncol=1, fontsize=9)
    fig.suptitle(f"What running ours at delta = d* instead buys in size, eps {eps:g}",
                 x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path)
    plt.close(fig)


def writeup_bars(ax, table, labels, eps, comparison, algorithm, metric):
    # The share of the baseline's pairs in each bucket of the ratio, the five levels pooled. The
    # line is the baseline's own; the grey bucket [100, 110) holds the ties.
    groups = [g for k, g in table.items()
              if k[0] == comparison and k[1] == algorithm and k[3] == eps]
    counts = [sum(g[metric]["buckets"][b] for g in groups) for b in labels]
    shares = [100.0 * c / sum(counts) for c in counts]
    ax.bar(range(len(labels)), shares, width=0.8, color=bucket_colors(labels), linewidth=0)
    ax.axvline(9.5, color=INK_2, linewidth=0.6)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylim(0, None)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2)
    return max(shares)


def writeup_histograms(table, labels, eps, comparison, metric, path):
    fig, axes = plt.subplots(len(WRITEUP_ORDER), 1, figsize=WRITEUP_SIZE, sharex=True)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.98, bottom=0.17, hspace=0.2)
    for ax, algorithm in zip(axes, WRITEUP_ORDER):
        writeup_bars(ax, table, labels, eps, comparison, algorithm, metric)
        ax.locator_params(axis="y", nbins=2)
        ax.text(0.99, 0.92, SHORT[algorithm], transform=ax.transAxes, ha="right", va="top",
                fontsize=7)
    axes[-1].set_xticks([-0.5, 4.5, 9.5, 14.5, 19.5])
    axes[-1].set_xticklabels(["0", "50", "100", "150", "200"])
    axes[-1].set_xlabel("ours' vertices / baseline's (%)" if metric == "size" else
                        "ours' distance / baseline's (%)")
    fig.text(0.015, 0.575, "pairs (%)", rotation=90, va="center", fontsize=8)
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


GRID_TITLES = {"fixed-rate": "1: fixed rate", "fixed-rate-unscaled": "2: rate, $\\delta = d^*$",
               "fixed-size": "3: fixed size"}


def writeup_grid(table, labels, eps, path):
    # All three experiments in the text's width: a row each, a column per baseline, each column
    # on one scale, the panels about square. The last row is a distance ratio, so a gap sets it
    # apart from the two sizes.
    comparisons = list(GRID_TITLES)
    fig = plt.figure(figsize=(6.5, 5.0))
    grid = fig.add_gridspec(len(comparisons) + 1, len(WRITEUP_ORDER),
                            height_ratios=[1, 1, 0.3, 1], left=0.1, right=0.99,
                            top=0.94, bottom=0.11, hspace=0.38, wspace=0.24)
    axes = [[fig.add_subplot(grid[row if row < 2 else row + 1, col])
             for col in range(len(WRITEUP_ORDER))] for row in range(len(comparisons))]
    tops = [[] for _ in WRITEUP_ORDER]
    for row, comparison in enumerate(comparisons):
        for col, algorithm in enumerate(WRITEUP_ORDER):
            ax = axes[row][col]
            tops[col].append(writeup_bars(ax, table, labels, eps, comparison, algorithm,
                                          "distance" if comparison == "fixed-size" else "size"))
            ax.set_xticks([-0.5, 4.5, 9.5, 14.5, 19.5])
            ax.set_xticklabels(["0", "50", "100", "150", "200"])
            ax.tick_params(labelsize=7.5)
        axes[row][0].set_ylabel(GRID_TITLES[comparison], fontsize=9)
    for col, algorithm in enumerate(WRITEUP_ORDER):
        axes[0][col].set_title(SHORT[algorithm], fontsize=9)
        for row in axes:
            row[col].set_ylim(0, 1.05 * max(tops[col]))
            row[col].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    for row, text in ((1, "ours' vertices / other's (%)"), (2, "ours' distance / other's (%)")):
        left, right = axes[row][0].get_position(), axes[row][-1].get_position()
        fig.text((left.x0 + right.x1) / 2, left.y0 - 0.068, text, ha="center", va="top",
                 fontsize=9)
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--buckets", default="data/ratio-buckets.json")
    ap.add_argument("--out", default="data/plots")
    ap.add_argument("--eps", type=float, default=0.1)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()

    with open(args.buckets, encoding="utf-8") as f:
        doc = json.load(f)
    labels = doc["buckets"]
    style()

    jobs = []
    for g in doc["groups"]:
        if g["dataset"] != "all":
            continue
        for metric in METRICS:
            d = os.path.join(args.out, "histograms", g["comparison"], metric, g["ours"])
            os.makedirs(d, exist_ok=True)
            jobs.append((g, labels, metric,
                         os.path.join(d, f"{g['algorithm']}_{g['level']}_eps{g['eps']:g}.png")))
    with ProcessPoolExecutor(args.jobs) as pool:
        done = list(pool.map(histogram, jobs, chunksize=8))
    print(f"wrote {len(done)} histograms under {os.path.join(args.out, 'histograms')}")

    table = pooled(doc["groups"])
    summary = os.path.join(args.out, "summary")
    os.makedirs(summary, exist_ok=True)
    win_tie_loss(table, args.eps, os.path.join(summary, "01_win_tie_loss.png"))
    lines_by_level(table, args.eps, "size", os.path.join(summary, "02_median_size_ratio.png"),
                   "Median size of ours against each baseline, eps {eps:g}")
    eps_effect(table, os.path.join(summary, "03_eps_effect.png"))
    overview(table, labels, args.eps, os.path.join(summary, "04_histograms_overview.png"))
    lines_by_level(table, args.eps, "distance",
                   os.path.join(summary, "05_median_distance_ratio.png"),
                   "Median Frechet distance of ours against each baseline, eps {eps:g}")
    if any(k[0] in UNSCALED for k in table):
        delta_rule(table, args.eps, os.path.join(summary, "06_delta_rule.png"))
        win_tie_loss(table, args.eps, os.path.join(summary, "07_win_tie_loss_unscaled.png"),
                     UNSCALED, "With delta = d*: how often ours needs fewer vertices, "
                               f"eps {args.eps:g} (ours allowed (1+eps) x the error)")
        lines_by_level(table, args.eps, "distance",
                       os.path.join(summary, "08_median_distance_ratio_unscaled.png"),
                       "With delta = d*: median Frechet distance of ours against each "
                       "baseline, eps {eps:g}", UNSCALED)
    if any(k[0] == "fixed-size" for k in table):
        lines_by_level(table, args.eps, "distance",
                       os.path.join(summary, "09_fixed_size_distance.png"),
                       "At the baseline's size: median distance of ours, eps {eps:g}",
                       ("fixed-size",))
    print(f"wrote the summary figures under {summary}")

    writeup = os.path.join(args.out, "writeup")
    os.makedirs(writeup, exist_ok=True)
    with plt.rc_context(WRITEUP):
        for k, comparison in enumerate(COMPARISONS, 1):
            if any(key[0] == comparison for key in table):
                writeup_histograms(table, labels, args.eps, comparison,
                                   "distance" if comparison == "fixed-size" else "size",
                                   os.path.join(writeup, f"{k:02d}_{comparison}"))
        if all(any(key[0] == c for key in table) for c in GRID_TITLES):
            writeup_grid(table, labels, args.eps, os.path.join(writeup, "04_all-experiments"))
    print(f"wrote the write-up figures under {writeup}")


if __name__ == "__main__":
    main()
