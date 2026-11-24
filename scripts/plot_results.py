"""Figures for step 11 of docs/running.md, from what scripts/ratio_buckets.py wrote.

    python scripts/plot_results.py [--buckets data/ratio-buckets.json] [--out data/plots]
                                   [--eps 0.1] [--jobs N]

The one script that needs a package: matplotlib (the venv has it). Writes

  <out>/histograms/<comparison>/<metric>/<ours>/<algorithm>_<level>_eps<E>.png
        one histogram per row of ratio_buckets.py's tables, pooled over the datasets
  <out>/summary/
        01_win_tie_loss.png        share of pairs where ours is smaller, the same, larger
        02_median_size_ratio.png   median size ratio by level, one line per baseline
        03_eps_effect.png          the same against eps, at one level per comparison
        04_histograms_overview.png the size histograms of one level per comparison
        05_median_distance_ratio.png  as 02, for the measured Frechet distance
        06_delta_rule.png          steps 7-8 against steps 9-10: ours at delta = error/(1+eps)
                                   beside delta = error, median size ratio by level
        07_win_tie_loss_unscaled.png  as 01, for steps 9 and 10
        08_median_distance_ratio_unscaled.png  as 05, for steps 9 and 10

01-05 are steps 7 and 8, where ours is held to the baseline's error; 06-08 appear once steps 9
and 10 have run. The summary figures are at --eps unless they vary it.
"""

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

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
    "fixed-error": ("Fixed error (step 8)", "baseline within T = f x diagonal, ours at delta = T/(1+eps)"),
    "fixed-rate-unscaled": ("Fixed rate, delta = d* (step 9)", "as step 7, ours at delta = d*"),
    "fixed-error-unscaled": ("Fixed error, delta = T (step 10)", "as step 8, ours at delta = T"),
}
MAIN = ("fixed-rate", "fixed-error")
UNSCALED = ("fixed-rate-unscaled", "fixed-error-unscaled")
RULE = "#898781"
REPRESENTATIVE = {"fixed-rate": "m4", "fixed-error": "f0.003"}
METRICS = {"size": "ours' vertices / baseline's", "distance": "ours' Frechet distance / baseline's"}


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
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.6), dpi=130, sharex=True)
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
    for ax in axes[1]:
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
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), dpi=130)
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
        ax.set_xlabel("compression rate 2^m" if comparison.startswith("fixed-rate")
                      else "error target, fraction of the diagonal")
    axes[0].set_ylabel(f"median {METRICS[metric]}, %")
    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="upper right", ncol=4, fontsize=9)
    fig.suptitle(title.format(eps=eps), x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def eps_effect(table, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), dpi=130, sharey=True)
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
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.2), dpi=130, sharey=True)
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
            if row == 1:
                ax.set_xlabel("ours / baseline vertices, %")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (SMALLER, SAME, LARGER)]
    fig.legend(handles, ["ours smaller", "100-110%", "ours larger"], loc="upper right", ncol=3,
               fontsize=9)
    fig.suptitle(f"Size ratio, pair by pair, eps {eps:g}", x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path)
    plt.close(fig)


def delta_rule(table, eps, path):
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.2), dpi=130, sharey=True)
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
    fig.suptitle(f"What running ours at delta = d* or T instead buys in size, eps {eps:g}",
                 x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path)
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
                     UNSCALED, "With delta = d* or T: how often ours needs fewer vertices, "
                               f"eps {args.eps:g} (ours allowed (1+eps) x the error)")
        lines_by_level(table, args.eps, "distance",
                       os.path.join(summary, "08_median_distance_ratio_unscaled.png"),
                       "With delta = d* or T: median Frechet distance of ours against each "
                       "baseline, eps {eps:g}", UNSCALED)
    print(f"wrote the summary figures under {summary}")


if __name__ == "__main__":
    main()
