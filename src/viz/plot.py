"""Render a trajectory or a simplification result to an image, with the statistics panel.

    python src/viz/plot.py data/synthetic/spiral2d.json -o data/renders/spiral2d.png

Handles both 2D and 3D results. The 3D case is drawn as a 3D line plot plus three
orthographic projections, since a single 3D view hides a lot of what the simplification
actually did.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # write files without needing a display

import matplotlib.pyplot as plt  # noqa: E402

from trajectory_io import load, points_of  # noqa: E402

TAU_STYLE = dict(color="#7a8ba3", lw=1.2, alpha=0.9, zorder=1)
SIGMA_STYLE = dict(color="#d1495b", lw=1.8, zorder=3)
SIGMA_PTS = dict(color="#d1495b", s=18, zorder=4, edgecolors="white", linewidths=0.5)


def stats_caption(doc: dict) -> str:
    st = doc.get("stats", {})
    pr = doc.get("params", {})
    lines = [
        f"algorithm : {doc.get('algorithm', '?')}  ({doc.get('mode', '?')})",
        f"eps / delta: {pr.get('eps', '?')} / {pr.get('delta', '?')}",
        f"input      : {st.get('input_size', '?')} points",
        f"output     : {st.get('output_size', '?')} points",
        f"compression: {st.get('compression', 0.0):.2f}x",
        f"avg latency: {st.get('avg_latency_us_per_point', 0.0):.3f} us/point",
        f"peak space : {st.get('peak_working_bytes', 0) / 1024.0:.1f} KB",
    ]
    fr = doc.get("frechet")
    if fr:
        ok = "OK" if fr.get("within_bound") else "VIOLATED"
        lines.append(
            f"frechet    : {fr.get('distance', 0.0):.4f} <= {fr.get('bound', 0.0):.4f}  [{ok}]"
        )
        lines.append(f"  computed by: {fr.get('computer', '?')}")
    return "\n".join(lines)


def plot_2d(doc: dict, out: str, title: str) -> None:
    tau = points_of(doc, "input_points") if "input_points" in doc else []
    sigma = points_of(doc, "points")

    fig = plt.figure(figsize=(12, 6.5))
    ax = fig.add_axes([0.05, 0.08, 0.62, 0.84])
    if tau:
        ax.plot([p[0] for p in tau], [p[1] for p in tau], label=f"input ({len(tau)})",
                **TAU_STYLE)
    ax.plot([p[0] for p in sigma], [p[1] for p in sigma], label=f"simplified ({len(sigma)})",
            **SIGMA_STYLE)
    ax.scatter([p[0] for p in sigma], [p[1] for p in sigma], **SIGMA_PTS)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="best", framealpha=0.9)
    ax.set_title(title)

    fig.text(0.70, 0.90, stats_caption(doc), va="top", ha="left", family="monospace",
             fontsize=10)
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_3d(doc: dict, out: str, title: str) -> None:
    tau = points_of(doc, "input_points") if "input_points" in doc else []
    sigma = points_of(doc, "points")

    fig = plt.figure(figsize=(13, 7.5))
    ax = fig.add_subplot(2, 3, 1, projection="3d")
    if tau:
        ax.plot([p[0] for p in tau], [p[1] for p in tau], [p[2] for p in tau], **TAU_STYLE)
    ax.plot([p[0] for p in sigma], [p[1] for p in sigma], [p[2] for p in sigma],
            **SIGMA_STYLE)
    ax.set_title("3D view", fontsize=10)

    # Orthographic projections: a single 3D view hides too much.
    for idx, (i, j, label) in enumerate([(0, 1, "xy"), (0, 2, "xz"), (1, 2, "yz")], start=2):
        axp = fig.add_subplot(2, 3, idx)
        if tau:
            axp.plot([p[i] for p in tau], [p[j] for p in tau], **TAU_STYLE)
        axp.plot([p[i] for p in sigma], [p[j] for p in sigma], **SIGMA_STYLE)
        axp.scatter([p[i] for p in sigma], [p[j] for p in sigma], **SIGMA_PTS)
        axp.set_title(f"{label} projection", fontsize=10)
        axp.set_aspect("equal", adjustable="datalim")
        axp.grid(alpha=0.25, lw=0.5)

    axt = fig.add_subplot(2, 3, 5)
    axt.axis("off")
    axt.text(0.0, 1.0, stats_caption(doc), va="top", ha="left", family="monospace",
             fontsize=9)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result", help="a trajectory, or a result document from a simplify run")
    ap.add_argument("-o", "--out", required=True, help="output image path")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    doc = load(args.result)
    # matplotlib will not create the directory, and data/renders/ is git-ignored, so a
    # fresh checkout would otherwise fail on the first plot.
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    dim = int(doc.get("dim", len(doc["points"][0])))
    title = args.title or f"{doc.get('algorithm', 'simplification')} - {args.result}"

    if dim == 2:
        plot_2d(doc, args.out, title)
    elif dim == 3:
        plot_3d(doc, args.out, title)
    else:
        raise SystemExit(f"cannot plot {dim}D trajectories")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
