"""Generate test trajectories as JSON.

    python src/viz/generate.py spiral2d      -n 600 -o data/synthetic/spiral2d.json
    python src/viz/generate.py shrink-spiral -n 75  -o data/synthetic/shrink_spiral3d.json

`shrink-spiral` is the 3D example asked for in the meeting: a spiral whose radius decays
fast, with a floor on the spacing between consecutive points so that they never collapse
onto each other.
"""

from __future__ import annotations

import argparse
import math
import random

from trajectory_io import make_trajectory, save


def spiral2d(n: int, turns: float = 3.0, growth: float = 0.30) -> list[list[float]]:
    pts = []
    for k in range(n):
        t = 2.0 * math.pi * turns * k / n
        r = growth * t
        pts.append([r * math.cos(t), r * math.sin(t)])
    return pts


def zigzag2d(n: int, step: float = 0.5, amp: float = 1.2) -> list[list[float]]:
    return [[step * k, 0.0 if k % 2 == 0 else amp] for k in range(n)]


def random_walk2d(n: int, seed: int = 42, scale: float = 1.0) -> list[list[float]]:
    rng = random.Random(seed)
    x = y = 0.0
    pts = []
    for _ in range(n):
        x += rng.uniform(-scale, scale)
        y += rng.uniform(-scale, scale)
        pts.append([x, y])
    return pts


def u_shape2d(n: int = 200, gap: float = 0.6, height: float = 12.0) -> list[list[float]]:
    """Two near-parallel legs joined by a bend.

    When `gap` is below 2*delta the two legs' neighbourhood balls overlap, so a segment
    straight up the middle meets every ball. Only the *ordered* stabbing condition rejects
    it -- the return leg would have to map backwards along such a segment. This is the
    sharpest small test of the ordering requirement.
    """
    per_leg = max(2, (n - 8) // 2)
    pts = [[0.0, height * k / per_leg] for k in range(per_leg + 1)]
    pts += [[gap * k / 8.0, height] for k in range(1, 8)]
    pts += [[gap, height * k / per_leg] for k in range(per_leg, -1, -1)]
    return pts


def o_shape2d(n: int = 200, radius: float = 6.0) -> list[list[float]]:
    """A closed loop: the last point repeats the first.

    Both endpoints are pinned to the same place, so the simplification cannot collapse --
    the Frechet matching still has to travel the whole way round.
    """
    return [[radius * math.cos(2.0 * math.pi * k / n),
             radius * math.sin(2.0 * math.pi * k / n)] for k in range(n + 1)]


def shrink_spiral3d(
    n: int = 75,
    turns: float = 4.0,
    r0: float = 12.0,
    decay: float = 0.965,
    rise: float = 0.55,
    min_gap: float = 0.9,
) -> list[list[float]]:
    """A 3D spiral that shrinks quickly but keeps consecutive points apart.

    The radius decays geometrically (r_k = r0 * decay**k). Left alone that would make the
    later points crowd together, so the angular step is widened whenever the chord between
    consecutive points would fall below `min_gap`. The z-rise per step is constant, which
    also puts a hard floor of `rise` on the spacing.
    """
    pts: list[list[float]] = []
    theta = 0.0
    base_step = 2.0 * math.pi * turns / n
    for k in range(n):
        r = r0 * (decay**k)
        z = rise * k
        pts.append([r * math.cos(theta), r * math.sin(theta), z])

        # Choose the next angular step so the chord is at least min_gap.
        r_next = r0 * (decay ** (k + 1))
        step = base_step
        for _ in range(64):
            dx = r_next * math.cos(theta + step) - r * math.cos(theta)
            dy = r_next * math.sin(theta + step) - r * math.sin(theta)
            if math.sqrt(dx * dx + dy * dy + rise * rise) >= min_gap:
                break
            step *= 1.15
        theta += step
    return pts


GENERATORS = {
    "spiral2d": lambda a: spiral2d(a.n, a.turns, a.growth),
    "zigzag2d": lambda a: zigzag2d(a.n),
    "walk2d": lambda a: random_walk2d(a.n, a.seed),
    "u-shape": lambda a: u_shape2d(a.n, a.gap),
    "o-shape": lambda a: o_shape2d(a.n, a.radius),
    "shrink-spiral": lambda a: shrink_spiral3d(a.n),
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=sorted(GENERATORS))
    ap.add_argument("-n", type=int, default=400, help="number of points")
    ap.add_argument("-o", "--out", required=True, help="output JSON path")
    ap.add_argument("--turns", type=float, default=3.0)
    ap.add_argument("--growth", type=float, default=0.30)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--gap", type=float, default=0.6, help="U-shape leg separation")
    ap.add_argument("--radius", type=float, default=6.0, help="O-shape radius")
    args = ap.parse_args()

    pts = GENERATORS[args.kind](args)
    save(args.out, make_trajectory(pts, name=args.kind))
    print(f"wrote {len(pts)} points ({len(pts[0])}D) to {args.out}")


if __name__ == "__main__":
    main()
