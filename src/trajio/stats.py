"""Counting: how many trajectories, how long they are."""

from __future__ import annotations

import math

from .core import TrajectorySource

__all__ = ["summarise", "LENGTH_BUCKETS"]

# Lengths span 1 to ~90 000 points, so linear bins would put everything in the first one.
LENGTH_BUCKETS = (2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)


def _percentile(values: list[int], q: float) -> float:
    if not values:
        return math.nan
    pos = q * (len(values) - 1)
    lo, hi = math.floor(pos), math.ceil(pos)
    if lo == hi:
        return float(values[lo])
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def summarise(source: TrajectorySource, *, limit: int | None = None) -> str:
    lengths: list[int] = []
    with_z = 0
    for track in source.tracks():
        if not track.points:
            continue
        lengths.append(len(track.points))
        if track.has_z():
            with_z += 1
        if limit is not None and len(lengths) >= limit:
            break
    if not lengths:
        return "no trajectories"

    lengths.sort()
    out = [
        f"trajectories : {len(lengths):,}",
        f"points       : {sum(lengths):,}",
        f"points/traj  : min {lengths[0]:,}  p50 {_percentile(lengths, .5):,.0f}  "
        f"p90 {_percentile(lengths, .9):,.0f}  p99 {_percentile(lengths, .99):,.0f}  "
        f"max {lengths[-1]:,}  mean {sum(lengths)/len(lengths):,.0f}",
    ]
    if 3 in source.dims:
        out.append(f"with altitude: {with_z:,} of {len(lengths):,} usable as 3D")

    counts = [0] * (len(LENGTH_BUCKETS) + 1)
    for n in lengths:
        i = 0
        while i < len(LENGTH_BUCKETS) and n >= LENGTH_BUCKETS[i]:
            i += 1
        counts[i] += 1
    peak = max(counts)
    out.append("")
    out.append("points/trajectory     count     share")
    for i, c in enumerate(counts):
        if not c:
            continue
        if i == 0:
            label = f"<{LENGTH_BUCKETS[0]}"
        elif i == len(LENGTH_BUCKETS):
            label = f">={LENGTH_BUCKETS[-1]}"
        else:
            label = f"{LENGTH_BUCKETS[i-1]}-{LENGTH_BUCKETS[i]-1}"
        out.append(f"{label:>12}  {c:>10,}  {100*c/len(lengths):5.1f}%  "
                   f"{'#' * max(1, round(44 * c / peak))}")
    return "\n".join(out)
