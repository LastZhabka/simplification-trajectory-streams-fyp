"""Writing tracks out as JSON.

    <out>/<dataset>-<ord>.json    {"dim": d, "name": ..., "t_unit": ..., "t": [...],
                                   "points": [[x, y], ...]}
    <out>/index.json              the dataset, and one entry per file written

This is the format the rest of the project reads: `src/algo/io/trajectory.hpp` on the C++
side, `src/viz/` on the Python side.

Coordinates and timestamps are written as the **source text**, placed directly in the JSON
number slot, so no float round-trip happens between the dataset file and ours. The one
exception is a value whose text is not legal JSON number syntax (a leading `+` or `.`, a
leading zero); those are re-rendered through a double, which preserves the value but not the
spelling.

`t` is present only when the source has a clock, and then holds one value per point. What
those values mean is declared by `t_unit`: `unix_ms` for an absolute epoch, `ms` for
milliseconds from the start of the sequence.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from .core import Track, TrajectorySource

__all__ = ["export"]

# JSON's number grammar (RFC 8259): optional -, int with no leading zeros, optional
# fraction, optional exponent. Anything matching can be emitted verbatim.
_JSON_NUMBER = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?\Z")


def _number(text: str) -> str:
    """The source text if JSON can carry it as-is, else a value-preserving rendering."""
    if _JSON_NUMBER.match(text):
        return text
    return repr(float(text))  # repr of a double round-trips exactly


def _write_trajectory(
    path: Path, name: str, dim: int, points: list[tuple], times: list[str], unit: str
) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write('{\n  "dim": %d,\n  "name": %s,\n' % (dim, json.dumps(name)))
        if times:
            fh.write('  "t_unit": %s,\n' % json.dumps(unit))
            fh.write('  "t": [%s],\n' % ", ".join(_number(t) for t in times))
        fh.write('  "points": [\n')
        last = len(points) - 1
        for k, (x, y, z) in enumerate(points):
            coords = (_number(x), _number(y), _number(z)) if dim == 3 else (
                _number(x), _number(y))
            fh.write("    [%s]%s\n" % (", ".join(coords), "" if k == last else ","))
        fh.write("  ]\n}\n")


def export(
    source: TrajectorySource,
    out_dir: str | Path,
    *,
    dims: int = 2,
    limit: int | None = None,
    digits: int = 6,
) -> list[tuple[str, str, int, int]]:
    """Write one JSON file per track. Returns the index rows.

    `dims=3` is honoured per track: a track missing any `z` is written 2D instead.
    """
    if dims not in source.dims:
        raise ValueError(
            f"{source.name} supports dims {list(source.dims)}, not {dims} ({source.axes})"
        )

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    index: list[tuple[str, str, int, int]] = []

    for track in source.tracks():
        if not track.points:
            continue
        n = len(index) + 1
        use_z = dims == 3 and track.has_z()
        name = f"{source.name}-{n:0{digits}d}.json"

        times = track.times if track.has_time() else []
        _write_trajectory(
            out / name, track.id, 3 if use_z else 2, track.points, times, source.time_unit
        )

        index.append((name, track.id, len(track.points), 3 if use_z else 2))
        if limit is not None and len(index) >= limit:
            break

    summary: dict = {"dataset": source.name, "count": len(index)}
    if source.time_unit:
        summary["t_unit"] = source.time_unit
    # The frame rate is a guess, not a reading, so record which one produced these files.
    if getattr(source, "fps", None) is not None:
        summary["fps"] = source.fps
    summary["trajectories"] = [
        {"file": f, "name": sid, "points": n, "dim": d} for f, sid, n, d in index
    ]

    with open(out / "index.json", "w", encoding="utf-8", newline="\n") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    return index
