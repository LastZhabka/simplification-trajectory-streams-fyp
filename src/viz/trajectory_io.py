"""Shared JSON trajectory I/O for the Python tooling.

The same documents `src/trajio` exports and `src/algo/io/trajectory.hpp` reads.

Trajectory format (input)::

    {"dim": 2, "name": "spiral", "points": [[x, y], ...]}

Result format (output of the `simplify` CLI)::

    {"dim": 2, "algorithm": "...", "mode": "...",
     "points": [...],            # the simplified curve, sigma
     "input_points": [...],      # the original curve, tau
     "params": {...}, "stats": {...}, "frechet": {...}}
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def save(path: str | Path, doc: dict[str, Any]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(doc, handle, indent=2)
        handle.write("\n")


def make_trajectory(points: list[list[float]], name: str = "") -> dict[str, Any]:
    if not points:
        raise ValueError("trajectory must have at least one point")
    dim = len(points[0])
    if any(len(p) != dim for p in points):
        raise ValueError("all points must have the same dimension")
    doc: dict[str, Any] = {"dim": dim, "points": [[float(c) for c in p] for p in points]}
    if name:
        doc["name"] = name
    return doc


def points_of(doc: dict[str, Any], key: str = "points") -> list[list[float]]:
    return [list(map(float, p)) for p in doc[key]]
