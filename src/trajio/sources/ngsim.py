"""NGSIM Vehicle Trajectories (US DOT), the `8ect-6jqj` export.

One flat CSV of 11.85 M rows at 10 Hz. There is no file-per-trajectory, so rows are grouped
by **`(Location, Vehicle_ID, Section_ID)`** -- `Vehicle_ID` alone is not a trajectory id: on
the arterial sites it is reused by distinct vehicles at the same instant, which welds two
cars hundreds of metres apart into one zigzagging track.

x = Global_X, y = Global_Y (State Plane NAD83 feet, as stored). No z.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from ..core import Option, Track, TrajectorySource, register_source
from ..grouping import group_by
from ._csvio import iter_csv_files, read_dict_rows


@register_source
class NGSIM(TrajectorySource):
    name = "ngsim"
    summary = "10 Hz vehicle tracks at 4 US road sites; 11.85 M rows in one 1.5 GB CSV"
    dims = (2,)
    axes = "x = Global_X, y = Global_Y (State Plane feet)"

    OPTIONS = {
        "pattern": Option(str, "*.csv", "glob when root is a directory"),
        "location": Option(str, "", "keep only this site: us-101, i-80, lankershim, peachtree"),
        "grouping": Option(str, "buffered", "buffered | external | contiguous"),
        "chunk_rows": Option(int, 1_000_000, "external: rows held in memory per run"),
    }

    def __init__(
        self,
        root: str,
        *,
        pattern: str = "*.csv",
        location: str = "",
        grouping: str = "buffered",
        chunk_rows: int = 1_000_000,
    ) -> None:
        super().__init__(root)
        self.pattern = pattern
        self.location = location.strip().lower()
        self.grouping = grouping
        self.chunk_rows = chunk_rows

    def tracks(self) -> Iterator[Track]:
        paths = iter_csv_files(Path(self.root), self.pattern)

        def rows() -> Iterator[tuple[str, str, str, str]]:
            for _path, head, row in read_dict_rows(paths):
                loc = (_cell(head, row, "location") or "unknown").strip()
                if self.location and loc.lower() != self.location:
                    continue
                vid = _cell(head, row, "vehicle_id")
                x = _cell(head, row, "global_x")
                y = _cell(head, row, "global_y")
                if not vid or not x or not y:
                    continue
                section = (_cell(head, row, "section_id") or "").strip()
                key = f"{loc}/{vid}" if not section or section.upper() == "NA" \
                    else f"{loc}/{vid}s{section}"
                yield key, (_cell(head, row, "global_time") or ""), x.strip(), y.strip()

        groups = group_by(
            self.grouping,
            rows(),
            key=lambda r: r[0],
            sort_key=lambda r: r[1],
            chunk_rows=self.chunk_rows,
        )
        for key, group in groups:
            group.sort(key=lambda r: _int(r[1]))
            yield Track(id=f"ngsim/{key}", points=[(r[2], r[3], None) for r in group])


def _cell(head: dict[str, int], row: list[str], name: str) -> str | None:
    i = head.get(name)
    if i is None or i >= len(row):
        return None
    return row[i]


def _int(text: str) -> int:
    try:
        return int(text)
    except ValueError:
        return 0
