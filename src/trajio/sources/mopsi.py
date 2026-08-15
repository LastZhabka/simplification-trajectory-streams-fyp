"""Mopsi Routes 2014 (University of Eastern Finland, SIPU).

    routes/<user 1..51>/<start_timestamp_ms>       one route per file, NO extension

One point per line, whitespace separated: `lat lon timestamp_ms altitude_m`.

x = lon, y = lat, z = altitude **in metres, as stored**. `-1.0` means "no reading" (13 % of
points, undocumented), so those points have no z and the track is written 2D.

The third field is already a Unix timestamp in milliseconds, so it is copied verbatim like
the coordinates. A line without it is not a point and is skipped.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from ..core import Option, Track, TrajectorySource, register_source

_ALT_MISSING = "-1.0"
# Route files have no extension, so the default glob would otherwise parse the downloaded
# archive itself as a route.
_NOT_ROUTES = {".zip", ".gz", ".7z", ".rar", ".tar", ".pdf", ".md", ".csv"}


@register_source
class Mopsi(TrajectorySource):
    name = "mopsi"
    summary = "6 779 routes, 51 users, Finland 2008-2014; has altitude (metres)"
    dims = (2, 3)
    axes = "x = lon, y = lat, z = altitude (metres)"
    time_unit = "unix_ms"

    OPTIONS = {
        "pattern": Option(str, "**/*", "glob for route files, relative to root"),
        "users": Option(str, "", "comma-separated user directory names to keep (default all)"),
    }

    def __init__(self, root: str, *, pattern: str = "**/*", users: str = "") -> None:
        super().__init__(root)
        self.pattern = pattern
        self.users = {u.strip() for u in users.split(",") if u.strip()}

    def tracks(self) -> Iterator[Track]:
        root = Path(root_str := self.root)
        files = sorted(
            p for p in root.glob(self.pattern)
            if p.is_file() and p.suffix.lower() not in _NOT_ROUTES
        )
        if not files:
            raise FileNotFoundError(f"no route files matching {self.pattern!r} under {root_str}")
        for path in files:
            if self.users and path.parent.name not in self.users:
                continue
            rows = list(self._read(path))
            if rows:
                yield Track(
                    id=f"mopsi/{path.relative_to(root).as_posix()}",
                    points=[point for point, _t in rows],
                    times=[t for _point, t in rows],
                )

    @staticmethod
    def _read(path: Path) -> Iterator[tuple[tuple[str, str, str | None], str]]:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 3:
                    continue
                alt = parts[3] if len(parts) >= 4 else None
                yield (parts[1], parts[0], None if alt == _ALT_MISSING else alt), parts[2]
