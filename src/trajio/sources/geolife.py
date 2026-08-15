"""GeoLife GPS Trajectories (Microsoft Research Asia), v1.3.

    Data/<uuu>/Trajectory/<yyyyMMddHHmmss>.plt      one trajectory per file

Six header lines, then `lat,lon,0,altitude_feet,days,date,time`.

x = lon, y = lat, z = altitude **in feet, as stored**. `-777` means "no reading", so those
points have no z and the track is written 2D.

Unlike the coordinates, the timestamp cannot be copied verbatim: the file stores
`2008-10-23,02:53:04`, which has to be converted to Unix milliseconds. The clock is read as
UTC. GeoLife does not state a zone and the traces are mostly Beijing, so absolute values may
be offset by hours -- intervals, which is what SED uses, are unaffected either way.
"""

from __future__ import annotations

from collections.abc import Iterator
from datetime import datetime, timezone
from pathlib import Path

from ..core import Option, Track, TrajectorySource, register_source

_ALT_MISSING = "-777"


def _epoch_ms(date: str, clock: str) -> str:
    stamp = datetime.strptime(f"{date} {clock}", "%Y-%m-%d %H:%M:%S")
    return str(int(stamp.replace(tzinfo=timezone.utc).timestamp() * 1000))


@register_source
class GeoLife(TrajectorySource):
    name = "geolife"
    summary = "182 users, 18 670 trajectories, Beijing 2007-2012; has altitude (feet)"
    dims = (2, 3)
    axes = "x = lon, y = lat, z = altitude (feet)"
    time_unit = "unix_ms"

    OPTIONS = {
        "users": Option(str, "", "comma-separated user ids to keep, e.g. 000,010 (default all)"),
    }

    def __init__(self, root: str, *, users: str = "") -> None:
        super().__init__(root)
        self.users = {u.strip() for u in users.split(",") if u.strip()}
        self._data = self._find_data_dir(Path(root))

    @staticmethod
    def _find_data_dir(root: Path) -> Path:
        """Accept the zip root, the `Geolife Trajectories 1.3` dir, or `Data` itself."""
        for candidate in (root / "Data", root):
            if candidate.is_dir() and any(candidate.glob("*/Trajectory")):
                return candidate
        for nested in sorted(root.glob("*/Data")):
            if any(nested.glob("*/Trajectory")):
                return nested
        raise FileNotFoundError(
            f"no GeoLife user directories under {root}; expected <root>/Data/000/Trajectory/*.plt"
        )

    def tracks(self) -> Iterator[Track]:
        for user_dir in sorted(p for p in self._data.iterdir() if p.is_dir()):
            if self.users and user_dir.name not in self.users:
                continue
            for plt in sorted((user_dir / "Trajectory").glob("*.plt")):
                rows = list(self._read(plt))
                if rows:
                    yield Track(
                        id=f"geolife/{user_dir.name}/{plt.stem}",
                        points=[point for point, _t in rows],
                        times=[t for _point, t in rows],
                    )

    @staticmethod
    def _read(path: Path) -> Iterator[tuple[tuple[str, str, str | None], str]]:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for _ in range(6):  # fixed-size header
                if not fh.readline():
                    return
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) < 7:
                    continue
                alt = parts[3].strip()
                try:
                    stamp = _epoch_ms(parts[5].strip(), parts[6].strip())
                except ValueError:
                    continue
                point = (parts[1].strip(), parts[0].strip(),
                         None if alt == _ALT_MISSING else alt)
                yield point, stamp
