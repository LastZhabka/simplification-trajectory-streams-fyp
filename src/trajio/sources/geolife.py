"""GeoLife GPS Trajectories (Microsoft Research Asia), v1.3.

    Data/<uuu>/Trajectory/<yyyyMMddHHmmss>.plt      one trajectory per file

Six header lines, then `lat,lon,0,altitude_feet,days,date,time`.

x = lon, y = lat, z = altitude **in feet, as stored**. `-777` means "no reading", so those
points have no z and the track is written 2D.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from ..core import Option, Track, TrajectorySource, register_source

_ALT_MISSING = "-777"


@register_source
class GeoLife(TrajectorySource):
    name = "geolife"
    summary = "182 users, 18 670 trajectories, Beijing 2007-2012; has altitude (feet)"
    dims = (2, 3)
    axes = "x = lon, y = lat, z = altitude (feet)"

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
                points = list(self._read(plt))
                if points:
                    yield Track(id=f"geolife/{user_dir.name}/{plt.stem}", points=points)

    @staticmethod
    def _read(path: Path) -> Iterator[tuple[str, str, str | None]]:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for _ in range(6):  # fixed-size header
                if not fh.readline():
                    return
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) < 4:
                    continue
                alt = parts[3].strip()
                yield parts[1].strip(), parts[0].strip(), None if alt == _ALT_MISSING else alt
