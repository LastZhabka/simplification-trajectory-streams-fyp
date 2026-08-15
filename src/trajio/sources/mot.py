"""MOT-format video object tracks (MOT17, MOT20, DanceTrack).

The archives from the Trajectory Simplify Benchmark are flattened -- one `.txt` per
sequence, no `gt/` directory:

    MOT17/MOT17-02-FRCNN.txt        frame, id, bb_left, bb_top, bb_width, bb_height, ...

One trajectory per `(sequence, track id)`. x, y are the **box centre** in pixels --
`left + width/2`, `top + height/2` -- so unlike every other source these two numbers are
computed, not copied from the file. No z.

There is no clock here, only a frame number, so the timestamp is **synthesised**:
`(frame - 1) * 1000 / fps`, milliseconds from the first frame of the sequence. `fps` is an
option rather than a constant because it is a guess -- the flattened archives have no
`seqinfo.ini`. DanceTrack is 20 fps; MOT17 is genuinely mixed (14, 25 and 30 across its
sequences), so one value per export run is wrong for some of them and any SED computed from
those is distorted accordingly. `index.json` records the value that was used.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from ..core import Option, Track, TrajectorySource, fmt, register_source


@register_source
class MOTTracks(TrajectorySource):
    name = "mot"
    summary = "video object tracks (MOT17/MOT20/DanceTrack); box centre in pixels"
    dims = (2,)
    axes = "x, y = bounding-box centre (pixels)"
    time_unit = "ms"

    OPTIONS = {
        "pattern": Option(str, "**/*.txt", "glob for annotation files; **/gt/gt.txt for MOTChallenge"),
        "fps": Option(float, 30.0, "frame rate assumed when synthesising timestamps; 20 for DanceTrack"),
    }

    def __init__(self, root: str, *, pattern: str = "**/*.txt", fps: float = 30.0) -> None:
        super().__init__(root)
        self.pattern = pattern
        if fps <= 0:
            raise ValueError(f"fps must be positive, got {fps}")
        self.fps = fps

    def tracks(self) -> Iterator[Track]:
        root = Path(self.root)
        files = sorted(p for p in root.glob(self.pattern) if p.is_file())
        if not files:
            raise FileNotFoundError(f"no annotation files matching {self.pattern!r} under {root}")
        for path in files:
            # MOTChallenge puts the sequence name on the grandparent directory; the
            # flattened archives put it on the file itself. Confusing the two collapses
            # every sequence into one name and collides their track ids.
            seq = path.parent.parent.name if path.parent.name in {"gt", "det"} else path.stem
            per_track: dict[str, list[tuple[int, tuple[str, str, str | None]]]] = {}
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) < 6:
                        continue
                    try:
                        frame = int(float(parts[0]))
                        left, top, width, height = (float(v) for v in parts[2:6])
                    except ValueError:
                        continue
                    point = (fmt(left + width / 2.0), fmt(top + height / 2.0), None)
                    per_track.setdefault(parts[1], []).append((frame, point))
            for track_id, rows in per_track.items():
                rows.sort(key=lambda r: r[0])
                yield Track(
                    id=f"mot/{seq}/{track_id}",
                    points=[p for _f, p in rows],
                    times=[str(round((f - 1) * 1000.0 / self.fps)) for f, _p in rows],
                )
