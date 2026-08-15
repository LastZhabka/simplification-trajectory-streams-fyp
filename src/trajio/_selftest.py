"""Every parser against a synthetic fixture in the dataset's real on-disk format.

No downloads needed:  python -m trajio selftest
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

from . import export, open_source

_FAILS: list[str] = []


def check(ok: bool, msg: str) -> None:
    if not ok:
        _FAILS.append(msg)


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


# -- fixtures -------------------------------------------------------------------------


def make_geolife(root: Path) -> None:
    d = root / "Data" / "000" / "Trajectory"
    d.mkdir(parents=True)
    (d / "20080412063127.plt").write_text(
        "\n".join(["hdr1", "hdr2", "hdr3", "hdr4", "hdr5", "hdr6",
                   "39.984702,116.318417,0,492,39550.2718402778,2008-04-12,06:31:27",
                   "39.984683,116.318450,0,492,39550.2718518519,2008-04-12,06:31:28",
                   "39.984686,116.318417,0,-777,39550.2718634259,2008-04-12,06:31:29"]) + "\n",
        encoding="utf-8")
    d2 = root / "Data" / "001" / "Trajectory"
    d2.mkdir(parents=True)
    (d2 / "20080412070000.plt").write_text(
        "\n".join(["hdr1", "hdr2", "hdr3", "hdr4", "hdr5", "hdr6",
                   "40.000000,116.400000,0,100,39550.29,2008-04-12,07:00:00",
                   "40.000100,116.400100,0,101,39550.30,2008-04-12,07:00:01"]) + "\n",
        encoding="utf-8")


def make_mopsi(root: Path) -> None:
    d = root / "routes" / "12"
    d.mkdir(parents=True)
    (d / "1216040000000").write_text(
        "62.601000 29.762000 1216040000000 85.0\n"
        "62.601100 29.762200 1216040005000 -1.0\n"
        "62.601200 29.762400 1216040010000 87.0\n", encoding="utf-8")
    d2 = root / "routes" / "13"
    d2.mkdir(parents=True)
    (d2 / "1216050000000").write_text(
        "62.700000 29.800000 1216050000000 90.0\n"
        "62.700100 29.800100 1216050005000 91.0\n", encoding="utf-8")


def make_ngsim(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    head = ("Vehicle_ID,Frame_ID,Global_Time,Local_X,Local_Y,Global_X,Global_Y,"
            "Lane_ID,Section_ID,Location")
    rows = [
        # interleaved, as in the real export
        "1,1,1113433135100,16.884,48.213,6451137.641,1873344.962,2,NA,us-101",
        "2,1,1113433135100,20.000,60.000,6451140.000,1873350.000,3,NA,us-101",
        "1,2,1113433135200,16.938,52.199,6451140.329,1873348.756,2,NA,us-101",
        "2,2,1113433135200,20.100,64.000,6451141.000,1873354.000,3,NA,us-101",
        # peachtree: ONE Vehicle_ID, two sections, same instant -- two distinct vehicles
        "7,1,1163483100000,-7.481,100.159,2225000.1,1375000.2,1,1,peachtree",
        "7,1,1163483100000,-14.121,700.411,2225600.9,1375000.5,1,2,peachtree",
        "7,2,1163483100100,-7.500,103.000,2225001.4,1375003.1,1,1,peachtree",
        "7,2,1163483100100,-14.200,703.000,2225601.7,1375003.4,1,2,peachtree",
    ]
    (root / "ngsim.csv").write_text(head + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


def make_mot(root: Path) -> None:
    flat = root / "MOT17"
    flat.mkdir(parents=True)
    (flat / "MOT17-02-FRCNN.txt").write_text(
        "1,1,100,200,50,80,1,1,1\n2,1,102,201,50,80,1,1,1\n1,2,600,300,40,60,1,1,1\n",
        encoding="utf-8")
    (flat / "MOT17-04-FRCNN.txt").write_text(
        "1,1,10,20,4,6,1,1,1\n2,1,12,22,4,6,1,1,1\n", encoding="utf-8")


# -- checks ---------------------------------------------------------------------------


def check_geolife(root: Path, out: Path) -> None:
    tracks = list(open_source("geolife", str(root)))
    check(len(tracks) == 2, f"geolife: expected 2 files -> 2 tracks, got {len(tracks)}")
    t = tracks[0]
    # x = lon, y = lat, both verbatim
    check(t.points[0] == ("116.318417", "39.984702", "492"),
          f"geolife: wrong point {t.points[0]}")
    check(t.points[2][2] is None, "geolife: -777 altitude not treated as missing")
    check(not t.has_z(), "geolife: track with a -777 point must not count as 3D")
    check(tracks[1].has_z(), "geolife: track with full altitude should be 3D-capable")

    idx = export(open_source("geolife", str(root)), out, dims=3)
    check([r[3] for r in idx] == [2, 3], f"geolife: dims per file wrong: {[r[3] for r in idx]}")
    first = _read(out / "geolife-000001.json")
    check(first["dim"] == 2, f"geolife: dim {first['dim']}")
    check(first["points"][0] == [116.318417, 39.984702], f"geolife: point {first['points'][0]}")
    second = _read(out / "geolife-000002.json")
    check(second["dim"] == 3, f"geolife: 3D dim {second['dim']}")
    check(second["points"][0] == [116.4, 40.0, 100], f"geolife: 3D point {second['points'][0]}")
    # the number tokens must be the source text, not a re-rendering of it
    raw = (out / "geolife-000001.json").read_text(encoding="utf-8")
    check("[116.318417, 39.984702]" in raw, "geolife: coordinates not written verbatim")

    # 2008-04-12 06:31:27 UTC, converted -- the one field that cannot be copied
    check(first["t_unit"] == "unix_ms", f"geolife: t_unit {first.get('t_unit')!r}")
    check(first["t"] == [1207981887000, 1207981888000, 1207981889000],
          f"geolife: timestamps {first['t']}")
    check(len(second["t"]) == len(second["points"]), "geolife: t/points length mismatch")


def check_mopsi(root: Path, out: Path) -> None:
    tracks = list(open_source("mopsi", str(root)))
    check(len(tracks) == 2, f"mopsi: expected 2 routes, got {len(tracks)}")
    t = tracks[0]
    check(t.points[0] == ("29.762000", "62.601000", "85.0"), f"mopsi: wrong point {t.points[0]}")
    check(t.points[1][2] is None, "mopsi: -1.0 altitude not treated as missing")
    check(not t.has_z(), "mopsi: track with a -1.0 point must not count as 3D")
    idx = export(open_source("mopsi", str(root)), out, dims=3)
    check([r[3] for r in idx] == [2, 3], f"mopsi: dims per file wrong: {[r[3] for r in idx]}")

    doc = _read(out / "mopsi-000001.json")
    check(doc["t_unit"] == "unix_ms", f"mopsi: t_unit {doc.get('t_unit')!r}")
    check(doc["t"] == [1216040000000, 1216040005000, 1216040010000],
          f"mopsi: timestamps {doc['t']}")
    raw = (out / "mopsi-000001.json").read_text(encoding="utf-8")
    check("1216040000000" in raw and '"1216040000000"' not in raw,
          "mopsi: epoch-ms timestamp not written verbatim as a number")


def check_ngsim(root: Path, out: Path) -> None:
    us = list(open_source("ngsim", str(root), location="us-101"))
    check(len(us) == 2, f"ngsim: expected 2 us-101 vehicles, got {len(us)}")
    check(us[0].points[0] == ("6451137.641", "1873344.962", None),
          f"ngsim: Global_X/Y not used verbatim: {us[0].points[0]}")

    # the arterial trap: one Vehicle_ID, two sections, two simultaneous vehicles
    pt = list(open_source("ngsim", str(root), location="peachtree"))
    check(len(pt) == 2, f"ngsim: Section_ID not in the key -- got {len(pt)} track(s), want 2")
    check(all(len(t) == 2 for t in pt), f"ngsim: section track lengths {[len(t) for t in pt]}")

    export(open_source("ngsim", str(root), location="us-101"), out, dims=2)
    doc = _read(out / "ngsim-000001.json")
    check(doc["dim"] == 2 and doc["points"][0] == [6451137.641, 1873344.962],
          f"ngsim: json {doc['dim']}, {doc['points'][:1]}")
    check(doc["t_unit"] == "unix_ms", f"ngsim: t_unit {doc.get('t_unit')!r}")
    check(doc["t"] == [1113433135100, 1113433135200], f"ngsim: timestamps {doc['t']}")
    # peachtree's Global_Time is milliseconds on a non-Unix origin, so it must not claim one
    check(open_source("ngsim", str(root), location="peachtree").time_unit == "ms",
          "ngsim: peachtree must not be labelled unix_ms")
    check(open_source("ngsim", str(root)).time_unit == "ms",
          "ngsim: an all-sites export mixes epochs and must not claim unix_ms")
    raw = (out / "ngsim-000001.json").read_text(encoding="utf-8")
    check("1113433135100" in raw and '"1113433135100"' not in raw,
          "ngsim: Global_Time not written verbatim as a number")

    try:
        export(open_source("ngsim", str(root)), out, dims=3)
        _FAILS.append("ngsim: dims=3 should be rejected (no altitude)")
    except ValueError:
        pass


def check_mot(root: Path, out: Path) -> None:
    tracks = list(open_source("mot", str(root)))
    check(len(tracks) == 3, f"mot: expected 3 tracks across 2 sequences, got {len(tracks)}")
    seqs = {t.id.split("/")[1] for t in tracks}
    check(seqs == {"MOT17-02-FRCNN", "MOT17-04-FRCNN"},
          f"mot: sequence not taken from the file stem: {seqs}")
    t = next(t for t in tracks if t.id.endswith("/1") and "02" in t.id)
    check(t.points[0] == ("125", "240", None), f"mot: box centre wrong: {t.points[0]}")
    export(open_source("mot", str(root)), out, dims=2)
    check((out / "mot-000001.json").exists(), "mot: no output written")

    # no clock in the data: t is synthesised from the assumed frame rate, and says so
    doc = _read(out / "mot-000001.json")
    check(doc["t_unit"] == "ms", f"mot: t_unit {doc.get('t_unit')!r}")
    check(doc["t"] == [0, 33], f"mot: frames 1,2 at 30 fps should be 0,33; got {doc['t']}")

    slow = out / "fps20"
    export(open_source("mot", str(root), fps=20), slow, dims=2)
    check(_read(slow / "mot-000001.json")["t"] == [0, 50],
          "mot: fps option not honoured")
    check(_read(slow / "index.json").get("fps") == 20,
          "mot: index.json must record the frame rate that was assumed")


def check_index(out: Path) -> None:
    idx = _read(out / "index.json")
    check({"dataset", "count", "trajectories"} <= set(idx), f"index: keys {sorted(idx)}")
    check(idx["count"] == len(idx["trajectories"]) > 0, "index: no entries")
    check(set(idx["trajectories"][0]) == {"file", "name", "points", "dim"},
          f"index: entry keys {sorted(idx['trajectories'][0])}")


def run() -> int:
    _FAILS.clear()
    with tempfile.TemporaryDirectory(prefix="trajio-selftest-") as tmp:
        base = Path(tmp)
        for name, make in [("geolife", make_geolife), ("mopsi", make_mopsi),
                           ("ngsim", make_ngsim), ("mot", make_mot)]:
            make(base / name)

        for name, fn in [("geolife", check_geolife), ("mopsi", check_mopsi),
                         ("ngsim", check_ngsim), ("mot", check_mot)]:
            before = len(_FAILS)
            try:
                fn(base / name, base / "out" / name)
            except Exception as exc:
                _FAILS.append(f"{name}: raised {type(exc).__name__}: {exc}")
            print(f"  {name:8s} {'ok' if len(_FAILS) == before else 'FAILED'}")
        try:
            check_index(base / "out" / "mopsi")
            print("  index    ok")
        except Exception as exc:
            _FAILS.append(f"index: raised {type(exc).__name__}: {exc}")
            print("  index    FAILED")

    if _FAILS:
        print(f"\n{len(_FAILS)} failure(s):")
        for m in _FAILS:
            print(f"  - {m}")
        return 1
    print("\nall checks passed")
    return 0
