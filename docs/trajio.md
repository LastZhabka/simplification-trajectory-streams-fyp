# trajio

A reader for public trajectory datasets. It parses each dataset's native format and writes
one JSON document per trajectory, with the coordinate values copied from the source files
unchanged.

Standard library only, Python 3.10 or later. No NumPy, no pandas. Lives in `src/trajio/`;
put `src/` on `PYTHONPATH` to run it. Paths below are relative to the repository root.

## What it does and does not do

It does: locate the records, group them into trajectories, select which stored fields become
`x`, `y` and `z`, and write them out.

It does not: project coordinates, convert units, resample, split on time gaps, remove
duplicate points, or reject outliers. The numbers in the output are the numbers in the input
file.

There are three exceptions, all of them documented rather than silent:

| Exception | Source | What happens |
|---|---|---|
| `x, y` are computed | MOT | the bounding box centre, `left + width/2` |
| the timestamp is converted | GeoLife | `2008-10-23,02:53:04` parsed as UTC into epoch ms |
| the timestamp is **synthesised** | MOT | `(frame - 1) * 1000 / fps`, from an assumed frame rate |

Otherwise this is deliberate. Any filtering is a modelling decision that belongs to the
experiment, not to the reader, and a reader that silently alters values makes results hard
to trace.

## Data model

```
Point = (x: str, y: str, z: str | None)
Track = (id: str, points: list[Point], times: list[str])
```

Coordinates and timestamps are carried as the source text, not as floats, so nothing is
re-rendered and no float round-trip artefacts appear in the output. `z` is `None` where the
dataset has no altitude for that point. `times` is empty where the source has no clock, and
otherwise holds one value per point.

A source yields `Track` objects. `export()` writes them. That is the whole pipeline.

## Time

The SED-based simplification algorithms (DOTS, SQUISH, TDTR, OPW-TR) need a timestamp per
point: they use it as the interpolation ratio `k = (t_i - t_prev) / (t_next - t_prev)`.
Douglas-Peucker and Frechet-bounded simplification do not. Time is carried so both kinds of
algorithm can run on the same files.

Each source declares what its timestamps mean:

| Source | `time_unit` | Where it comes from |
|---|---|---|
| Mopsi | `unix_ms` | field 3, **already epoch ms** -- copied verbatim |
| NGSIM | `unix_ms` / `ms` | `Global_Time`, **already ms** -- copied verbatim; see below |
| GeoLife | `unix_ms` | `date` + `time` columns, parsed as UTC |
| MOT | `ms` | synthesised from `--opt fps=`; milliseconds from frame 1 |

`unix_ms` is an absolute epoch; `ms` is milliseconds on an origin that is not the Unix epoch.
MOT is deliberately **not** labelled `unix_ms` -- its origin is invented, and claiming a wall
clock it does not have would be worse than saying so.

Three caveats worth carrying into any result:

- **NGSIM's epoch is not uniform across its sites.** On us-101, i-80 and lankershim
  `Global_Time` is 13 digits and decodes to the real collection dates. On **peachtree** it is
  10 digits and decodes to January 1970, while still stepping by 100 per 10 Hz sample -- so
  the unit is milliseconds there too, but the origin is not Unix. Peachtree exports, and any
  export that does not pin `--opt location=` to a single Unix site, are labelled `ms`.
  Intervals are unaffected, so SED is unaffected.

- **MOT's frame rate is a guess.** The flattened archives carry no `seqinfo.ini`. DanceTrack
  is 20 fps; MOT17 is genuinely mixed (14, 25 and 30 across its sequences), so a single
  `--opt fps=` per export run is wrong for some of them, and every SED computed from those
  is distorted in proportion. `index.json` records the value used.
- **GeoLife's zone is assumed.** The dataset does not state one and the traces are mostly
  Beijing. Reading them as UTC may offset absolute times by hours; intervals, which is all
  SED uses, are unaffected.

## Dimensionality

Each source declares which dimensionalities it can supply:

```python
class GeoLife(TrajectorySource):
    dims = (2, 3)
    axes = "x = lon, y = lat, z = altitude (feet)"
```

`--dims 3` on a source that declares `(2,)` is rejected with an error rather than silently
writing 2D. Where a source declares `(2, 3)`, the decision is made per trajectory: if every
point in the track carries a `z`, the file is written `"dim": 3`; if any point lacks one,
the file is written `"dim": 2`. Both cases can appear in a single export, and `index.json`
records which is which.

Altitude sentinels are the reason this matters. GeoLife writes `-777` and Mopsi writes `-1.0`
to mean "no reading". These are mapped to `None`, so they never enter the `z` column as if
they were measurements.

## Grouping

Datasets stored as one file per trajectory need no grouping. NGSIM is a flat table of 11.85
million rows with the trajectories interleaved, so its rows must be collected by key.
`grouping.py` provides three strategies:

| Strategy | Memory | Use for |
|---|---|---|
| `contiguous` | O(1) | rows already ordered by object |
| `buffered` | O(objects) | interleaved, object count fits in memory (the NGSIM default) |
| `external` | O(chunk) | interleaved and too large; merge sort through temporary files |

For NGSIM the key is `(Location, Vehicle_ID, Section_ID)`. `Vehicle_ID` alone is not a
trajectory identifier: on the arterial sites it is reused by different vehicles at the same
instant, so grouping on it alone produces tracks that jump between two cars hundreds of
metres apart.

## Command line

```powershell
python -m trajio sources                  # datasets, their dims, one-line summaries
python -m trajio describe geolife         # options, axes, format notes
python -m trajio selftest                 # all parsers against fixtures, no downloads
python -m trajio stats   --source mopsi --root data\downloads\mopsi
python -m trajio export  --source mopsi --root data\downloads\mopsi --dims 2 `
                         --out data\trajectories\mopsi
```

Shared flags: `--opt KEY=VALUE` (repeatable, source specific) and `--limit N`.

`stats` prints trajectory and point counts, length percentiles, how many tracks have a usable
altitude, and a length histogram. It reads without writing anything.

## Output

`export` writes `<out>/<dataset>-<ord>.json`, numbered from 1 in the order the source
yields them, plus `<out>/index.json`:

```json
{
  "dim": 2,
  "name": "geolife/000/20081023025304",
  "t_unit": "unix_ms",
  "t": [1224730384000, 1224730390000],
  "points": [
    [116.318417, 39.984702],
    [116.31845, 39.984683]
  ]
}
```

`t` is present only when the source has a clock, and then holds exactly one value per point.

```json
{"dataset": "geolife", "count": 18670, "t_unit": "unix_ms",
 "trajectories": [
   {"file": "geolife-000001.json", "name": "geolife/000/20081023025304",
    "points": 1047, "dim": 2}
 ]}
```

The index also carries `"fps"` when the source synthesised its timestamps, so the assumption
that produced a directory of files is recorded alongside them.

`name` traces a file back to the record it came from, which is what makes a result
reproducible without re-running the export.

Coordinates and timestamps are written as the source text, placed straight into the JSON
number slot, so nothing is re-rendered on the way out. A value whose text is not legal JSON
number syntax — a leading `+` or `.`, a leading zero — is the one exception: it goes through
a double, which preserves the value but not the spelling. None of the four datasets produces
one.

Epoch milliseconds survive the round trip through the C++ reader's `double` exactly: 1.2e12
is far below 2^53.

The C++ reader in `src/algo/io/trajectory.hpp` takes these files directly. `dim` comes from
the field when present and is inferred from the first point otherwise; either way every
point must match it, and a `t` that is present must have one entry per point.

```cpp
auto doc = ssk::io::read_trajectory_file("data/trajectories/mopsi/mopsi-000002.json");
```

Note that a tolerance `delta` is expressed in the file's own units, which differ per
dataset. See [datasets.md](datasets.md).

## Layout

```
src/trajio/
  core.py        Point, Track, TrajectorySource, the registry
  export.py      JSON writer and index
  stats.py       counts, percentiles, length histogram
  grouping.py    contiguous / buffered / external merge sort
  cli.py         python -m trajio ...
  _selftest.py   fixtures in each dataset's real on-disk format
  sources/       geolife.py, mopsi.py, ngsim.py, mot.py, _csvio.py
```

`_csvio.py` streams CSV out of `.zip`, `.gz` and `.zst` archives without unpacking them.

## Testing

`python -m trajio selftest` builds a synthetic fixture for every source in that dataset's
real on-disk format, then asserts the parsed result. It needs no downloads and runs in about
a second.

The fixtures encode the traps, so a regression fails the test rather than appearing quietly
in the data:

- GeoLife six line header, `-777` altitude, and the lon/lat column order
- Mopsi files with no extension, and the `-1.0` altitude sentinel
- NGSIM rows interleaved across vehicles, and one `Vehicle_ID` used by two vehicles in
  different sections at the same instant
- MOT flattened sequence files, where the sequence name comes from the file stem

It also checks that `--dims 3` is rejected on 2D-only sources, that a track with a missing
`z` falls back to 2D, that coordinates and epoch-ms timestamps reach the file verbatim, that
GeoLife's clock converts to the expected epoch and MOT's synthesises at the requested frame
rate, and that `index.json` is well formed.

## Adding a dataset

1. Subclass `TrajectorySource` in `trajio/sources/`, setting `name`, `summary`, `dims`,
   `axes` and `OPTIONS`.
2. Implement `tracks()`, yielding `Track(id, [(x, y, z_or_None), ...])` with coordinates as
   source strings.
3. Decorate with `@register_source` and import it in `sources/__init__.py`.
4. Add a fixture and checks to `_selftest.py`.

No projection, no filtering and no dimension logic belong in a source.
