# Datasets

Four public trajectory datasets, how they were obtained, what is taken from each, and how
the results are stored.

All figures below were measured on this machine, not quoted from the dataset pages. Where
the two disagree, the discrepancy is noted.

All paths in this document are relative to the repository root, which is where the commands
are meant to be run.

## Summary

| Dataset | Source | Download | Trajectory unit | dims |
|---|---|---|---|---|
| GeoLife | Microsoft Research Asia | 313 164 406 B, 33 s | one `.plt` file | 2, 3 |
| Mopsi | University of Eastern Finland | 81 091 051 B, 59 s | one route file | 2, 3 |
| NGSIM | US DOT / FHWA | 1 532 183 381 B, 38 min | one (site, vehicle, section) | 2 |
| MOT | Trajectory Simplify Benchmark | 30 881 984 B, about 5 s | one (sequence, track id) | 2 |

Total downloaded: about 1.96 GB.

## GeoLife GPS Trajectories

GPS traces collected by Microsoft Research Asia between April 2007 and August 2012, mostly in
Beijing. Research use licence bundled in the archive. Cite Zheng et al., WWW 2009, UbiComp
2008, and IEEE Data Engineering Bulletin 2010.

### Download

```powershell
curl.exe -L -C - -o data\downloads\geolife\geolife.zip "https://download.microsoft.com/download/F/4/8/F4894AA5-FDBC-481E-9285-D5F8C4C4F039/Geolife%20Trajectories%201.3.zip"
Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.IO.Compression.ZipFile]::ExtractToDirectory("data\downloads\geolife\geolife.zip", "data\downloads\geolife")
```

Use `ZipFile` rather than `Expand-Archive`, which takes minutes on 18 670 small files.

### Layout and format

```
data/downloads/geolife/Geolife Trajectories 1.3/Data/<uuu>/Trajectory/<yyyyMMddHHmmss>.plt
```

Six header lines, then one point per line:

```
39.984702,116.318417,0,492,39744.1201851852,2008-10-23,02:53:04
lat,       lon,       0, alt, days since 1899-12-30, date, time
```

### What is taken

`x = lon` (field 1), `y = lat` (field 0), `z = altitude` (field 3). The altitude is in feet
as stored and is not converted. The value `-777` means "no reading" and is treated as a
missing `z`, so a track containing one is written 2D.

### Version discrepancy

The archive contains 182 user directories, 18 670 `.plt` files, 69 `labels.txt` files and
24 876 978 points, and bundles `User Guide-1.3.pdf`. The frequently cited figures of 17 621
trajectories and 73 labelled users come from the v1.2 guide and do not describe this
download. Report counts as measured.

## Mopsi Routes 2014

GPS routes collected by the SIPU group at the University of Eastern Finland between 2008 and
2014, mostly around Joensuu. Citation required: Mariescu-Istodor and Fränti, *Grid-based
method for GPS route analysis for retrieval*, ACM TSAS 3(3), 2017.

### Download

```powershell
curl.exe -L -C - -o data\downloads\mopsi\MopsiRoutes2014.zip "http://cs.uef.fi/mopsi/routes/dataset/MopsiRoutes2014.zip"
```

The archive is served from the dataset page's own directory (`/mopsi/routes/dataset/`), not
from the parent path.

### Layout and format

```
data/downloads/mopsi/routes/<user 1..51>/<start_timestamp_ms>
```

One route per file, named after its start timestamp, with **no file extension**. One point
per line, whitespace separated:

```
62.61478 29.74249 1216461503656 -1.0
lat      lon      timestamp_ms   altitude_m
```

### What is taken

`x = lon`, `y = lat`, `z = altitude` in metres as stored. The value `-1.0` means "no reading"
and is treated as a missing `z`. This sentinel is not documented by the dataset page, which
states only that altitude is "not always available"; in the files every line carries all four
fields and about 13 percent of them hold `-1.0`.

Counts match the published figures exactly: 6 779 routes, 51 users, 7 850 387 points.

## NGSIM Vehicle Trajectories

Vehicle positions extracted from synchronised video at 10 Hz at four US road sites in 2005
and 2006. US Government public domain.

### Download

```powershell
curl.exe -L --retry 3 -o data\downloads\ngsim\ngsim.csv "https://data.transportation.gov/api/views/8ect-6jqj/rows.csv?accessType=DOWNLOAD"
```

The Socrata endpoint generates the export on request, sends no `Content-Length`, and held
about 650 kB/s, so the download took 38 minutes and reported no progress percentage.

### Layout and format

One flat CSV, 11 850 526 data rows, distributed across the sites as:

| Site | Rows | Type |
|---|---:|---|
| `us-101` | 4 802 933 | freeway |
| `i-80` | 4 566 387 | freeway |
| `lankershim` | 1 607 319 | arterial |
| `peachtree` | 873 887 | arterial |

Column names use inconsistent capitalisation (`v_length` next to `v_Width`), so header
matching is case insensitive.

### What is taken

`x = Global_X`, `y = Global_Y`, State Plane NAD83 feet as stored. No altitude, so this
source is 2D only.

### Grouping

There is no file per trajectory, so rows are grouped by `(Location, Vehicle_ID, Section_ID)`.

`Vehicle_ID` alone is not a trajectory identifier. There are only 8 899 distinct
`(Location, Vehicle_ID)` pairs for 11.85 million rows, and 58 percent of them span more
wall-clock time than 10 Hz sampling allows. Two separate causes:

- On the freeway sites, `Section_ID` is `NA` and the identifier is reused across the
  15-minute recording extracts.
- On the arterial sites, it is additionally reused across road sections by *different
  vehicles at the same instant*. 98 percent (lankershim) and 89 percent (peachtree) of
  identifiers appearing in more than one section have overlapping time ranges.

An example row pair, one identifier, one timestamp, two positions 135 m apart:

```
peachtree  Vehicle_ID=701  t=1163483100
  A: Local_X=-14.121  Local_Y=542.411  Section_ID=2
  B: Local_X= -7.481  Local_Y=100.159  Section_ID=1
```

Including `Section_ID` separates them. It has no effect on the freeway sites, where the
column is `NA`.

## MOT video object tracks

Bounding box annotations from MOT17, MOT20 and DanceTrack, as distributed by the Trajectory
Simplify Benchmark repository. Cite the underlying tracking dataset used.

### Download

`motchallenge.net` was unreachable from this machine (connection timeout on every request,
including plain range GETs), so the repository's Google Drive copies were used. They contain
annotations only, no images, which is 31 MB rather than about 5 GB.

```powershell
git clone --depth 1 https://github.com/foreverYoungGitHub/trajectory-simplify-benchmark
python -m gdown 1AjiqAP2AGR_Qk8M0t2y388LvH7EDpZok -O MOT17.zip       # 3.97 MB
python -m gdown 1xkpnUaM54dzwBfakVUQMlG5qaQdyVQZc -O MOT20.zip       # 11.5 MB
python -m gdown 1mOb1g-ptPX9h9Djlj-xVvn7MdEerqWLX -O DanceTrack.zip  # 4.93 MB
```

### Layout and format

The archives are flattened, one file per sequence, with no `gt/` directory and no
`seqinfo.ini`:

```
data/downloads/mot/dataset/MOT17/MOT17-02-FRCNN.txt        21 sequences
data/downloads/mot/dataset/MOT20/MOT20-01.txt               4 sequences
data/downloads/mot/dataset/DanceTrack/dancetrack0001.txt   65 sequences
```

Comma separated, no header:

```
frame, id, bb_left, bb_top, bb_width, bb_height, conf, class, visibility
```

The reader takes the sequence name from the file stem for this layout, and from the
grandparent directory for the standard MOTChallenge `<seq>/gt/gt.txt` layout. Using the
parent directory in the flattened case would collapse all 65 DanceTrack sequences into one
name and collide their track identifiers.

### What is taken

`x, y` is the bounding box centre in pixels: `left + width/2`, `top + height/2`. These are
the only values in any dataset that are computed rather than copied. No altitude, so this
source is 2D only.

## Storage

### Downloads

```
data/downloads/geolife/Geolife Trajectories 1.3/Data/<uuu>/Trajectory/*.plt
data/downloads/mopsi/routes/<user>/<start_ms>
data/downloads/ngsim/ngsim.csv
data/downloads/mot/dataset/{MOT17,MOT20,DanceTrack}/<seq>.txt
data/downloads/mot/trajectory-simplify-benchmark/    upstream baseline code
```

`data/downloads/` is git-ignored. Archives are left compressed where the reader can stream them.

### Exported trajectories

One JSON document per trajectory, plus an index, per dataset directory. `data/trajectories/` is
git-ignored.

```
data/trajectories/mopsi/mopsi-000001.json
data/trajectories/mopsi/index.json
```

File contents:

```json
{
  "dim": 2,
  "name": "mopsi/routes/1/1216461503656",
  "t_unit": "unix_ms",
  "t": [1216461503656, 1216461508749],
  "points": [
    [29.74249, 62.61478],
    [29.742665, 62.614714]
  ]
}
```

Index contents:

```json
{"dataset": "mopsi", "count": 6779, "t_unit": "unix_ms",
 "trajectories": [
   {"file": "mopsi-000001.json", "name": "mopsi/routes/1/1216461503656",
    "points": 47, "dim": 2}
 ]}
```

### Export commands

```powershell
python -m trajio export --source geolife --root "data\downloads\geolife\Geolife Trajectories 1.3" --dims 2 --out data\trajectories\geolife
python -m trajio export --source mopsi   --root data\downloads\mopsi --dims 2 --out data\trajectories\mopsi
python -m trajio export --source mot     --root data\downloads\mot\dataset\MOT17 --opt fps=30 --dims 2 --out data\trajectories\mot-mot17
python -m trajio export --source mot     --root data\downloads\mot\dataset\DanceTrack --opt fps=20 --dims 2 --out data\trajectories\mot-dancetrack
python -m trajio export --source ngsim   --root data\downloads\ngsim --opt location=us-101 --opt grouping=external --dims 2 --out data\trajectories\ngsim-us-101
```

NGSIM is exported per site, because the four sites are separate coordinate frames. MOT is
exported per archive, and `--opt fps=` has to be set per archive too: DanceTrack is 20 fps,
MOT17 and MOT20 are taken as 30. MOT17's sequences are actually mixed (14, 25, 30), so its
synthesised timestamps are approximate for some of them.

### Export size

Measured on the CSV export this replaced. The point and file counts are unchanged by the
move to JSON; the byte sizes are not — one point per line as `[x, y],` runs roughly 1.3 to
1.4 times the width of `x,y`, and the parallel `t` array adds roughly 16 bytes per point on
top, so budget about 2 GB in total rather than 972 MB.

| Directory | Files | Points | Size (as CSV) |
|---|---:|---:|---:|
| `data/trajectories/geolife` | 18 670 | 24 876 978 | 529.6 MB |
| `data/trajectories/mopsi` | 6 779 | 7 850 387 | 148.0 MB |
| `data/trajectories/ngsim-us-101` | 2 847 | 4 802 933 | 109.0 MB |
| `data/trajectories/ngsim-i-80` | 3 001 | 4 566 387 | 103.7 MB |
| `data/trajectories/ngsim-lankershim` | 6 712 | 1 607 319 | 36.8 MB |
| `data/trajectories/ngsim-peachtree` | 4 495 | 873 887 | 20.0 MB |
| `data/trajectories/mot-mot17` | 2 388 | 614 103 | 6.1 MB |
| `data/trajectories/mot-mot20` | 2 332 | 1 336 920 | 12.8 MB |
| `data/trajectories/mot-dancetrack` | 692 | 574 078 | 5.7 MB |
| **Total** | **47 916** | **47 102 992** | **972 MB** |

## Units

Coordinates are never converted, so each dataset keeps its own units and `delta` means
something different in each:

| Dataset | Units | Example working `delta` |
|---|---|---|
| GeoLife, Mopsi | degrees | 0.0005 |
| NGSIM | feet (State Plane) | 3 |
| MOT | pixels | 5 |

Verified against `simplify`, all three within the Fréchet bound:

| Input | delta | Result | Fréchet |
|---|---|---|---|
| `mopsi-000002.json` | 0.0005 | 103 to 2 points | 0.000595 within 0.0006 |
| `ngsim-000001.json` (us-101) | 3 | 990 to 12 points | 3.572 within 3.6 |
| `mot-000001.json` (DanceTrack) | 5 | 697 to 62 points | 5.953 within 6 |
