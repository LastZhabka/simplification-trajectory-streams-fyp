# data/downloads/

Four public trajectory datasets, exactly as fetched — each in its own on-disk format,
untouched. Nothing here is generated, and nothing here is committed.

`trajio` reads these and writes `data/trajectories/`, our format: one JSON document per
trajectory, coordinates copied from the source files verbatim. No projection, no unit
conversion, no filtering, no splitting. Standard library only, Python 3.10 or later.

Full detail: [`../../docs/datasets.md`](../../docs/datasets.md) for the datasets and their
formats, [`../../docs/trajio.md`](../../docs/trajio.md) for the reader, and
[`../README.md`](../README.md) for what our format is.

## Converting

From the repository root, with `src/` on the import path:

```powershell
$env:PYTHONPATH = "$PWD\src"

python -m trajio sources                 # what is available
python -m trajio describe geolife        # options, axes, format notes
python -m trajio selftest                # every parser, no downloads needed

python -m trajio export --source mopsi --root data\downloads\mopsi --dims 2 `
                        --out data\trajectories\mopsi
python -m trajio stats  --source mopsi --root data\downloads\mopsi

# MOT has no clock: state the frame rate its timestamps are synthesised from
python -m trajio export --source mot --root data\downloads\mot\dataset\DanceTrack `
                        --opt fps=20 --dims 2 --out data\trajectories\mot-dancetrack
```

Export writes `<out>/<dataset>-<ord>.json` plus `<out>/index.json`.

## The datasets

| Dataset | `--source` | dims | x, y | time | Trajectory = |
|---|---|---|---|---|---|
| GeoLife | `geolife` | 2, 3 | lon, lat (deg); z = altitude in **feet** | date+time, read as UTC | one `.plt` file |
| Mopsi | `mopsi` | 2, 3 | lon, lat (deg); z = altitude in **metres** | epoch ms, verbatim | one route file |
| NGSIM | `ngsim` | 2 | `Global_X`, `Global_Y` (State Plane **feet**) | ms, verbatim; epoch only on 3 of 4 sites | one (site, vehicle, section) |
| MOT | `mot` | 2 | bounding-box centre (**pixels**) | **synthesised** from `--opt fps=` | one (sequence, track id) |

`--dims 3` is rejected on 2D-only datasets. Where a dataset supports 3, a trajectory with a
missing altitude on any point is written 2D instead. GeoLife's `-777` and Mopsi's `-1.0`
mean "no reading".

## Download

```powershell
# GeoLife  313 164 406 B
curl.exe -L -C - -o data\downloads\geolife\geolife.zip "https://download.microsoft.com/download/F/4/8/F4894AA5-FDBC-481E-9285-D5F8C4C4F039/Geolife%20Trajectories%201.3.zip"

# Mopsi  81 091 051 B
curl.exe -L -C - -o data\downloads\mopsi\MopsiRoutes2014.zip "http://cs.uef.fi/mopsi/routes/dataset/MopsiRoutes2014.zip"

# NGSIM  1 532 183 381 B, ~38 min (no Content-Length, generated on the fly)
curl.exe -L --retry 3 -o data\downloads\ngsim\ngsim.csv "https://data.transportation.gov/api/views/8ect-6jqj/rows.csv?accessType=DOWNLOAD"

# MOT  31 MB total, annotations only (motchallenge.net is unreachable from here)
python -m gdown 1AjiqAP2AGR_Qk8M0t2y388LvH7EDpZok -O MOT17.zip
python -m gdown 1xkpnUaM54dzwBfakVUQMlG5qaQdyVQZc -O MOT20.zip
python -m gdown 1mOb1g-ptPX9h9Djlj-xVvn7MdEerqWLX -O DanceTrack.zip
```

Unzip with `[System.IO.Compression.ZipFile]::ExtractToDirectory`, since `Expand-Archive`
takes minutes on GeoLife's 18 670 files. Archives are left compressed where the reader can
stream them.

Roots after extraction — this is what `--root` expects:

```
data\downloads\mopsi                                routes\<user>\<start_ms>  (no extension)
data\downloads\geolife\Geolife Trajectories 1.3     Data\<uuu>\Trajectory\*.plt
data\downloads\ngsim                                ngsim.csv
data\downloads\mot\dataset\MOT17 | MOT20 | DanceTrack
```

## Notes

- **NGSIM `Vehicle_ID` is not a trajectory id.** On the arterial sites (lankershim,
  peachtree) it is reused by distinct vehicles at the same instant, so the grouping key
  includes `Section_ID`. Without that, two cars hundreds of metres apart become one track.
- **NGSIM `Global_Time` is not on one epoch.** us-101, i-80 and lankershim decode to their
  real 2005 collection dates; peachtree's values are 1000× smaller and decode to January
  1970, though they still step by 100 per 10 Hz sample. Peachtree is exported as `ms`, not
  `unix_ms`. Sampling intervals are correct everywhere, so SED is unaffected.
- **MOT `x,y` are computed** (`left + width/2`), and its timestamps are synthesised from an
  assumed frame rate — the only values not copied from the files. The flattened archives
  have no `seqinfo.ini`, and MOT17's sequences genuinely run at 14, 25 and 30 fps, so one
  `--opt fps=` per export is wrong for some of them. `index.json` records what was used.
- **GeoLife's timestamps are converted**, not copied: `2008-10-23,02:53:04` read as UTC.
  The dataset states no zone and the traces are mostly Beijing, so absolute times may be
  offset; intervals are not.
- **NGSIM ships as one 1.5 GB CSV** and GeoLife/Mopsi as their own line formats. Those are
  the datasets' formats and are read as they are; the CSV in this repository is all
  *input*, never output.
- **δ is in the file's own units**: degrees for GeoLife/Mopsi, feet for NGSIM, pixels for
  MOT. A δ that works on one dataset means nothing on another.
