# From zero to Fréchet distances

Every step from a fresh clone to a measured error for every simplification. **Each step gives
the PowerShell and the POSIX shell form**, since the project is developed on Windows but
nothing in it is Windows-only — see *Portability* at the end for what actually differs.

Run everything from the repository root. Each step says roughly what it costs, because two of
them are measured in hours.

---

## 1. Toolchain

A C++20 compiler and CMake 3.20+. Everything under `src/algo` and `src/pipelines` is standard
library only, so there is nothing to fetch and the build works offline. Python 3.10+ is needed
for ingestion and plotting.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

`requirements.txt` is matplotlib and numpy, and they are needed **only** by `src/viz/`.
`trajio` deliberately uses nothing outside the standard library.

## 2. Build and test

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
cmake --build build -j
.\build\tests\ssk_tests.exe
```

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/ssk_tests
```

The generator is the only build difference: MinGW needs naming, the default works elsewhere.
Use `-G Ninja` on either platform if you have it.

63 tests, a few seconds. Run them before anything long — they cover the simplifiers against
their reference implementations, the Fréchet distance against three, and both pipelines'
document formats.

A fast check that ingestion works, needing no downloads at all:

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m trajio selftest
```

```sh
export PYTHONPATH=$PWD/src
python -m trajio selftest
```

## 3. Download the datasets

About **1.96 GB**, and the NGSIM download alone takes ~40 minutes. Full detail on what each
dataset contains is in [datasets.md](datasets.md).

`curl` will not create the directory it writes into, and only `data/downloads/` itself is in
the repository, so make the four first. MOT is fetched with `gdown`, which is not in
`requirements.txt` because nothing in the project imports it — install it for this step alone:

```powershell
mkdir data\downloads\geolife, data\downloads\mopsi, data\downloads\ngsim, data\downloads\mot
python -m pip install gdown
```

```sh
mkdir -p data/downloads/{geolife,mopsi,ngsim,mot}
python -m pip install gdown
```

```powershell
curl.exe -L -C - -o data\downloads\geolife\geolife.zip "https://download.microsoft.com/download/F/4/8/F4894AA5-FDBC-481E-9285-D5F8C4C4F039/Geolife%20Trajectories%201.3.zip"
curl.exe -L -C - -o data\downloads\mopsi\MopsiRoutes2014.zip "http://cs.uef.fi/mopsi/routes/dataset/MopsiRoutes2014.zip"
curl.exe -L --retry 3 -o data\downloads\ngsim\ngsim.csv "https://data.transportation.gov/api/views/8ect-6jqj/rows.csv?accessType=DOWNLOAD"
```

```sh
curl -L -C - -o data/downloads/geolife/geolife.zip "https://download.microsoft.com/download/F/4/8/F4894AA5-FDBC-481E-9285-D5F8C4C4F039/Geolife%20Trajectories%201.3.zip"
curl -L -C - -o data/downloads/mopsi/MopsiRoutes2014.zip "http://cs.uef.fi/mopsi/routes/dataset/MopsiRoutes2014.zip"
curl -L --retry 3 -o data/downloads/ngsim/ngsim.csv "https://data.transportation.gov/api/views/8ect-6jqj/rows.csv?accessType=DOWNLOAD"
```

`curl.exe` is spelled out on Windows because bare `curl` is a PowerShell alias for
`Invoke-WebRequest`, which does not take these flags.

MOT comes from the Trajectory Simplify Benchmark's Google Drive copies, because
`motchallenge.net` was unreachable — see [datasets.md](datasets.md). The three archives are
annotations only, 31 MB rather than about 5 GB, and unpack into `mot/dataset/`:

```powershell
python -m gdown 1AjiqAP2AGR_Qk8M0t2y388LvH7EDpZok -O data\downloads\mot\dataset\MOT17.zip
python -m gdown 1xkpnUaM54dzwBfakVUQMlG5qaQdyVQZc -O data\downloads\mot\dataset\MOT20.zip
python -m gdown 1mOb1g-ptPX9h9Djlj-xVvn7MdEerqWLX -O data\downloads\mot\dataset\DanceTrack.zip
```

```sh
mkdir -p data/downloads/mot/dataset
python -m gdown 1AjiqAP2AGR_Qk8M0t2y388LvH7EDpZok -O data/downloads/mot/dataset/MOT17.zip
python -m gdown 1xkpnUaM54dzwBfakVUQMlG5qaQdyVQZc -O data/downloads/mot/dataset/MOT20.zip
python -m gdown 1mOb1g-ptPX9h9Djlj-xVvn7MdEerqWLX -O data/downloads/mot/dataset/DanceTrack.zip
```

### Unpack

NGSIM is a single CSV and is read as-is. The other three are archives, and each must extract
to the directory the export commands in step 4 name — `Geolife Trajectories 1.3` under
`geolife/`, `routes/` under `mopsi/`, and `MOT17`, `MOT20`, `DanceTrack` under `mot/dataset/`:

```powershell
Add-Type -AssemblyName System.IO.Compression.FileSystem
[IO.Compression.ZipFile]::ExtractToDirectory("data\downloads\geolife\geolife.zip", "data\downloads\geolife")
[IO.Compression.ZipFile]::ExtractToDirectory("data\downloads\mopsi\MopsiRoutes2014.zip", "data\downloads\mopsi")
foreach ($z in "MOT17", "MOT20", "DanceTrack") {
  [IO.Compression.ZipFile]::ExtractToDirectory("data\downloads\mot\dataset\$z.zip", "data\downloads\mot\dataset")
}
```

```sh
unzip -q data/downloads/geolife/geolife.zip -d data/downloads/geolife
unzip -q data/downloads/mopsi/MopsiRoutes2014.zip -d data/downloads/mopsi
for z in MOT17 MOT20 DanceTrack; do
  unzip -q "data/downloads/mot/dataset/$z.zip" -d data/downloads/mot/dataset
done
```

Use `ZipFile` rather than `Expand-Archive` on Windows: GeoLife is 18 670 small files and
`Expand-Archive` takes minutes on it.

## 4. Export to the project's format

Turns each dataset's native format into one JSON document per trajectory. **~10 minutes for
all nine**, producing about 2.1 GB.

```powershell
$env:PYTHONPATH = "$PWD\src"

python -m trajio export --source mopsi   --root data\downloads\mopsi --dims 2 --out data\trajectories\mopsi
python -m trajio export --source geolife --root "data\downloads\geolife\Geolife Trajectories 1.3" --dims 2 --out data\trajectories\geolife
python -m trajio export --source mot     --root data\downloads\mot\dataset\MOT17      --opt fps=30 --dims 2 --out data\trajectories\mot-mot17
python -m trajio export --source mot     --root data\downloads\mot\dataset\MOT20      --opt fps=30 --dims 2 --out data\trajectories\mot-mot20
python -m trajio export --source mot     --root data\downloads\mot\dataset\DanceTrack --opt fps=20 --dims 2 --out data\trajectories\mot-dancetrack

foreach ($site in "us-101", "i-80", "lankershim", "peachtree") {
  python -m trajio export --source ngsim --root data\downloads\ngsim `
    --opt location=$site --opt grouping=external --dims 2 --out data\trajectories\ngsim-$site
}
```

```sh
export PYTHONPATH=$PWD/src

python -m trajio export --source mopsi   --root data/downloads/mopsi --dims 2 --out data/trajectories/mopsi
python -m trajio export --source geolife --root "data/downloads/geolife/Geolife Trajectories 1.3" --dims 2 --out data/trajectories/geolife
python -m trajio export --source mot     --root data/downloads/mot/dataset/MOT17      --opt fps=30 --dims 2 --out data/trajectories/mot-mot17
python -m trajio export --source mot     --root data/downloads/mot/dataset/MOT20      --opt fps=30 --dims 2 --out data/trajectories/mot-mot20
python -m trajio export --source mot     --root data/downloads/mot/dataset/DanceTrack --opt fps=20 --dims 2 --out data/trajectories/mot-dancetrack

for site in us-101 i-80 lankershim peachtree; do
  python -m trajio export --source ngsim --root data/downloads/ngsim \
    --opt location=$site --opt grouping=external --dims 2 --out data/trajectories/ngsim-$site
done
```

Expect **47 916 trajectories and 47 102 992 points**. Each dataset prints its own counts, and
they should match the table in [datasets.md](datasets.md) exactly — if they do not, something
is wrong with the download rather than the reader.

`trajio export` writes into the directory without clearing it, so re-exporting after a format
change leaves the old files behind. Empty the directory first if the format has moved on.

## 5. Simplify

Every baseline, every trajectory, at `ceil(N / 2^m)` vertices for m = 1…6. Produces **786 642
documents, about 9 GB**.

```powershell
foreach ($ds in Get-ChildItem data\trajectories -Directory) {
  .\build\src\pipelines\simplify-sweep\ssk_simplify.exe --in $ds.FullName --out data\simplified-trajectories
}
```

```sh
for ds in data/trajectories/*/; do
  ./build/src/pipelines/simplify-sweep/ssk_simplify --in "$ds" --out data/simplified-trajectories
done
```

**Budget 2–3 core-hours**, dominated entirely by DOTS' threshold search. Eight of the nine
datasets take 4–36 minutes each. **GeoLife needs sharding** — a handful of its 60k-point
tracks dominate the whole dataset, because DOTS' cost grows faster than linearly:

```powershell
0..7 | ForEach-Object -Parallel {
  .\build\src\pipelines\simplify-sweep\ssk_simplify.exe --in data\trajectories\geolife `
    --out data\simplified-trajectories --shard $_ --shards 8
} -ThrottleLimit 8
```

```sh
seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/simplify-sweep/ssk_simplify \
  --in data/trajectories/geolife --out data/simplified-trajectories --shard {} --shards 8
```

`ForEach-Object -Parallel` needs PowerShell 7; on Windows PowerShell 5.1 use `Start-Job` in a
loop instead. Why the shards stride rather than take blocks, and where that still goes wrong,
is in [pipeline.md](pipeline.md).

One trajectory is knowingly excluded from DOTS: **`geolife-013552` does not finish**, four
attempts of 90+ CPU-minutes each. The sharded run above will sit on it forever, so stop that
worker when the other seven have finished and give the file its own pass with `--no-dots`,
which produces its Douglas–Peucker and SQUISH results in under a second. The driver takes a
directory, so put the one file in one — named `geolife`, since the dataset name comes from the
directory:

```powershell
mkdir -Force tmp\geolife | Out-Null
Copy-Item data\trajectories\geolife\geolife-013552.json tmp\geolife\
.\build\src\pipelines\simplify-sweep\ssk_simplify.exe --in tmp\geolife `
    --out data\simplified-trajectories --no-dots
Remove-Item -Recurse tmp
```

```sh
mkdir -p tmp/geolife
cp data/trajectories/geolife/geolife-013552.json tmp/geolife/
./build/src/pipelines/simplify-sweep/ssk_simplify --in tmp/geolife \
    --out data/simplified-trajectories --no-dots
rm -r tmp
```

### Bad documents cannot be produced

`DotsSimplifier`'s path decode steps backwards on some real inputs — the defect is in the
original implementation, and our port reproduces it faithfully — so a run can come back whose
points are **not a subsequence of the input**, which makes it not a simplification of anything.
It affected 1 034 documents of a corpus of 786 636 before it was caught, 96% of them at m1 and
90% of them in NGSIM.

**The driver now refuses to write such a run.** Every run is checked for strictly increasing
indices before it is serialised; a failing one is dropped, reported on stderr, and counted in
the summary:

```
ngsim-us-101: 1 trajectories, 1762 points -> 17 documents in 0.2 s
  1 runs dropped: indices not increasing, so not a subsequence
```

So a corpus produced by the current driver is clean by construction, and nothing needs cleaning
up afterwards. What this costs is coverage rather than correctness: DOTS ends up with slightly
fewer operating points than DPn and SQUISH, unevenly — 80.1% at m1 on `ngsim-us-101`, its worst
cell. That gap has to be stated in any results table, and
[comparison.md](comparison.md) says why.

Both incidents are written up in `archive/`.

### Verifying a corpus you did not just produce

The guard covers anything the current driver writes. To check a corpus produced by an older
driver, or after changing a port, `ssk_audit` reads what is actually on disk and matches each
result against its input on **position and timestamp together** — coordinates alone give false
failures wherever a trace repeats a fix:

```sh
for ds in data/trajectories/*/; do
  ./build/src/pipelines/audit/ssk_audit --in "$ds" --results data/simplified-trajectories
done
```

```powershell
foreach ($ds in Get-ChildItem data\trajectories -Directory) {
  .\build\src\pipelines\audit\ssk_audit.exe --in $ds.FullName --results data\simplified-trajectories
}
```

It exits non-zero when it finds anything, so it drops straight into a script. `--paths` prints
only the offending paths, which makes removal a pipe — the tool never deletes anything itself:

```sh
./build/src/pipelines/audit/ssk_audit --in data/trajectories/ngsim-us-101 \
    --results data/simplified-trajectories --paths | xargs rm
```

```powershell
.\build\src\pipelines\audit\ssk_audit.exe --in data\trajectories\ngsim-us-101 `
    --results data\simplified-trajectories --paths | Remove-Item
```

It is not fast — it reads every document — so budget a couple of hours for the whole corpus,
and run it per dataset in parallel as with the sweep.

## 6. Measure the error

The eight smaller datasets first — about 25 core-minutes between them:

```powershell
foreach ($ds in Get-ChildItem data\trajectories -Directory -Exclude geolife) {
  .\build\src\pipelines\frechet-distance\ssk_frechet.exe --in $ds.FullName `
    --simplified data\simplified-trajectories --out data\frechet-distances
}
```

```sh
for ds in data/trajectories/*/; do
  [ "$(basename "$ds")" = geolife ] && continue
  ./build/src/pipelines/frechet-distance/ssk_frechet --in "$ds" \
    --simplified data/simplified-trajectories --out data/frechet-distances
done
```

Then GeoLife, which is **77% of all the work** and has to be split within itself — running it
as one process would take five hours while seven cores idle:

```powershell
0..7 | ForEach-Object -Parallel {
  .\build\src\pipelines\frechet-distance\ssk_frechet.exe --in data\trajectories\geolife `
    --simplified data\simplified-trajectories --out data\frechet-distances `
    --shard $_ --shards 8
} -ThrottleLimit 8
```

```sh
seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/frechet-distance/ssk_frechet \
  --in data/trajectories/geolife --simplified data/simplified-trajectories \
  --out data/frechet-distances --shard {} --shards 8
```

**Budget about 6.6 core-hours of compute, but roughly three hours of wall clock**: the job
reads every result document and its input, and that I/O dominates, exactly as it did for the
corpus audit. Cost is `n × m` cells per pair, so it is quadratic in trajectory length at a
fixed rate and **GeoLife alone is 77% of the total** — sharding by dataset achieves nothing on
its own, the split has to be within GeoLife. `--shard` works here exactly as above.

A whole-corpus pass runs for hours, so if one dies, add **`--skip-existing`** when restarting:
a measurement is a pure function of the two documents it reads, so one already written need
not be computed again.

Unlike the sweep, this cost is known in advance: cells are `n × m` exactly at a measured
34.45 ns each, and the worst single trajectory in the corpus is 15 minutes on one core. Nothing
can stall — there is no search and no data-dependent branching. The per-dataset breakdown is in
[frechet-distance.md](frechet-distance.md).

The result is one document per simplification at
`data/frechet-distances/<algorithm>/<dataset>/m<rate>/<name>.json`, carrying the distance
beside the `params` and `stats` that produced it. That document is self-contained: a results
table does not need to open the simplification or the input to know what it is looking at.

## 7. Read the results

Everything needed for a table is in the `frechet-distances` tree. A rate-by-algorithm summary
is a short script over it, and this one is the same on both platforms:

```python
import json, pathlib, collections
rows = collections.defaultdict(list)
for f in pathlib.Path("data/frechet-distances").rglob("*.json"):
    d = json.loads(f.read_text())
    rows[(d["algorithm"], d["stats"]["m"])].append(d["frechet"]["distance"])
for key in sorted(rows):
    v = rows[key]
    print(f"{key[0]:20} m{key[1]}  n={len(v):>7}  mean {sum(v)/len(v):.3e}")
```

Two things to know before drawing conclusions from it, both in
[comparison.md](comparison.md) and [pipeline.md](pipeline.md): the baselines are being scored
against an objective **none of them optimises**, and **error is not monotone in the vertex
budget** — a kept vertex forces the matching through it, so a finer simplification can score
worse than a coarser one. Means behave; individual series do not.

---

## Visualising

`src/viz/plot.py` renders any document in the project's format — an input trajectory, or a
simplification — with a statistics panel beside it. It uses matplotlib's `Agg` backend, so it
writes files and needs no display, which also means it works over SSH and in a container.

```powershell
python src\viz\plot.py data\trajectories\mopsi\mopsi-000003.json -o data\renders\input.png
python src\viz\plot.py data\simplified-trajectories\dots\mopsi\m3\mopsi-000003.json -o data\renders\simplified.png
```

```sh
python src/viz/plot.py data/trajectories/mopsi/mopsi-000003.json -o data/renders/input.png
python src/viz/plot.py data/simplified-trajectories/dots/mopsi/m3/mopsi-000003.json \
                       -o data/renders/simplified.png
```

Because a simplification is a trajectory document in its own right, it draws directly with no
conversion step. 3D documents render as a 3D line plot plus the `xy`, `xz` and `yz`
projections, since one 3D view hides most of what a simplification did.

### Input and simplification on the same axes

The result documents deliberately do **not** carry a copy of the input — that would repeat the
whole trajectory in each of the ~14 documents per trajectory. `plot.py` draws the input
underneath when the document has an `input_points` array, so build a throwaway document that
has one. This also fills in the statistics panel's Fréchet line. Python, so identical on both
platforms:

```python
import json, pathlib

name, algo, rate, ds = "mopsi-000003.json", "dots", "m3", "mopsi"

src = json.loads(pathlib.Path(f"data/trajectories/{ds}/{name}").read_text())
res = json.loads(pathlib.Path(f"data/simplified-trajectories/{algo}/{ds}/{rate}/{name}").read_text())
fr = json.loads(pathlib.Path(f"data/frechet-distances/{algo}/{ds}/{rate}/{name}").read_text())

res["input_points"] = src["points"]
res["frechet"] = {"distance": fr["frechet"]["distance"],
                  "bound": fr["frechet"]["distance"],
                  "within_bound": True,
                  "computer": fr["frechet"]["computer"]}
pathlib.Path("data/renders").mkdir(parents=True, exist_ok=True)
pathlib.Path("data/renders/overlay.json").write_text(json.dumps(res))
```

```powershell
python src\viz\plot.py data\renders\overlay.json -o data\renders\overlay.png
```

```sh
python src/viz/plot.py data/renders/overlay.json -o data/renders/overlay.png
```

The input appears in grey behind the simplification in red, with its kept vertices marked, and
the panel reports the algorithm, the compression achieved and the measured Fréchet distance.
That picture is the quickest way to see whether a number is believable — a distance that looks
wrong usually looks wrong on the plot too.

### Synthetic curves

For working on an algorithm without touching a dataset:

```powershell
python src\viz\generate.py spiral2d      -n 600 -o data\synthetic\spiral2d.json
python src\viz\generate.py shrink-spiral -n 75  -o data\synthetic\shrink_spiral3d.json
```

```sh
python src/viz/generate.py spiral2d      -n 600 -o data/synthetic/spiral2d.json
python src/viz/generate.py shrink-spiral -n 75  -o data/synthetic/shrink_spiral3d.json
```

`spiral2d`, `zigzag2d`, `walk2d`, `u-shape`, `o-shape` in 2D and `shrink-spiral` in 3D.
`u-shape` and `o-shape` are the two that stress a simplifier — near-parallel legs whose
tolerance balls overlap, and a closed loop with both endpoints pinned. They take the same
pipelines as a real dataset: point `--in` at `data/synthetic`.

---

## Portability

**No source file is platform-specific.** The C++ uses `std::filesystem` throughout with no
`_WIN32`, no `windows.h` and no `_MSC_VER`; CMake has no `WIN32` branch and never spells
`.exe`; every Python file is opened with an explicit `encoding="utf-8"` and `newline="\n"`, so
the bytes are identical on both platforms rather than depending on the system codepage; and
`.gitattributes` pins `eol=lf`. Nothing in `src/` shells out, so no shell is required at run
time.

One detail that keeps the *data* portable too: the `source` and `simplified` fields in result
documents are built by string concatenation with `/`, not by `std::filesystem::path`. Had they
been paths, Windows would have written backslashes into the JSON and the documents would not
move between machines.

What differs is only ever the shell:

| | Windows (PowerShell) | Linux / macOS |
|---|---|---|
| activate the venv | `.\.venv\Scripts\Activate.ps1` | `source .venv/bin/activate` |
| set `PYTHONPATH` | `$env:PYTHONPATH = "$PWD\src"` | `export PYTHONPATH=$PWD/src` |
| configure CMake | add `-G "MinGW Makefiles"` | no generator flag needed |
| run a binary | `.\build\...\ssk_simplify.exe` | `./build/.../ssk_simplify` |
| download | `curl.exe` (bare `curl` is an alias) | `curl` |
| loop over datasets | `foreach ($ds in Get-ChildItem ... -Directory)` | `for ds in data/trajectories/*/; do … done` |
| run N shards at once | `ForEach-Object -Parallel` (PS 7) or `Start-Job` (5.1) | `xargs -P N` |
| line continuation | backtick `` ` `` | backslash `\` |

Two Windows-only annoyances worth knowing, neither of which affects the code: MinGW-built
binaries do not run from Git Bash because its `PATH` lacks the runtime DLLs — use PowerShell —
and a filename cannot contain `:`, which is why the incident reports in `archive/` use a dash.

The project has been built and run end to end on Windows with MinGW g++ 15.1. Nothing in it
should stop it building on Linux with GCC 11+ or Clang 14+, but that has not been tried, so
treat it as "no known obstacles" rather than "verified".
