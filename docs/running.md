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

56 tests, a few seconds. Run them before anything long — they cover the simplifiers against
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

About **1.96 GB**, and the NGSIM download alone takes ~40 minutes. Full commands, checksums
and what each dataset contains are in [datasets.md](datasets.md); the short version:

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

MOT comes from the Trajectory Simplify Benchmark's Google Drive copies — see
[datasets.md](datasets.md), which explains why `motchallenge.net` is not used.

Unpack the two archives in place. GeoLife and Mopsi are read from their extracted trees;
NGSIM is one 1.5 GB CSV read as-is.

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

One trajectory is knowingly excluded from DOTS: `geolife-013552` does not finish. Run it with
`--no-dots` so it still gets Douglas–Peucker and SQUISH results.

## 6. Measure the error

```powershell
foreach ($ds in Get-ChildItem data\trajectories -Directory) {
  .\build\src\pipelines\frechet-distance\ssk_frechet.exe --in $ds.FullName `
    --simplified data\simplified-trajectories --out data\frechet-distances
}
```

```sh
for ds in data/trajectories/*/; do
  ./build/src/pipelines/frechet-distance/ssk_frechet --in "$ds" \
    --simplified data/simplified-trajectories --out data/frechet-distances
done
```

**Budget about 6.6 core-hours**, under an hour on eight. Cost is `n × m` cells per pair, so it
is quadratic in trajectory length at a fixed rate and GeoLife dominates again — `--shard`
works here exactly as above.

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
