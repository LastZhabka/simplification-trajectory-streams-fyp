# From zero to Fréchet distances

Every step from a fresh clone to a measured error for every simplification, to ours run at
each of those errors (step 7), to both sides held to the same fixed errors (step 8), to
both again with ours run at the error itself rather than inside it (steps 9 and 10), and to
ours searched to each baseline's size (step 11).
**Each step gives the PowerShell and the POSIX shell form**, since the
project is developed on Windows but nothing in it is Windows-only — see *Portability* at the
end for what actually differs.

Run everything from the repository root. Each step says roughly what it costs, because three of
them are measured in hours.

---

## 0. Starting over

No stage clears the directory it writes into, so a rerun onto an old tree mixes the two. Remove
what you are about to regenerate. `data/error-matched` and `data/target-matched` also hold the
reports step 12 wrote and the shard logs; move those out first if you want them.

```powershell
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue `
  data\simplified-trajectories, data\frechet-distances, data\error-matched, `
  data\target-simplified, data\target-matched, data\error-matched-unscaled, `
  data\target-matched-unscaled, data\size-matched, data\ratio-buckets.json, data\plots
```

```sh
rm -rf data/simplified-trajectories data/frechet-distances data/error-matched \
  data/target-simplified data/target-matched data/error-matched-unscaled \
  data/target-matched-unscaled data/size-matched data/ratio-buckets.json data/plots
```

**Let a delete finish before anything writes.** Step 7's tree alone is 4.5 million files and
takes over 10 minutes to remove; a run started meanwhile writes into a directory the delete has
not reached yet, and loses what it wrote.

Steps 3 and 4 only need repeating if `data/downloads` or `data/trajectories` is gone or the
export format changed; to redo step 4, empty each `data/trajectories/<dataset>` first. A
`mopsi-sample` directory left there from earlier work is a subset of `mopsi`: delete it, or
every loop below processes those 40 trajectories twice.

## 1. Toolchain

A C++20 compiler and CMake 3.20+. Everything under `src/algo` and `src/pipelines` is standard
library only, so there is nothing to fetch and the build works offline. Python 3.10+ is needed
for ingestion and the statistics, and it too is standard library only: nothing to install.

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

86 tests, a few seconds. Run them before anything long — they cover the simplifiers against
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
the repository, so make the four first. MOT is fetched with `gdown`, which nothing in the
project imports — install it for this step alone:

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

Every baseline — DPn, SQUISH, SQUISH-E and DOTS — on every trajectory, at `ceil(N / 2^m)`
vertices for m = 2…6 (`--rates A..B` picks others): **857 461 documents**, predicted per
dataset in [datasets.md](datasets.md). At m = 1…6, the three without SQUISH-E produced
786 642 documents, about 9 GB. The eight smaller datasets first:

```powershell
foreach ($ds in Get-ChildItem data\trajectories -Directory -Exclude geolife) {
  .\build\src\pipelines\simplify-sweep\ssk_simplify.exe --in $ds.FullName --out data\simplified-trajectories
}
```

```sh
for ds in data/trajectories/*/; do
  [ "$(basename "$ds")" = geolife ] && continue
  ./build/src/pipelines/simplify-sweep/ssk_simplify --in "$ds" --out data/simplified-trajectories
done
```

**Budget 2–3 core-hours**, dominated entirely by DOTS' threshold search. Eight of the nine
datasets take 4–36 minutes each. **GeoLife needs sharding** — a handful of its 60k-point
tracks dominate the whole dataset, because DOTS' cost grows faster than linearly.

One trajectory is knowingly excluded from DOTS: **`geolife-013552` does not finish**, four
attempts of 90+ CPU-minutes each. A worker that reaches it never moves on, so every trajectory
behind it in that shard would go unsimplified. Take it out of the directory before the sharded
run, give it its own pass with `--no-dots` — its other three baselines take seconds — and
put it back, because step 6 measures those three. The driver takes a directory,
and the dataset name comes from the directory, so the one file goes in one named `geolife`:

```powershell
mkdir -Force tmp\geolife | Out-Null
Move-Item data\trajectories\geolife\geolife-013552.json tmp\geolife\

$procs = 0..7 | ForEach-Object {
  Start-Process .\build\src\pipelines\simplify-sweep\ssk_simplify.exe -NoNewWindow -PassThru `
    -ArgumentList "--in data\trajectories\geolife --out data\simplified-trajectories --shard $_ --shards 8"
}
$procs | Wait-Process

.\build\src\pipelines\simplify-sweep\ssk_simplify.exe --in tmp\geolife `
    --out data\simplified-trajectories --no-dots
Move-Item tmp\geolife\geolife-013552.json data\trajectories\geolife\
Remove-Item -Recurse tmp
```

```sh
mkdir -p tmp/geolife
mv data/trajectories/geolife/geolife-013552.json tmp/geolife/

seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/simplify-sweep/ssk_simplify \
  --in data/trajectories/geolife --out data/simplified-trajectories --shard {} --shards 8

./build/src/pipelines/simplify-sweep/ssk_simplify --in tmp/geolife \
    --out data/simplified-trajectories --no-dots
mv tmp/geolife/geolife-013552.json data/trajectories/geolife/
rm -r tmp
```

`Start-Process` works on Windows PowerShell 5.1 and on 7. Why the shards stride rather than
take blocks, and where that still goes wrong, is in [pipeline.md](pipeline.md).

The sweep has no `--skip-existing`: if a dataset's run dies, run that dataset again.

### Bad documents cannot be produced

Upstream's `DotsSimplifier` can join a point to a frontier vertex that comes after it, so its
decoded path steps backwards and the points are **not a subsequence of the input**. That
affected 1 034 documents of a corpus of 786 636 before it was caught, 96% of them at m1 and 90%
of them in NGSIM. **Our port no longer does it** — the one place it departs from upstream's
algorithm, described in [simplification.md](simplification.md).

The driver's guard stays regardless. Every run is checked for strictly increasing indices
before it is serialised; a failing one is dropped, reported on stderr, and counted in the
summary as `N runs dropped: indices not increasing, so not a subsequence`. With the fix that
line should not appear, and if it does, it is a bug.

DOTS results from before the fix are not comparable with new ones wherever the old code made
a backwards link, so regenerate the DOTS tree rather than mixing the two. Both incidents are
written up in `archive/`.

## 6. Measure the error

The eight smaller datasets first — about 90 core-minutes of compute between them, over half of
it Mopsi, so run them concurrently rather than in the loop below if you have the cores:

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
$procs = 0..7 | ForEach-Object {
  Start-Process .\build\src\pipelines\frechet-distance\ssk_frechet.exe -NoNewWindow -PassThru `
    -ArgumentList "--in data\trajectories\geolife --simplified data\simplified-trajectories --out data\frechet-distances --shard $_ --shards 8"
}
$procs | Wait-Process
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
not be computed again. Every pipeline writes a document as `<name>.json.part` and renames it
once it is whole, so a killed run never leaves a truncated document for `--skip-existing` to
take as done. A leftover `.part` is harmless and can be deleted.

Unlike the sweep, this cost is known in advance: cells are `n × m` exactly at a measured
34.45 ns each, and the worst single trajectory in the corpus is 15 minutes on one core. Nothing
can stall — there is no search and no data-dependent branching. The per-dataset breakdown is in
[frechet-distance.md](frechet-distance.md).

The result is one document per simplification at
`data/frechet-distances/<algorithm>/<dataset>/m<rate>/<name>.json`, carrying the distance
beside the `params` and `stats` that produced it. That document is self-contained: a results
table does not need to open the simplification or the input to know what it is looking at.

## 7. Run ours at each baseline's error

`ssk_match` runs ours once per baseline document from step 6, at
`δ = d*/(1+ε)` where `d*` is that baseline's measured error. Ours is then guaranteed no worse
than the baseline, and the question is how many vertices it needs. It measures ours and checks
the bound exactly in the same pass. What it computes and why is in [algorithm.md](algorithm.md);
this is how to run it.

It needs the build and step 6's `data/frechet-distances/`; a trajectory with no baseline
document is skipped. **`--eps` is required**, and each value is its own run under its own root.

**`--ours` picks the implementation**: `greedy-stabbing-shared`, the default, O(ε⁻¹) per
vertex ([convex-shadow-32.md](convex-shadow-32.md) §11); `greedy-stabbing-naive`, the same
without the sharing, O(ε^(−3/2)) per vertex (§10 there); or `greedy-stabbing-naive-paper` and
`greedy-stabbing-naive-paper-boundary`, Algorithm 1 with the paper's `conv(G_v)` and grid
starts; or `greedy-stabbing-no-decomp`, `greedy-stabbing-no-decomp-paper` and
`greedy-stabbing-no-decomp-paper-boundary`, the same three stepped as the paper does, a polygon
per start (§12).
Each is also its own run, under
`data/error-matched/<ours>/`. The commands below set it once at the top as `$ours`; to run
another, change that line and run them again.

One dataset at all three ε, 8 processes at a time:

```powershell
$ours = "greedy-stabbing-shared"
$ds = "mot-dancetrack"
New-Item -ItemType Directory -Force data\error-matched | Out-Null
foreach ($eps in "0.5", "0.3", "0.1") {
  $procs = 0..7 | ForEach-Object {
    Start-Process .\build\src\pipelines\error-match\ssk_match.exe -NoNewWindow -PassThru `
      -ArgumentList "--in data\trajectories\$ds --ours $ours --eps $eps --shard $_ --shards 8 --skip-existing" `
      -RedirectStandardOutput "data\error-matched\$ours-$ds-eps$eps-shard$_.log"
  }
  $procs | Wait-Process
}
Get-Content data\error-matched\$ours-$ds-eps*-shard*.log
```

```sh
ours=greedy-stabbing-shared
ds=mot-dancetrack
mkdir -p data/error-matched
for eps in 0.5 0.3 0.1; do
  seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/error-match/ssk_match \
    --in data/trajectories/$ds --ours $ours --eps $eps --shard {} --shards 8 --skip-existing
done 2>&1 | tee -a data/error-matched/$ours-$ds.log
```

`Start-Process` works on PowerShell 5.1 and 7. Prefer it to `Start-Job` here: each job is a whole
PowerShell host of about 70 MB, while `ssk_match` itself peaks at 7 MB.

For the corpus, wrap that in a loop over the datasets, one ε at a time so the cheap ones
finish first:

```powershell
$ours = "greedy-stabbing-shared"
New-Item -ItemType Directory -Force data\error-matched | Out-Null
foreach ($eps in "0.5", "0.3", "0.1") {
  foreach ($ds in (Get-ChildItem data\trajectories -Directory).Name) {
    $procs = 0..7 | ForEach-Object {
      Start-Process .\build\src\pipelines\error-match\ssk_match.exe -NoNewWindow -PassThru `
        -ArgumentList "--in data\trajectories\$ds --ours $ours --eps $eps --shard $_ --shards 8 --skip-existing" `
        -RedirectStandardOutput "data\error-matched\$ours-$ds-eps$eps-shard$_.log"
    }
    $procs | Wait-Process
  }
}
```

```sh
ours=greedy-stabbing-shared
mkdir -p data/error-matched
for eps in 0.5 0.3 0.1; do
  for ds in data/trajectories/*/; do
    ds=$(basename "$ds")
    seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/error-match/ssk_match \
      --in data/trajectories/$ds --ours $ours --eps $eps --shard {} --shards 8 --skip-existing \
      2>&1 | tee -a data/error-matched/$ours-$ds.log
  done
done
```

Each implementation is another full step 7. If time is short, run `greedy-stabbing-shared` on
the whole corpus, and `greedy-stabbing-naive` on one dataset at one ε: it gives
`greedy-stabbing-shared`'s output to the bit on every run checked so far (§11 of
[convex-shadow-32.md](convex-shadow-32.md)), and that dataset is enough to confirm it.

It is resumable as a whole: run it again and finished trajectories are skipped.

Each shard ends with a line like this one, from a run of the former `convex-shadow`:

```
convex-shadow mot-dancetrack eps=0.5: 87 trajectories -> 1553 matches over 1337209 points in 1156.1 s
  out: data/error-matched/convex-shadow/eps0.5/{<algorithm>,summary}/mot-dancetrack/
```

On a resume it adds `N already done, skipped`. A run that breaks the bound adds
`K over (1+eps)*delta -- a broken guarantee, look at these`. If it does, it is a bug. It
appeared once, on 2026-11-24, at ε = 0.5: 11 of 857 456 matches, in two Mopsi and three MOT20
trajectories, up to twice the bound. `geometry::closest` returned the query point as inside a
shadow that had shrunk to a sliver far from it; fixed the same day, with a regression test in
`tests/test_greedy_stabbing.cpp`, and again on 2026-11-26 for such a sliver far from the first
vertex, where the first fix did not reach ([results.md](results.md)). A resumed run rewrites the logs of the shards it skips, so
check the result documents (`within_bound` in each summary) as well as the logs.

**The precision floor.** A baseline whose error d* equals its measurement tolerance has error 0
within measurement: on Mopsi and NGSIM a budget can drop exact repeats, and on NGSIM runs of
fixes collinear to the bit, at no error at all. Ours is then asked for δ = d*/(1+ε) (step 7)
or d* (step 9): a few ulps of the coordinates. `GreedyStabber::run` builds each vertex's region
for δ less what its arithmetic can err there, from what has arrived, as a stream must. Each run
computes relative to its own first vertex o, so at a vertex T the margin is
`Δ = DBL_EPSILON · (‖o‖∞ + 2δ + 64 ‖T − o‖∞)`: writing p and q back rounds by ulps of the
absolute coordinates, and every step errs by ulps of how far T is from o. 64 is a safety factor.
The regions keep their normals, so only the offsets scale, and a run's starts are built for
δ − Δ at o. On NGSIM near the floor Δ is 1.4·10⁻⁹ ft. Where Δ ≥ δ no new point can be written
within δ of T, and `run` outputs T itself, dropping exact repeats, which is exact. Until
2026-11-27 the margin was the whole curve's, from its largest coordinate and its extent, which a
stream cannot know.

Fixed twice on 2026-11-24:

- Under the margin `run` returned the input whole, repeats and all: 55 pairs in 13
  trajectories over Theorem 3's size bound. Now it drops the repeats.
- The margin was 16 · DBL_EPSILON · max|coord|, all of it charged to the absolute coordinates:
  2.3·10⁻⁸ ft on NGSIM, 25 ulps. Step 9 then found ngsim-002153, whose DPn drops runs
  collinear to the bit at error 0; ours at δ = 3.4·10⁻⁹ ft fell back and kept 170 vertices
  against 2 · 61 − 2. With the margin split, ours stabs those runs (regression test in
  `tests/test_greedy_stabbing.cpp`).

What is left under the floor: 98 pairs per ε in step 7 and in step 9 (168–208 and 164 under the
old margin), all in 7 trajectories that barely move, with extents of 0.03–0.9 ft on NGSIM
Lankershim and 7·10⁻⁶ degrees on GeoLife. Every baseline's error there is 0 within
measurement; ours keeps their 3 to 7 distinct points and is within both bounds. Steps 8 and
10's targets never come near the floor.

**Cost, measured on `mot-dancetrack`** (692 trajectories, 12 335 baseline documents; 8
processes on 8 cores):

| ε | wall clock | core-hours |
|---|---:|---:|
| 0.5 | 19.6 min | 2.5 |
| 0.3 | 30.4 min | 3.9 |
| 0.1 | 143 min | ≈ 19.5 |

That is about **1.3× the estimate in [algorithm.md](algorithm.md) §7**. Scale its corpus figures
(189 / 300 / 1 356 core-hours) by about 1.3 before planning a corpus run, and do ε = 0.5 and 0.3
first.

On disk, each ε of `mot-dancetrack` is 13 027 files and 55 MB: one result per baseline
document plus one summary per trajectory.

**Resuming.** The summary is written last, so a trajectory that has one is finished. Rerunning
the same command (`--skip-existing` is already in it) picks up where it stopped, losing at most
the trajectory each shard was on. The shard count may change between runs.

**Check before reading the results:**

- **Summaries:** each implementation at each ε has one per trajectory in the dataset.
- **Logs:** no `over (1+eps)*delta` line.
- **Coverage:** in step 12's coverage table, each `<ours> eps=E` cell is the sum of the four
  baselines in its row.

```powershell
(Get-ChildItem data\error-matched\greedy-stabbing-shared\eps0.5\summary\mot-dancetrack).Count
Select-String -Path data\error-matched\*-shard*.log -Pattern "over \(1\+eps\)"
```

```sh
ls data/error-matched/greedy-stabbing-shared/eps0.5/summary/mot-dancetrack | wc -l
grep "over (1+eps)" data/error-matched/*.log
```

The first prints the trajectory count (692 for `mot-dancetrack`); the second prints nothing.

## 8. Fix the error: each baseline's best compression, and ours, within the same target

Step 7 fixes a compression rate and asks how small ours gets at the baseline's error. This is
the other protocol of [comparison.md](comparison.md), *fix error, compare compression*: fix
an error target T, give every baseline the fewest vertices it can find within T, and run ours
within the same T.

**The targets** are T = f · the trajectory's bounding-box diagonal, for
f ∈ {0.0003, 0.001, 0.003, 0.01, 0.03}. The diagonal makes a level mean the same in degrees,
feet and pixels. The five fractions are half a decade apart and span the medians of step 6's
d*/diagonal at m2…m6 (8·10⁻⁴ to 2.5·10⁻², all four baselines pooled), so the levels play the
part the rates play in step 7. A trajectory with no extent, 1.3% of a sample, has no targets.

**A baseline at a target** is the fewest vertices its knob was found to give within T, each
step decided exactly by `Frechet::within`. DPn and SQUISH bisect on their budget, SQUISH-E on
a budget k run as `cmp_ratio = ⌈(N − 1)/k⌉`, DOTS on log(`lssd_threshold`). The budget's
midpoint is geometric, because a decision costs N·k and the targets sit at large compression.
The five targets are solved smallest first, and each answer brackets the next, as it is within
every larger target. **None of these is guaranteed monotone** — a kept vertex forces the
matching through it, so fewer vertices can be closer — so the bisection stops at a boundary:
the size is the fewest found, not proven the fewest. Where no setting tried is within T, the
baseline gets the input itself (`"input": true`), compression 1. The tight levels do this on
noisy traces: SQUISH-E cannot drop even one point of a fast wiggle without breaking f = 0.0003.

**Ours at a target** runs at `delta = T/(1 + ε)`, so it is within T by its guarantee, and is
checked exactly. It does not depend on the baseline, so it runs once per level, not once per
pair.

Two passes, like steps 5–7: `ssk_target_search` is the expensive one and runs once;
`ssk_target_match` is ours per ε, cheap to repeat after a change to ours. The search needs
only step 4's trajectories; the match needs the search.

```powershell
$procs = foreach ($ds in (Get-ChildItem data\trajectories -Directory).Name) {
  0..7 | ForEach-Object {
    Start-Process .\build\src\pipelines\target-match\ssk_target_search.exe -NoNewWindow -PassThru `
      -ArgumentList "--in data\trajectories\$ds --shard $_ --shards 8 --skip-existing"
  }
}
$procs | Wait-Process
```

```sh
for ds in data/trajectories/*/; do
  seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/target-match/ssk_target_search \
    --in "$ds" --shard {} --shards 8 --skip-existing
done
```

**Trajectories over 30 000 points are searched without DOTS**: taken out first, given their
own `--no-dots` pass, and put back. That is 31 trajectories, 0.06% of the corpus, all GeoLife
and Mopsi, `geolife-013552` among them. DOTS' search costs about 3 minutes at 30 669 points,
over 4 at 45 993, and more than 7 at geolife-001523's 92 645, where one DOTS run that
collapses the curve to 6 points took 900 s. Below 30 000 it stays: 6 s at 9 981 points, 31 s at
19 960. Their other three baselines take at most a couple of minutes.

Length is a proxy, not the cause, and it missed one: **mopsi-003388** (27 774 points) ran DOTS'
search for over 7 minutes without finishing, where its other three baselines took 10.5 s. It
was given the same `--no-dots` pass, so 32 trajectories have no DOTS result in step 8. If
another shard stalls on DOTS, find its trajectory (the shard's next file without a summary),
kill the shard, and do the same: the rest of the shard resumes with `--skip-existing`.

```powershell
foreach ($ds in "geolife", "mopsi") {
  New-Item -ItemType Directory -Force tmp\$ds | Out-Null
  (Get-Content data\trajectories\$ds\index.json -Raw | ConvertFrom-Json).trajectories |
    Where-Object { $_.points -gt 30000 } |
    ForEach-Object { Move-Item data\trajectories\$ds\$($_.file) tmp\$ds\ }
}
# ... the search above ...
foreach ($ds in "geolife", "mopsi") {
  .\build\src\pipelines\target-match\ssk_target_search.exe --in tmp\$ds --no-dots --skip-existing
  Move-Item tmp\$ds\*.json data\trajectories\$ds\
}
Remove-Item -Recurse tmp
```

```sh
for ds in geolife mopsi; do
  mkdir -p tmp/$ds
  python3 -c "import json,sys; [print(t['file']) for t in json.load(open(sys.argv[1]))['trajectories'] if t['points'] > 30000]" \
    data/trajectories/$ds/index.json | while read f; do mv data/trajectories/$ds/$f tmp/$ds/; done
done
# ... the search above ...
for ds in geolife mopsi; do
  ./build/src/pipelines/target-match/ssk_target_search --in tmp/$ds --no-dots --skip-existing
  mv tmp/$ds/*.json data/trajectories/$ds/
done
rm -r tmp
```

Step 12's tables then count 32 fewer pairs for DOTS than for the other three at each level.

Then ours, one ε at a time:

```powershell
$ours = "greedy-stabbing-shared"
foreach ($eps in "0.5", "0.3", "0.2", "0.1", "0.05") {
  $procs = foreach ($ds in (Get-ChildItem data\trajectories -Directory).Name) {
    0..7 | ForEach-Object {
      Start-Process .\build\src\pipelines\target-match\ssk_target_match.exe -NoNewWindow -PassThru `
        -ArgumentList "--in data\trajectories\$ds --ours $ours --eps $eps --shard $_ --shards 8 --skip-existing" `
        -RedirectStandardOutput "data\target-matched\$ours-$ds-eps$eps-shard$_.log"
    }
  }
  $procs | Wait-Process
}
```

```sh
ours=greedy-stabbing-shared
mkdir -p data/target-matched
for eps in 0.5 0.3 0.2 0.1 0.05; do
  for ds in data/trajectories/*/; do
    seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/target-match/ssk_target_match \
      --in "$ds" --ours $ours --eps $eps --shard {} --shards 8 --skip-existing
  done
done 2>&1 | tee -a data/target-matched/$ours.log
```

The PowerShell form needs `data\target-matched` to exist first (`New-Item -ItemType Directory
-Force data\target-matched`).

**What it writes:**

| path | what |
|---|---|
| `data/target-simplified/<algorithm>/<dataset>/f<fraction>/<name>.json` | a baseline's best within the target: its points, `params` (the knob), `stats` (`output_size`, `algorithm_runs`, `input`), `target` and its measured `frechet` |
| `data/target-simplified/summary/<dataset>/<name>.json` | the targets and every baseline's result, written last |
| `data/target-matched/<ours>/eps<E>/results/<dataset>/f<fraction>/<name>.json` | ours at a target |
| `data/target-matched/<ours>/eps<E>/summary/<dataset>/<name>.json` | ours beside each baseline at each level, the shape of step 7's summary with `level` where it has `m` |

Both passes resume with `--skip-existing`, as the summary is written last.

**Check before reading the results:** each pass has one summary per trajectory; no log has a
`results dropped` line (a baseline that did not keep a subsequence) or an `over the target`
line (ours broke its guarantee); and step 12's tables count `n` the same for each baseline at a
level, short trajectories apart.

## 9 and 10. The same, with ours at δ = the error

Steps 7 and 8 run ours at δ = error/(1 + ε), so that it is within the baseline's error. Theorem
3's size bound is then against κ(error/(1 + ε)), the optimum at a stricter error than the
baseline is held to, not against κ(error). Ours is bicriteria — run at δ, within (1 + ε)δ and
at most 2κ(δ) − 2 vertices — so no setting puts both at the baseline's error. Steps 9 and 10
take the other side of that choice: ours at **δ = d\*** (step 9) and **δ = T** (step 10), its
size bound against κ of the baseline's own error, and its error allowed up to (1 + ε) times
the baseline's. Wherever the baseline's error is within δ, which is every pair here, its own
size bounds κ, so Theorem 3 reads directly: **ours ≤ 2 · the baseline's size − 2**. Ours uses
its whole bound (a median 1.00 of it in step 8), so it really does take the extra error: read
these steps for what the size guarantee does, and steps 7–8 for who is smaller at the same
error.

Both are steps 7 and 8's drivers with `--unscaled-delta`, which also moves their output to
`data/error-matched-unscaled` and `data/target-matched-unscaled`. They read what steps 6 and
8 read, and nothing of steps 7 or 8's output.

```powershell
$ours = "greedy-stabbing-shared"
New-Item -ItemType Directory -Force data\error-matched-unscaled, data\target-matched-unscaled | Out-Null
foreach ($eps in "0.5", "0.3", "0.2", "0.1", "0.05") {
  $procs = foreach ($ds in (Get-ChildItem data\trajectories -Directory).Name) {
    0..7 | ForEach-Object {
      Start-Process .\build\src\pipelines\target-match\ssk_target_match.exe -NoNewWindow -PassThru `
        -ArgumentList "--in data\trajectories\$ds --ours $ours --eps $eps --unscaled-delta --shard $_ --shards 8 --skip-existing" `
        -RedirectStandardOutput "data\target-matched-unscaled\$ours-$ds-eps$eps-shard$_.log"
    }
  }
  $procs | Wait-Process
}
foreach ($eps in "0.5", "0.3", "0.2", "0.1", "0.05") {
  $procs = foreach ($ds in (Get-ChildItem data\trajectories -Directory).Name) {
    0..7 | ForEach-Object {
      Start-Process .\build\src\pipelines\error-match\ssk_match.exe -NoNewWindow -PassThru `
        -ArgumentList "--in data\trajectories\$ds --ours $ours --eps $eps --unscaled-delta --shard $_ --shards 8 --skip-existing" `
        -RedirectStandardOutput "data\error-matched-unscaled\$ours-$ds-eps$eps-shard$_.log"
    }
  }
  $procs | Wait-Process
}
```

```sh
ours=greedy-stabbing-shared
mkdir -p data/error-matched-unscaled data/target-matched-unscaled
for eps in 0.5 0.3 0.2 0.1 0.05; do
  for ds in data/trajectories/*/; do
    seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/target-match/ssk_target_match \
      --in "$ds" --ours $ours --eps $eps --unscaled-delta --shard {} --shards 8 --skip-existing
  done
done 2>&1 | tee -a data/target-matched-unscaled/$ours.log
for eps in 0.5 0.3 0.2 0.1 0.05; do
  for ds in data/trajectories/*/; do
    seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/error-match/ssk_match \
      --in "$ds" --ours $ours --eps $eps --unscaled-delta --shard {} --shards 8 --skip-existing
  done
done 2>&1 | tee -a data/error-matched-unscaled/$ours.log
```

Step 10 costs what step 8's ours pass does, step 9 what step 7 does. Measured on the corpus with
24 processes on 24 cores (2026-11-24): step 10 took 1.3 min per ε at 0.5 to 2.1 at 0.05, step
9 4.6 to 8.9 min. Every result document and summary carries `"delta_rule"`: `d*/(1+eps)` or `d*`
for steps 7 and 9, `T/(1+eps)` or `T` for 8 and 10; step 12 labels each tree by it.

**Check:** one summary per trajectory per ε; no `over (1+eps)*delta` line. Step 10's logs say
`N over the target, which delta = T allows` — that is the point of the step, not a fault (89% of
its runs). In step 12, `over 2kappa-2` must be 0 for both: here it tests ours against twice the
baseline, in every pair. On the corpus it is 0 in all 9 041 580 pairs of steps 9 and 10.

## 11. Fix the size: ours at each baseline's size, compare the error

Steps 7–10 hold the error and compare sizes. This step reads step 6's operating points the
other way round, the protocol *fix size, compare error* of [comparison.md](comparison.md):
ours gets **the baseline's own number of vertices**, and the comparison is of error, ours'
measured Fréchet distance against the baseline's d\*.

**Ours cannot be set by size**, so `ssk_size_match` bisects δ until it fits. For each baseline
document of step 6, k is its vertex count: ⌈N/2^m⌉, or a point or two under for SQUISH-E and
DOTS, so ours never gets more than the baseline had. The search keeps a bracket [lo, hi], with
ours at hi within k vertices and at lo not. It starts at lo = 0, never run, and hi = the
diagonal of the trajectory's bounding box. Every segment in the box is within that of the
trajectory, so κ is 2 and ours has 2 vertices; hi doubles if the margin `run` takes off δ ever
says otherwise. It halves the bracket until hi − lo is 10⁻⁹ of the diagonal, the tolerance
every distance is measured to, about 31 runs of ours. Ours at hi is then measured as step 7
measures it, and its bound (1 + ε)·hi is checked exactly.

The search depends on the trajectory and k, not on the baseline, so the budgets of one
trajectory share a cache of runs. Two baselines at the same size get the same run of ours, and
budgets a point apart share all of the bisection above the δ where their sizes part.
**Bisection assumes ours' size never grows with δ, which is not a theorem.** Where it fails, the
search stops at a boundary, so hi is the top of a bracket, not proven the least δ at which ours
fits. ε is ours alone and no part of the search; the comparison runs at ε = 0.1.

It needs the build and step 6's `data/frechet-distances/`, and nothing of steps 7–10.

```powershell
$ours = "greedy-stabbing-shared"
New-Item -ItemType Directory -Force data\size-matched | Out-Null
$procs = foreach ($ds in (Get-ChildItem data\trajectories -Directory).Name) {
  0..7 | ForEach-Object {
    Start-Process .\build\src\pipelines\size-match\ssk_size_match.exe -NoNewWindow -PassThru `
      -ArgumentList "--in data\trajectories\$ds --ours $ours --eps 0.1 --shard $_ --shards 8 --skip-existing" `
      -RedirectStandardOutput "data\size-matched\$ours-$ds-eps0.1-shard$_.log"
  }
}
$procs | Wait-Process
```

```sh
ours=greedy-stabbing-shared
mkdir -p data/size-matched
for ds in data/trajectories/*/; do
  seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/size-match/ssk_size_match \
    --in "$ds" --ours $ours --eps 0.1 --shard {} --shards 8 --skip-existing
done 2>&1 | tee -a data/size-matched/$ours.log
```

**What it writes**, in step 7's layout:

| path | what |
|---|---|
| `data/size-matched/<ours>/eps<E>/<algorithm>/<dataset>/m<rate>/<name>.json` | ours at a baseline's size: its points, `params` (the δ found), `stats` (`budget`, `output_size`, `algorithm_runs`), its measured `frechet` with `bound` and `within_bound`, and the baseline document it was fitted to under `matched` |
| `data/size-matched/<ours>/eps<E>/summary/<dataset>/<name>.json` | ours beside each baseline run, the shape of step 7's summary, `delta_rule` `size` |

It resumes with `--skip-existing`, as the summary is written last. Each shard ends with a line
like

```
greedy-stabbing-shared geolife eps=0.1: 778 trajectories -> 14405 fits, 285241 runs of ours, in 785.2 s
```

**Cost**, measured on the corpus with 24 processes on 24 cores (2026-11-27): **24 minutes**,
GeoLife 14 of them. 857 456 fits took 16.4 million runs of ours, 19.1 per fit rather than 31,
as budgets share their searches. On disk it is step 7's shape at one ε: one result per baseline
document and one summary per trajectory.

**Check before reading the results:** one summary per trajectory; no `over (1+eps)*delta`
line; ours' `output_size` at most the baseline's in every pair. In step 12, `over 2kappa-2` must
be 0 here too. `over the baseline's error` is the result, not a fault.

**The precision floor.** In 5 925 pairs the baseline's d\* is within its tolerance, an error of
0 within measurement (step 7's *precision floor*: dropped repeats, or runs collinear to the
bit). Ours ends within the two tolerances in 5 839 of them. A ratio of two numbers under the tolerance is noise, so step 12 counts those pairs as
the same distance, 100%; the other 86 keep their ratio, over 200%.

## 12. Read the results

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

### Ours against the baselines

`scripts/comparison_stats.py` reads only step 7's summaries and uses only the standard library,
so it needs no venv and runs in seconds per dataset:

```powershell
python scripts\comparison_stats.py --matched data\error-matched | Out-File -Encoding utf8 data\error-matched\comparison.md
```

```sh
python scripts/comparison_stats.py --matched data/error-matched > data/error-matched/comparison.md
```

**It prints Markdown, per dataset:**

- coverage;
- each baseline's median error at a fixed rate;
- per ε, ours at each baseline's error: the size ratio, the error ratio and dominance in both
  directions.

It then prints the matched table pooled over datasets, and **Theorem 3, checked**. Everything,
plus the rate–distortion curves, also goes as JSON to `data/error-matched/comparison.json`
(change it with `--out`). Every column is described in [algorithm.md](algorithm.md) §6.

**Read the Theorem 3 table first.** Its `over (1+eps)delta` and `over 2kappa-2` columns must be
0, and anything else is a bug.

**Then read "smaller", not "ours dominates".** Both errors are bisection tops, so at ε = 0.1
and 0.5, where ours' error reaches the baseline's exactly, the measured comparison is a coin
toss: on `mot-dancetrack`, "ours dominates" ranges from 31% to 98%. Every pair has already
passed the exact check `d_F(ours) ≤ d*`, so "smaller" is the dominance figure. At ε = 0.3 the
grid's farthest corner is at 0.9959·(1+ε)δ, so ours' error stays below the baseline's and the
two columns agree.

The `mot-dancetrack` results, bucketed by size ratio, are in
`archive/[2026-10-21] Meeting Material - Benchmark data.md`.

Step 9's tree has step 7's shape, so the same script reads it; its Theorem 3 table is the size
bound against κ(d*):

```powershell
python scripts\comparison_stats.py --matched data\error-matched-unscaled | Out-File -Encoding utf8 data\error-matched-unscaled\comparison.md
```

```sh
python scripts/comparison_stats.py --matched data/error-matched-unscaled > data/error-matched-unscaled/comparison.md
```

### Pair by pair, in buckets: all five comparisons

`scripts/ratio_buckets.py` reads the summaries of steps 7–11 and, for every pair, takes two
ratios of ours to the baseline, in percent:

- **size**: ours' vertices over the baseline's. Under 100, ours is smaller.
- **distance**: ours' measured Fréchet distance over the baseline's. How much of the error the
  baseline spent ours spent; in step 7 it cannot pass 100 by much, as ours is held to d*. In
  step 11 it is the result, and the size ratio is at most 100 by construction.

Each is counted into the buckets [0, 10), [10, 20), …, [190, 200) and 200+, separately for every
(algorithm, level, ε): the level is m for steps 7, 9 and 11 (*fixed-rate*,
*fixed-rate-unscaled*, *fixed-size*), f for steps 8 and 10 (*fixed-error*,
*fixed-error-unscaled*), each tree labelled by its summaries' `delta_rule`. A ratio of exactly 100, the same size, falls in
[100, 110).

It first prints **a check table** per comparison and ε: pairs over ours' bound (1+ε)δ, which
must be 0; pairs over the baseline's error, 0 in steps 7–8, nearly all in 9–10, and in 11
whatever the comparison found; and pairs
over Theorem 3's size bound 2κ − 2, which must be 0, κ bounded by the input and by every
baseline of the same summary within δ.

```powershell
python scripts\ratio_buckets.py | Out-File -Encoding utf8 data\ratio-buckets.md
```

```sh
python scripts/ratio_buckets.py > data/ratio-buckets.md
```

Then the ten tables pooled over datasets — size and distance, for each comparison — with
each row's `n`, median and share under 100%. `data/ratio-buckets.json` (`--out`) has every
bucket's count as well as its share, per dataset as well as pooled, and the three check counts.

### Percentiles, one file per experiment

`scripts/ratio_percentiles.py` reads the same summaries with the same ratios and writes
[docs/tables/](tables/) — `experiment-1-fixed-rate.md` to `experiment-5-fixed-size.md`, one per
step 7–11. Each has a table per ε and level, a row per baseline, with the p5, p10, p20, …, p90,
p95, p99 and worst ratio (nearest rank), pooled over the datasets: the size ratio for steps
7–10, the distance ratio for step 11. Only what is between each file's two markers is
rewritten; the description, data sources, goal and caveats around them are written by hand.
Standard library only; about 30 seconds on 24 cores.

```powershell
python scripts\ratio_percentiles.py
```

```sh
python scripts/ratio_percentiles.py
```

### Figures

`scripts/plot_results.py` draws every row of those tables as a histogram, 840 of them, and
nine summary figures that carry the result in fewer pictures: 01–05 for steps 7–8, 06 setting
them beside 9–10, 07–08 for 9–10, and 09 for 11, the median distance ratio by rate. The
write-up's figures, `data/plots/writeup/01_fixed-rate.pdf` to `05_fixed-size.pdf`, are a
histogram per step pooled over the levels, of the size ratio for steps 7–10 and of the distance
ratio for step 11; `06_all-experiments.pdf` puts the five side by side in the text's width, a
column per step and a row per baseline. It is the one script that needs a package, matplotlib,
which the venv has; seconds on 24 cores.

```powershell
.\.venv\Scripts\python.exe scripts\plot_results.py [--eps 0.1]
```

```sh
.venv/bin/python scripts/plot_results.py [--eps 0.1]
```

Histograms go to `data/plots/histograms/<comparison>/<size|distance>/<ours>/`, the summary to
`data/plots/summary/`. What they show, which to present, and what to say alongside them is
[results.md](results.md).

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
| set `PYTHONPATH` | `$env:PYTHONPATH = "$PWD\src"` | `export PYTHONPATH=$PWD/src` |
| configure CMake | add `-G "MinGW Makefiles"` | no generator flag needed |
| run a binary | `.\build\...\ssk_simplify.exe` | `./build/.../ssk_simplify` |
| download | `curl.exe` (bare `curl` is an alias) | `curl` |
| loop over datasets | `foreach ($ds in Get-ChildItem ... -Directory)` | `for ds in data/trajectories/*/; do … done` |
| run N shards at once | `Start-Process ... -PassThru`, then `Wait-Process` | `xargs -P N` |
| line continuation | backtick `` ` `` | backslash `\` |

Two Windows-only annoyances worth knowing, neither of which affects the code: MinGW-built
binaries do not run from Git Bash because its `PATH` lacks the runtime DLLs — use PowerShell —
and a filename cannot contain `:`, which is why the incident reports in `archive/` use a dash.

The project has been built and run end to end on Windows with MinGW g++ 15.1. Nothing in it
should stop it building on Linux with GCC 11+ or Clang 14+, but that has not been tried, so
treat it as "no known obstacles" rather than "verified".
