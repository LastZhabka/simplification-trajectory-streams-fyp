# Running the pipelines

`src/pipelines/` holds the experiment drivers, one directory each, and each keeps its logic in
a library beside a thin `main` so the tests can drive it without running the binary:

| Directory | Binary | What |
|---|---|---|
| `simplify-sweep/` | `ssk_simplify` | every baseline's simplification of every trajectory, at every compression rate |
| `frechet-distance/` | `ssk_frechet` | how far each of those simplifications is from the trajectory it came from |
| `audit/` | `ssk_audit` | checks result documents on disk against the trajectories they came from |

The first two run in that order and the second mirrors the first's output tree one for one;
`ssk_audit` is verification, not a stage. What the experiment is *for* is
[comparison.md](comparison.md); this is how to run it and what it costs.

## The sweep

```sh
./build/src/pipelines/simplify-sweep/ssk_simplify \
    --in data/trajectories/mopsi --out data/simplified-trajectories
```

```powershell
.\build\src\pipelines\simplify-sweep\ssk_simplify.exe `
    --in data\trajectories\mopsi --out data\simplified-trajectories
```

| Flag | Default | What |
|---|---|---|
| `--in DIR` | required | a dataset directory of trajectory documents; `index.json` is skipped |
| `--out DIR` | `data/simplified-trajectories` | root of the result tree |
| `--rates N` | 6 | deepest compression rate, `1/2^N` |
| `--limit N` | all | stop after N trajectories |
| `--shard I --shards N` | 0 / 1 | process only every Nth trajectory, starting at I |
| `--no-dots` | off | run DPn and SQUISH only — see *The one trajectory DOTS does not cover* |

Output lands at `<out>/<algorithm>/<dataset>/m<rate>/<name>.json`, one document per
trajectory per algorithm per rate. The dataset name comes from the input directory, so
`--in data/trajectories/mopsi` writes under `mopsi/`.

## Sharding

The driver is single-threaded, and trajectories are independent, so the way to use more than
one core is to run several processes over disjoint slices of the file list. `--shard I
--shards N` gives a worker positions `I, I+N, I+2N, …` of the sorted list:

```
position:  0  1  2  3  4  5  6  7  8  9 10 11
shard 0:   ●        ●        ●        ●
shard 1:      ●        ●        ●        ●
shard 2:         ●        ●        ●        ●
shard 3:            ●        ●        ●        ●
```

Every worker derives the same sorted list from the same directory and computes its slice
arithmetically, so there is **no queue, no lock and no master process**, and since each
trajectory writes to its own path, two workers can never touch the same file.

**Shards stride rather than take contiguous blocks**, and that is the whole point. Cost is
not uniform across trajectories — DOTS approaches `n²` in trajectory length, and GeoLife
ranges from 3 points to 92 645 — and the files are ordered by user and folder, so one person
who recorded long tracks produces a *run* of consecutive expensive files. A block split hands
that run to a single worker, which is then still going when the others have finished.

### Subdividing one shard

If one worker falls behind, its slice can be split further without any new flag, because a
shard of N is exactly a union of shards of any multiple of N. Every position that is
`7 mod 64` is also `7 mod 8`, and so is every position that is `15, 23, 31, 39, 47, 55` or
`63 mod 64`, so:

```
shard 7 of 8  =  shards {7, 15, 23, 31, 39, 47, 55, 63} of 64
```

Those eight cover shard 7's files completely and disjointly. Running them as eight processes
subdivides that one shard eight ways using a flag that only knows about the whole dataset.

### Where striding fails

Striding spreads cost *on average*, but it is a fixed arithmetic pattern: two expensive
trajectories collide whenever the gap between their positions is a multiple of the stride.
That happened on the first corpus run. `geolife-013552` (64 483 points) sits at position
13 551 and `geolife-014000` (33 361 points) at 13 999 — **448 apart, which divides by both 8
and 64** — so they shared a worker in the 8-way split and shared one again in the 64-way
re-split, and that worker was the critical path both times.

A shared work queue, or dealing longest-first by the point counts already in `index.json`,
would not have this failure mode. Striding is four lines and needs no coordination; that is
the trade being made.

## What a corpus run costs

Measured on this machine (MinGW g++ 15.1 `-O2`, 8 cores), one pass over all nine datasets at
m = 1…6:

| Dataset | Trajectories | Points | Documents | Wall time |
|---|---:|---:|---:|---:|
| `geolife` | 18 670 | 24 876 978 | 315 620 | sharded 8 ways |
| `mopsi` | 6 779 | 7 850 387 | 114 029 | 23 min |
| `ngsim-us-101` | 2 847 | 4 802 933 | 51 244 | 36 min |
| `ngsim-i-80` | 3 001 | 4 566 387 | 54 007 | 34 min |
| `ngsim-lankershim` | 6 712 | 1 607 319 | 99 281 | 22 min |
| `ngsim-peachtree` | 4 495 | 873 887 | 66 512 | 15 min |
| `mot-mot17` | 2 388 | 614 103 | 33 831 | 3.7 min |
| `mot-mot20` | 2 332 | 1 336 920 | 39 783 | 8.4 min |
| `mot-dancetrack` | 692 | 574 078 | 12 335 | 8.6 min |
| **total** | **47 916** | **47 102 992** | **786 642** | |

Those timings are single-process per dataset, several datasets at once, so they include disk
contention. **DOTS is essentially the whole cost** — measured separately, a 6-rate sweep over
778 766 points takes DPn 0.4 s, SQUISH 2.9 s, and DOTS 76 s once its threshold search is
included.

**GeoLife needs sharding, the others do not.** Its median trajectory is 506 points but five
exceed 56 000, and at DOTS' scaling those five dominate the dataset: run as one process it
managed roughly one trajectory every four minutes once it reached them. See the length
distributions in [datasets.md](datasets.md) for why this is a GeoLife problem specifically.

The document count is predictable rather than emergent: a rate is only produced when its
budget stays at or above 2 vertices, and 5 for SQUISH, so summing the qualifying rates over
all three algorithms gives 786 642 exactly — which every dataset has matched to the document.

## The one trajectory DOTS does not cover

**`geolife/geolife-013552.json`, 64 483 points, has DPn and SQUISH results but no DOTS
results.** The corpus is otherwise complete: DOTS covers 18 669 of GeoLife's 18 670
trajectories and every trajectory in the other eight datasets.

It is not a size limit. Fifteen of the sixteen trajectories over 40 000 points completed
normally, including `geolife-001523` at **92 645 points** — 43% larger — which came through in
the first unsharded run. `geolife-014000` at 33 361 points takes seconds. Whatever makes
`013552` expensive is a property of its geometry keeping DOTS' frontier alive, not its length,
and it is the only trajectory in 47 916 that behaves this way. Left running, its budget search
consumed over 90 CPU-minutes without finishing, four separate times.

To reproduce the state on disk:

```sh
./build/src/pipelines/simplify-sweep/ssk_simplify \
    --in <dir containing only geolife-013552.json, named geolife> \
    --out data/simplified-trajectories --no-dots
```

```powershell
.\build\src\pipelines\simplify-sweep\ssk_simplify.exe `
    --in <dir containing only geolife-013552.json, named geolife> `
    --out data\simplified-trajectories --no-dots
```

**What is understood, and what is not.** The budget search had a real defect, since fixed: not
every budget is reachable, because output size is a step function of the threshold, and here
the m = 1 budget of 32 242 sits above the 31 197 points DOTS emits at *any* threshold — it
always takes some shortcut. The search could therefore never hit the target exactly and
halved its bracket down to a width of `1e-9` in log space, which is meaningless precision for
an integer-valued function. It now stops once the bracket's two ends have gone three halvings
without moving, which is the point at which it is only refining the location of a step it has
already bracketed. `sweep/dots gives up on a budget it cannot reach` pins that, and
`sweep/dots reaches the same budgets as an independent search` pins that the results did not
change.

**That fix did not make `013552` tractable**, so it is not the root cause. Individual DOTS
runs on this trajectory measure 0.2–2.7 s across the thresholds sampled, and the search is
bounded well under a hundred evaluations per rate, which does not add up to 90 minutes. The
remaining explanation has to be that the search visits thresholds far more expensive than any
sampled, or that some evaluation does not terminate in the time expected. **This is not
diagnosed.** Anyone picking it up should trace every evaluation — threshold, resulting size
and elapsed time — rather than trusting the arithmetic above.

## Measuring the error

```sh
./build/src/pipelines/frechet-distance/ssk_frechet --in data/trajectories/mopsi \
    --simplified data/simplified-trajectories --out data/frechet-distances
```

```powershell
.\build\src\pipelines\frechet-distance\ssk_frechet.exe --in data\trajectories\mopsi `
    --simplified data\simplified-trajectories --out data\frechet-distances
```

Writes `data/frechet-distances/<algorithm>/<dataset>/m<rate>/<name>.json`, **mirroring the
simplified tree exactly**, one document per simplification:

```json
{ "algorithm": "dots", "dataset": "mopsi", "dim": 2, "mode": "budget",
  "name": "mopsi/routes/1/1216464589688",
  "frechet": { "distance": 0.00019507700706403058,
               "tolerance": 2.537329637630661e-11, "computer": "dv-gis-cup-2017" },
  "params": { "lssd_threshold": 1.6915935776057353e-07 },
  "stats": { "input_size": 130, "output_size": 17, "compression": 7.647, "m": 3,
             "rate": 8, "target": 17, "algorithm_runs": 4 },
  "source": "mopsi/mopsi-000003.json",
  "simplified": "dots/mopsi/m3/mopsi-000003.json" }
```

`params` and `stats` are carried over from the simplification, so a results table needs only
this document — it does not have to open two more to learn which threshold produced the
number. The geometry is not repeated; `simplified` and `source` name where it lives.

**The dataset is part of the path**, and it has to be: trajectory files are named per source,
not globally, so `ngsim-000001.json` exists in all four NGSIM sites and `mot-000001.json` in
all three MOT sets. Flattening the dataset out of the path would silently overwrite about
19 000 trajectories' results.

**The tolerance is relative** — `--tol`, default `1e-9` — taken against the input's
bounding-box diagonal. It has to be, because the coordinate unit differs by six orders of
magnitude across the corpus: one absolute tolerance is either meaningless on degrees or
ruinous on State Plane feet. The absolute value used is recorded in each document.

`--shard I --shards N` works exactly as in the sweep, and matters for the same reason: the
free-space diagram is `n × m` cells, so cost is quadratic in trajectory length at a fixed rate
and the long GeoLife tracks dominate again.

### Reading the numbers

**Error is not monotone in the vertex budget.** Measured over 40 Mopsi trajectories, 22 of 120
(algorithm, trajectory) series had a *deeper* rate score *lower* than the one above it — for
Douglas–Peucker by as much as 27%, far above the tolerance. That is a property of the measure,
not a defect: a kept vertex forces the matching to pass through it, so a 3-point simplification
can be further from the input than the 2-point chord that skips it. A table of mean error
against rate is well behaved; individual series are not, and nothing downstream should assume
they are.

## Checking a run

Result documents are trajectory documents, so they can be read back with
`io::read_trajectory` and drawn with `viz/plot.py` directly:

```sh
python src/viz/plot.py \
  data/simplified-trajectories/dots/mopsi/m3/mopsi-000003.json -o data/renders/check.png
```

```powershell
python src\viz\plot.py `
  data\simplified-trajectories\dots\mopsi\m3\mopsi-000003.json -o data\renders\check.png
```

Beyond that, a result is verifiable against its input, and `source` names it: the points must
be a subsequence of the input matched on point *and* timestamp — matching on coordinates
alone gives false failures wherever a GPS trace repeats a fix — the endpoints must survive,
`stats.output_size` must not exceed `stats.target`, and for DPn it must equal it.
