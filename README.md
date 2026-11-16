# Streaming simplification of trajectory streams

Final-year project. Streaming δ-simplification of polygonal curves under the continuous
Fréchet distance, and the experiment that compares it with four published simplifiers on
47 916 real trajectories.

What is here:

- **Our algorithm**, `ConvexShadow`: Algorithm 1 of Cheng, Huang & Jiang, *Simplification of
  Trajectory Streams* (SoCG 2025). Its output is within `(1+ε)·δ` of the input and has at most
  `2·opt − 2` vertices.
- **Four baselines** ported from their reference implementations: Douglas–Peucker (DPn),
  SQUISH, SQUISH-E and DOTS.
- **The continuous Fréchet distance**, ported from the GIS Cup 2017 winner, to measure all of
  them the same way.
- **The pipelines** that run the comparison, and the tooling around them: dataset ingestion,
  one JSON format, and visualisation.

## Layout

```
src/
  trajio/       Python   dataset formats  ->  our format          (ingestion)
  algo/         C++      the library `ssk`; reads our format
    io/                  json.{hpp,cpp}, trajectory.{hpp,cpp}
    simplify/            the Simplifier interface, the four baselines, and ours (convex_shadow,
                         and convex_shadow_far_points, its first implementation)
    frechet/             the continuous Fréchet distance
    validate/            checks of our algorithm's invariants, counted as it runs
  pipelines/    C++      the experiment drivers
    simplify-sweep/      ssk_simplify   every baseline at every compression rate
    frechet-distance/    ssk_frechet    the error of each of those simplifications
    error-match/         ssk_match      ours at each baseline's measured error
    audit/               ssk_audit      checks result documents against their inputs
    corpus.{hpp,cpp}     what the four share: file list, sharding, atomic document write
  viz/          Python   renders trajectories and results; the compression-rate report

scripts/        comparison_stats.py, the statistics over what ssk_match wrote
tests/          C++ unit tests, own harness
data/           downloads, trajectories and every result tree (git-ignored)
docs/           one document per topic, listed at the end

CMakeLists.txt   the C++ build
requirements.txt the Python dependencies -- matplotlib and numpy, for src/viz only
UPDATES.md       what changed, and why (local, git-ignored)
CLAUDE.md        conventions, for agents working in this repo
```

One format between all of it: **trajio writes it, `algo/io` reads it, the pipelines extend it,
`viz` draws it.**

### Our format

One JSON document per curve — the only trajectory format in this repository:

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

`dim` is 2 or 3 and every point has exactly that many coordinates; it is inferred from the
first point when the field is absent. `name` is provenance and is optional. `t` is optional
too — one timestamp per point, in `unix_ms` (absolute epoch) or `ms` (from the start of the
sequence), as declared by `t_unit`. It is carried because the SED-based baselines need a
clock, while Fréchet-bounded and perpendicular-distance simplification do not.

A *result* document is a trajectory document of the simplified points, plus `algorithm`,
`mode`, `params`, `stats` and `source`. A measured one adds `frechet`, and a validated run of
ours adds `validation`. It does not repeat the input: `source` names it. `viz/plot.py` renders
either kind.

More on the data directories and how to refill them: [`data/README.md`](data/README.md).

## Setup

Prerequisites: a C++20 compiler, CMake ≥ 3.20, and Python 3.10 or later.

### 1. Python environment

`src/viz` needs matplotlib and numpy; `src/trajio` and `scripts/` need nothing but the
standard library. Use a virtual environment rather than a system or Anaconda Python —
matplotlib and NumPy have to be a matched pair, and a stale pair fails to import at all.

```powershell
# Windows / PowerShell, from the repository root
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

```sh
# Linux / macOS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

`.venv/` is git-ignored. If PowerShell refuses to run the activation script, either
`Set-ExecutionPolicy -Scope Process RemoteSigned` for that session, or skip activation and
call `.\.venv\Scripts\python.exe` directly.

`trajio` is imported as a package, so `src/` has to be on the import path:

```powershell
$env:PYTHONPATH = "$PWD\src"      # PowerShell
export PYTHONPATH="$PWD/src"      # bash
```

Check it:

```sh
python -m trajio selftest      # every dataset parser, no downloads needed
```

### 2. C++ build

No third-party libraries; the build works offline.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/ssk_tests
```

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
cmake --build build -j
.\build\tests\ssk_tests.exe
```

The generator is the only difference: MinGW needs naming, the default works elsewhere. 82
tests, a few seconds.

### 3. End to end

A picture, needing no dataset:

```powershell
python src\viz\generate.py spiral2d -n 600 -o data\synthetic\spiral2d.json
python src\viz\plot.py data\synthetic\spiral2d.json -o data\renders\spiral2d.png
```

```sh
python src/viz/generate.py spiral2d -n 600 -o data/synthetic/spiral2d.json
python src/viz/plot.py data/synthetic/spiral2d.json -o data/renders/spiral2d.png
```

The whole experiment, from downloads to the comparison tables, is
[`docs/running.md`](docs/running.md): every step in PowerShell and POSIX form, with what each
costs.

## The pieces

### `src/trajio/` — datasets in

Parses GeoLife, Mopsi, NGSIM and MOT in their native formats and writes one JSON document
per trajectory, values copied verbatim: no projection, no unit conversion, no resampling.
Run it with `src/` on the import path:

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m trajio sources
python -m trajio export --source mopsi --root data\downloads\mopsi --dims 2 `
                        --out data\trajectories\mopsi
python -m trajio selftest
```

```sh
export PYTHONPATH=$PWD/src
python -m trajio sources
python -m trajio export --source mopsi --root data/downloads/mopsi --dims 2 \
                        --out data/trajectories/mopsi
python -m trajio selftest
```

[`docs/trajio.md`](docs/trajio.md), [`docs/datasets.md`](docs/datasets.md),
[`data/downloads/README.md`](data/downloads/README.md).

### `src/algo/` — the library

Builds as `ssk`, with `src/algo` as the include root, so headers are included by their path
within the library.

```cpp
#include "io/trajectory.hpp"

auto doc = ssk::io::read_trajectory_file("data/trajectories/mopsi/mopsi-000002.json");
// doc.dim, doc.name, doc.points  --  points is vector<vector<double>>
```

**`io/`** is a small hand-written JSON value, parser and writer, so the build stays offline,
and the trajectory document layer on top of it.

**`simplify/`** holds the `Simplifier<D>` interface, the four baselines and ours:

```cpp
#include "simplify/convex_shadow.hpp"
#include "simplify/douglas_peucker.hpp"

const auto curve = ssk::simplify::curve_of<2>(doc);

ssk::simplify::DouglasPeucker<2> dp({{"count", 32}});
const auto kept = dp.indices(curve);            // positions into `curve`

ssk::simplify::ConvexShadow ours({{"delta", 1e-4}, {"eps", 0.1}});
const auto simplified = ours.simplify(curve);   // points; ours invents its vertices
```

| | Source | Parameter |
|---|---|---|
| Douglas–Peucker (DPn) | psimpl v7 | `count`, a vertex budget |
| SQUISH | caoweiquan322/dots | `buffer_size`, a vertex budget |
| SQUISH-E | uestc-db/traj-compression | `cmp_ratio` and `sed_error` |
| DOTS | caoweiquan322/dots | `lssd_threshold` |
| ours, `ConvexShadow` | Cheng, Huang & Jiang, SoCG 2025 | `delta` and `eps` |
| ours, `ConvexShadowFarPoints` | its first implementation, kept for comparison | `delta` and `eps` |

The ports follow their sources statement for statement and are kept in that form. The one
change to an upstream algorithm is in DOTS, which could join a point to a vertex that comes
after it and return points out of order; the port does not.
[`docs/simplification.md`](docs/simplification.md) has the interface, what was changed from
each upstream and how each port was checked against it.
[`docs/algorithm.md`](docs/algorithm.md) describes ours.

**`frechet/`** is the continuous Fréchet distance. The decision `within(a, b, r)` is the
primitive, and `distance` bisects it. [`docs/frechet-distance.md`](docs/frechet-distance.md).

**`validate/`** checks the invariants our guarantee rests on while the algorithm runs — the
grid polygon, every shadow, every emitted segment — and counts the outcomes instead of
throwing, so a long run reports how many checks failed. It is off unless a validator is
passed, because it costs 2–3.5×.

### `src/pipelines/` — the experiment

Four binaries, run in this order, each reading the tree the one before wrote:

| Binary | Reads | Writes |
|---|---|---|
| `ssk_simplify` | `data/trajectories/<dataset>` | `data/simplified-trajectories/<algorithm>/<dataset>/m<rate>/<name>.json` |
| `ssk_frechet` | the two above | `data/frechet-distances/<algorithm>/<dataset>/m<rate>/<name>.json` |
| `ssk_match` | trajectories and distances | `data/error-matched/eps<E>/{<algorithm>,summary}/<dataset>/...` |
| `ssk_audit` | trajectories and simplifications | nothing; reports documents that are not a subsequence of their input |

```sh
./build/src/pipelines/simplify-sweep/ssk_simplify   --in data/trajectories/mopsi
./build/src/pipelines/frechet-distance/ssk_frechet  --in data/trajectories/mopsi
./build/src/pipelines/error-match/ssk_match         --in data/trajectories/mopsi --eps 0.5 --validate
python scripts/comparison_stats.py --matched data/error-matched
```

```powershell
.\build\src\pipelines\simplify-sweep\ssk_simplify.exe   --in data\trajectories\mopsi
.\build\src\pipelines\frechet-distance\ssk_frechet.exe  --in data\trajectories\mopsi
.\build\src\pipelines\error-match\ssk_match.exe         --in data\trajectories\mopsi --eps 0.5 --validate
python scripts\comparison_stats.py --matched data\error-matched
```

**`ssk_simplify`** runs every baseline on every trajectory at a fixed compression rate of
`1/2^m` for m = 1…6. DPn and SQUISH take the budget directly, SQUISH-E takes the rate, and
DOTS has its threshold searched for. Each result is a trajectory document in its own right.

**`ssk_frechet`** measures how far each simplification is from its input, in a tree that
mirrors the first one document for document.

**`ssk_match`** runs ours once per baseline document, at `δ = d*/(1+ε)` where `d*` is that
baseline's measured error. Ours is then no worse than the baseline, and the question is how
many vertices it needs. It measures ours in the same pass, checks the bound exactly, and with
`--validate` counts the invariant checks into every result.

**`scripts/comparison_stats.py`** turns the summaries into the tables: coverage, each
baseline's error at a fixed rate, ours against each baseline, Theorem 3 checked, and the
validation counts. `src/viz/compression.py` draws the distribution of our compression rate.

The first three take `--shard I --shards N` to split a dataset across processes and
`--limit N` to stop early; `ssk_frechet` and `ssk_match` take `--skip-existing` to resume. A document is
written under a temporary name and renamed, so one on disk is always whole.

The flags and layout are in [`docs/pipeline.md`](docs/pipeline.md). Why a fixed compression
rate rather than a fixed tolerance, and what the numbers can and cannot claim, is
[`docs/comparison.md`](docs/comparison.md).

### `src/viz/` — pictures out

```sh
python src/viz/generate.py spiral2d -n 600 -o data/synthetic/spiral2d.json
python src/viz/plot.py data/synthetic/spiral2d.json -o data/renders/spiral2d.png
```

```powershell
python src\viz\generate.py spiral2d -n 600 -o data\synthetic\spiral2d.json
python src\viz\plot.py data\synthetic\spiral2d.json -o data\renders\spiral2d.png
```

`plot.py` renders with matplotlib on the `Agg` backend, so it writes files without needing
a display. 2D is one panel; 3D is a 3D line plot plus the `xy`, `xz` and `yz` projections,
since a single 3D view hides most of what a simplification did. The statistics panel lists
the document's `params` and `stats`, and its `frechet` block when it has one, so a bare
trajectory renders fine.

`generate.py` builds the synthetic test curves: `spiral2d`, `zigzag2d`, `walk2d`,
`u-shape`, `o-shape` (2D) and `shrink-spiral` (3D). `u-shape` and `o-shape` are the two
that stress a simplifier — near-parallel legs whose tolerance balls overlap, and a closed
loop with both endpoints pinned.

## Documents

| | |
|---|---|
| [`docs/running.md`](docs/running.md) | every step from a fresh clone to the comparison tables, with costs |
| [`docs/algorithm.md`](docs/algorithm.md) | our algorithm, `ssk_match`, the statistics and the validation |
| [`docs/comparison.md`](docs/comparison.md) | the experimental protocol and what it can claim |
| [`docs/simplification.md`](docs/simplification.md) | the baselines: interface, ports, deviations, verification |
| [`docs/frechet-distance.md`](docs/frechet-distance.md) | the Fréchet port and its measurements |
| [`docs/pipeline.md`](docs/pipeline.md) | the drivers' flags, sharding and costs |
| [`docs/datasets.md`](docs/datasets.md), [`docs/trajio.md`](docs/trajio.md) | the datasets and how they are read |
| [`docs/references.md`](docs/references.md) | the papers and the code this builds on |
