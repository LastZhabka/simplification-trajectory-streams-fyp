# Streaming simplification of trajectory streams

Final-year project. Streaming δ-simplification of polygonal curves under the continuous
Fréchet distance, and the experiment that compares it with four published simplifiers on
47 916 real trajectories.

What is here:

- **Our algorithms**, the greedy stabbers of *Faster δ-simplification of Trajectory Streams*,
  after Cheng, Huang & Jiang, *Simplification of Trajectory Streams* (SoCG 2025). Their output is
  within `(1+ε)·δ` of the input and has at most `2·opt − 2` vertices. `GreedyStabbingShared`
  takes O(ε⁻¹) time per vertex.
- **Four baselines** ported from their reference implementations: Douglas–Peucker (DPn),
  SQUISH, SQUISH-E and DOTS.
- **The continuous Fréchet distance**, ported from the GIS Cup 2017 winner, to measure all of
  them the same way.
- **The pipelines** that run the comparison, and the tooling around them: dataset ingestion
  and one JSON format.

## Layout

```
src/
  trajio/       Python   dataset formats  ->  our format          (ingestion)
  algo/         C++      the library `ssk`; reads our format
    io/                  json.{hpp,cpp}, trajectory.{hpp,cpp}
    geometry/            plane geometry for ours: convex polygons, cones, clipping, the hull
    simplify/            the Simplifier interface, the four baselines, and ours: the greedy
                         stabbers, greedy_stabber.hpp and greedy_stabbing.{hpp,cpp}
    frechet/             the continuous Fréchet distance
    validate/            invariant checks, counted as they run; no stabber is watched yet
  pipelines/    C++      the experiment drivers
    simplify-sweep/      ssk_simplify   every baseline at every compression rate
    frechet-distance/    ssk_frechet    the error of each of those simplifications
    error-match/         ssk_match      ours at each baseline's measured error
    audit/               ssk_audit      checks result documents against their inputs
    corpus.{hpp,cpp}     what the four share: file list, sharding, atomic document write

scripts/        comparison_stats.py, the statistics over what ssk_match wrote
tests/          C++ unit tests, own harness
data/           downloads, trajectories and every result tree (git-ignored)
docs/           one document per topic, listed at the end

CMakeLists.txt   the C++ build
UPDATES.md       what changed, and why (local, git-ignored)
CLAUDE.md        conventions, for agents working in this repo
```

One format between all of it: **trajio writes it, `algo/io` reads it, the pipelines extend
it.**

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
ours adds `validation`. It does not repeat the input: `source` names it.

More on the data directories and how to refill them: [`data/README.md`](data/README.md).

## Setup

Prerequisites: a C++20 compiler, CMake ≥ 3.20, and Python 3.10 or later.

### 1. Python

`src/trajio` and `scripts/` need nothing but the standard library, so there is nothing to
install. `trajio` is imported as a package, so `src/` has to be on the import path:

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

The generator is the only difference: MinGW needs naming, the default works elsewhere. 85
tests, a few seconds.

### 3. End to end

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
#include "simplify/douglas_peucker.hpp"
#include "simplify/greedy_stabbing.hpp"

const auto curve = ssk::simplify::curve_of<2>(doc);

ssk::simplify::DouglasPeucker<2> dp({{"count", 32}});
const auto kept = dp.indices(curve);            // positions into `curve`

ssk::simplify::GreedyStabbingShared ours({{"delta", 1e-4}, {"eps", 0.1}});
const auto simplified = ours.simplify(curve);   // points; ours invents its vertices
```

| | Source | Parameter |
|---|---|---|
| Douglas–Peucker (DPn) | psimpl v7 | `count`, a vertex budget |
| SQUISH | caoweiquan322/dots | `buffer_size`, a vertex budget |
| SQUISH-E | uestc-db/traj-compression | `cmp_ratio` and `sed_error` |
| DOTS | caoweiquan322/dots | `lssd_threshold` |
| ours, `GreedyStabbingShared` | Algorithm 2 of the write-up: a regular m-gon, starts on the circle of radius δ, computations shared, O(ε⁻¹) per vertex | `delta` and `eps` |
| ours, `GreedyStabbingNaive` | Algorithm 1: the same without the sharing, O(ε^(−3/2)) per vertex | `delta` and `eps` |
| ours, `GreedyStabbingNaivePaper` | Algorithm 1 with the paper's `conv(G_v)` and every grid point as a start | `delta` and `eps` |
| ours, `GreedyStabbingNaivePaperBoundary` | the same with the boundary cells' grid points only | `delta` and `eps` |

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
| `ssk_match` | trajectories and distances | `data/error-matched/<ours>/eps<E>/{<algorithm>,summary}/<dataset>/...` |
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
`1/2^m` for m = 3…8 (`--rates`). DPn and SQUISH take the budget directly, SQUISH-E takes the rate, and
DOTS has its threshold searched for. Each result is a trajectory document in its own right.

**`ssk_frechet`** measures how far each simplification is from its input, in a tree that
mirrors the first one document for document.

**`ssk_match`** runs ours once per baseline document, at `δ = d*/(1+ε)` where `d*` is that
baseline's measured error. Ours is then no worse than the baseline, and the question is how
many vertices it needs. It measures ours in the same pass, checks the bound exactly, and with
`--validate` counts the invariant checks into every result. `--ours` picks the implementation,
`greedy-stabbing-shared` by default, `greedy-stabbing-naive`, `greedy-stabbing-naive-paper` or
`greedy-stabbing-naive-paper-boundary`, and each gets its own tree. With `--validate` the
counts are zero for now: no stabber is watched yet.

**`scripts/comparison_stats.py`** turns the summaries into the tables: coverage, each
baseline's error at a fixed rate, ours against each baseline, Theorem 3 checked, and the
validation counts, with each implementation of ours at each ε as its own run.

The first three take `--shard I --shards N` to split a dataset across processes and
`--limit N` to stop early; `ssk_frechet` and `ssk_match` take `--skip-existing` to resume. A document is
written under a temporary name and renamed, so one on disk is always whole.

The flags and layout are in [`docs/pipeline.md`](docs/pipeline.md). Why a fixed compression
rate rather than a fixed tolerance, and what the numbers can and cannot claim, is
[`docs/comparison.md`](docs/comparison.md).

## Documents

| | |
|---|---|
| [`docs/running.md`](docs/running.md) | every step from a fresh clone to the comparison tables, with costs |
| [`docs/algorithm.md`](docs/algorithm.md) | our algorithm, `ssk_match`, the statistics and the validation |
| [`docs/convex-shadow-32.md`](docs/convex-shadow-32.md) | the m-gon region, and why it has O(ε^(−1/2)) sides; the greedy stabbers in §10–11 |
| [`docs/comparison.md`](docs/comparison.md) | the experimental protocol and what it can claim |
| [`docs/simplification.md`](docs/simplification.md) | the baselines: interface, ports, deviations, verification |
| [`docs/frechet-distance.md`](docs/frechet-distance.md) | the Fréchet port and its measurements |
| [`docs/pipeline.md`](docs/pipeline.md) | the drivers' flags, sharding and costs |
| [`docs/datasets.md`](docs/datasets.md), [`docs/trajio.md`](docs/trajio.md) | the datasets and how they are read |
| [`docs/references.md`](docs/references.md) | the papers and the code this builds on |
