# Streaming simplification of trajectory streams

Final-year project. Streaming δ-simplification of polygonal curves under the continuous
Fréchet distance, plus the tooling around it: dataset ingestion, an I/O layer, and
visualisation.

The algorithm itself is being rewritten. What is in `src/` today is everything that does
not depend on which algorithm is used.

## Layout

```
src/
  trajio/     Python   dataset formats  ->  our format          (ingestion)
  algo/       C++      the algorithm library; reads our format
    io/                json.{hpp,cpp}, trajectory.{hpp,cpp}
    simplify/          the Simplifier interface + Douglas-Peucker, SQUISH, DOTS
  viz/        Python   renders trajectories and results to an image

data/            downloads, converted trajectories, generated curves, images (git-ignored)
docs/            datasets.md, trajio.md, simplification.md
scripts/         entry points and experiment drivers
tests/           C++ unit tests

CMakeLists.txt   the C++ build
requirements.txt the Python dependencies -- matplotlib, for src/viz only
UPDATES.md       what changed, and why
CLAUDE.md        conventions, for agents working in this repo
```

Three pieces, one format between them: **trajio writes it, `algo/io` reads it, `viz` draws
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
sequence), as declared by `t_unit`. It is carried because the SED-based baselines in the
literature need a clock, while Fréchet-bounded and perpendicular-distance simplification do
not. A *result* document is the same plus `algorithm`, `mode`, `params`, `stats`,
`input_points` (the original curve, with `points` now the simplified one) and, when checked,
`frechet`. `viz/plot.py` renders either — the statistics panel fills in from whatever is
present.

More on the data directories and how to refill them: [`data/README.md`](data/README.md).

## Setup

Prerequisites: a C++20 compiler, CMake ≥ 3.20, and Python 3.10 or later.

### 1. Python environment

`src/viz` needs matplotlib; `src/trajio` needs nothing but the standard library. Use a
virtual environment rather than a system or Anaconda Python — matplotlib and NumPy have to
be a matched pair, and a stale pair fails to import at all.

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

On Windows with MinGW add `-G "MinGW Makefiles"`, and the test binary is
`.\build\tests\ssk_tests.exe`.

### 3. End to end

```powershell
python src\viz\generate.py spiral2d -n 600 -o data\synthetic\spiral2d.json
python src\viz\plot.py data\synthetic\spiral2d.json -o data\renders\spiral2d.png
```

That writes a PNG under `data/renders/`. If you have the datasets downloaded, the same
works on a converted trajectory — see [`data/README.md`](data/README.md).

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

[`docs/trajio.md`](docs/trajio.md), [`docs/datasets.md`](docs/datasets.md),
[`data/downloads/README.md`](data/downloads/README.md).

### `src/algo/` — the algorithm library

Builds as `ssk`, with `src/algo` as the include root, so headers are included by their path
within the library. Today it holds only `io/`; the algorithm tracks land beside it.

```cpp
#include "io/trajectory.hpp"

auto doc = ssk::io::read_trajectory_file("data/trajectories/mopsi/mopsi-000002.json");
// doc.dim, doc.name, doc.points  --  points is vector<vector<double>>
```

`io/json` is a small hand-written JSON value, parser and writer, so the build stays offline
instead of pulling a third-party library. `io/trajectory` is the document layer on top of
it: dimension-agnostic, and strict about a point that disagrees with `dim`.

`simplify/` holds the `Simplifier<D>` interface and three baselines ported from their
reference implementations — Douglas–Peucker, SQUISH and DOTS:

```cpp
#include "simplify/douglas_peucker.hpp"

const auto curve = ssk::simplify::curve_of<2>(doc);
ssk::simplify::DouglasPeucker<2> dp({{"count", 32}});
const auto kept = dp.indices(curve);      // positions into `curve`
```

See [`docs/simplification.md`](docs/simplification.md) for the interface, the parameter
scales, and what was changed from each upstream.

### `src/pipelines/` — experiments

```sh
./build/src/pipelines/simplify-sweep/ssk_simplify \
    --in data/trajectories/mopsi --out data/simplified-trajectories

./build/src/pipelines/frechet-distance/ssk_frechet \
    --in data/trajectories/mopsi --out data/frechet-distances
```

Simplifies every trajectory in a dataset directory with all three baselines, at a fixed
compression rate of 1/2^m for m = 1…6, and writes
`data/simplified-trajectories/<algorithm>/<dataset>/m<rate>/<name>.json` — one document per
trajectory per algorithm per rate.

Each result is a **trajectory document in its own right**: the kept points and their
timestamps, so `io/trajectory` reads it as a curve and `viz/plot.py` draws it with no
rehydration step. On top of that it carries `algorithm`, `mode`, the `params` that produced
it — `count`, `buffer_size`, or the threshold DOTS' search resolved to — `stats`, and
`source` naming the input it came from. `--rates` changes the deepest rate, `--limit` stops
after N trajectories, and `--shard I --shards N` splits one dataset across processes — worth
it on GeoLife, where a handful of 60k-point tracks dominate the run.

`ssk_frechet` then measures how far each simplification is from the trajectory it came from,
writing `data/frechet-distances/<algorithm>/<dataset>/m<rate>/<name>.json` — the same tree,
one document per measurement, carrying the distance alongside the `params` and `stats` that
produced it, so a results table needs only that document.

Budget: about 4.3× the input on disk, and roughly 14 documents per trajectory. How to run a
full corpus sweep, what it costs and how sharding works is
[`docs/pipeline.md`](docs/pipeline.md).

Why a fixed compression rate rather than a fixed tolerance, and what the resulting numbers
can and cannot claim, is [`docs/comparison.md`](docs/comparison.md).

### `src/viz/` — pictures out

```sh
python src/viz/generate.py spiral2d -n 600 -o data/synthetic/spiral2d.json
python src/viz/plot.py data/synthetic/spiral2d.json -o data/renders/spiral2d.png
```

`plot.py` renders with matplotlib on the `Agg` backend, so it writes files without needing
a display. 2D is one panel; 3D is a 3D line plot plus the `xy`, `xz` and `yz` projections,
since a single 3D view hides most of what a simplification did. The statistics panel on the
right fills in from `params`, `stats` and `frechet` and shows `?` for whatever is absent, so
a bare trajectory renders fine.

`generate.py` builds the synthetic test curves: `spiral2d`, `zigzag2d`, `walk2d`,
`u-shape`, `o-shape` (2D) and `shrink-spiral` (3D). `u-shape` and `o-shape` are the two
that stress a simplifier — near-parallel legs whose tolerance balls overlap, and a closed
loop with both endpoints pinned.
