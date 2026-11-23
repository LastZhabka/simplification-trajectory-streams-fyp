# CLAUDE.md

Streaming δ-simplification of trajectory streams under the continuous Fréchet distance.
FYP. Ours are the greedy stabbers in `src/algo/simplify/greedy_stabbing.{hpp,cpp}`, each a
`GreedyStabber` (`greedy_stabber.hpp`), whose `run` is the greedy loop they share:
`GreedyStabbingNaive` and `GreedyStabbingShared` are Algorithms 1 and 2 of the write-up, O(ε^(−3/2))
and O(ε⁻¹) per vertex, and `GreedyStabbingNaivePaper` and `GreedyStabbingNaivePaperBoundary` are
Algorithm 1 with the paper's `conv(G_v)` and grid starts (`docs/convex-shadow-32.md` §10–11). The
rest is the baselines they are compared with and the pipelines that run the comparison. The
convex-shadow implementations are gone (UPDATES 65): `docs/algorithm.md` and
`docs/convex-shadow-32.md` §1–9 describe them as of commit 902a3c5.

## Style

- Short, simple, direct. No cleverness, no layers that do not earn their place.
- **No comments.** Names carry the meaning. Comment only a non-obvious *why*.
- **No defensive `try`/`except`.** Catch only what you can actually handle. Let it crash.
- No speculative abstraction — no options, hooks or base classes without a second caller.
- Match the file you are editing.
- Log every change in `UPDATES.md` as a new numbered entry.

## Map

| Path | What | Why it is like this |
|---|---|---|
| `src/trajio/` | Python. Dataset formats → our format. | Standard library only. Readers must run anywhere, and dataset parsing must not depend on the plotting stack. Values are copied verbatim — filtering is a modelling decision, not a reader's. |
| `src/algo/` | C++20 library `ssk`. Reads our format. | No third-party libraries; the build works offline, hence the hand-rolled JSON in `io/json`. `src/algo` is the include root, so headers are `#include "io/..."`. Algorithm tracks land beside `io/`. |
| `src/viz/` | Python. Renders trajectories and results. | The only part that needs installed packages (matplotlib, `Agg` backend). |
| `data/` | `downloads/` `trajectories/` `simplified-trajectories/` `frechet-distances/` `error-matched/` `synthetic/` `renders/` | All git-ignored. Gigabytes, and licence-restricted. Commands are committed, data is not. |
| `tests/` | C++ unit tests. | Own ~60-line harness — GoogleTest via FetchContent would break an offline configure. |
| `src/algo/frechet/` | Continuous Frechet distance, `Frechet<D>`. | Ported from the GIS Cup 2017 winner (MIT, standard library only), cross-verified against two independent implementations. The decision `within()` is the primitive and is exact; `distance()` bisects it. |
| `src/algo/simplify/` | `Simplifier<D>` + Douglas-Peucker, SQUISH, SQUISH-E, DOTS, and ours, the greedy stabbers. | Ports of the reference implementations, not reinventions — match upstream behaviour and record any deviation in `docs/simplification.md`. **Keep them in upstream's form: no restructuring or clean-up of the ports or of `frechet/`.** Dimension is a template parameter; `simplify()` returns points because a future algorithm may invent them, and `SubsetSimplifier` adds `indices()` for the ones that don't. |
| `src/algo/geometry/` | Plane geometry for the greedy stabbers: `Vec`, `ConvexPolygon` as unit normals and offsets, cones from a point, clipping, the hull. | The stabbers depend on it and on `GreedyStabber`, never on `validate` or on one another. |
| `src/algo/validate/` | `validate::Validator`: checks a shadow-keeping algorithm's invariants as it runs and counts the outcomes. Written for the convex-shadow implementations; no stabber is watched yet, so `--validate` counts zero. | Counts, never throws — a corpus run must report how many checks failed, not stop at the first. Its predicates are its own, not the algorithm's, so the two can disagree. Off unless a validator is passed (`ssk_match --validate`), because it costs 2–3.5×. |
| `docs/` | `datasets.md`, `trajio.md`, `simplification.md`, `algorithm.md`, `convex-shadow-32.md`, `comparison.md`, `pipeline.md`, `frechet-distance.md`, `references.md`, `running.md` | Measured numbers, not quoted ones. |
| `src/pipelines/` | Experiment drivers, one directory each: `simplify-sweep` (`ssk_simplify`), `frechet-distance` (`ssk_frechet`), `error-match` (`ssk_match`), `audit` (`ssk_audit`). `corpus.{hpp,cpp}` is what they share: file list, sharding, atomic document write. | The algorithms are C++, so the driver is too — a Python one would spawn a process per trajectory. Entry points live with their piece under `src/`, as `viz` and `trajio` already do. Each keeps its logic in a library beside a thin `main`, so the tests drive it without running the binary. |

**One format everywhere:** `{"dim": 2, "name": "...", "t_unit": "unix_ms", "t": [...],
"points": [[x, y], ...]}`. trajio writes it, `algo/io` reads it, `viz` draws it. `t` is
optional (one per point; `unix_ms` absolute or `ms` from the sequence start) and exists for
the SED-based baselines. A *result* document adds `algorithm`, `mode`, `params`, `stats`,
`source`, and `frechet` once measured; ours adds `validation` when run with `--validate`. No CSV output anywhere — the CSV in the repo is dataset input
only.

## Running

Setup once: `python -m venv .venv`, activate, `python -m pip install -r requirements.txt`.
`trajio` needs `PYTHONPATH=src`.

```powershell
$env:PYTHONPATH = "$PWD\src"

python -m trajio sources | describe <name> | selftest
python -m trajio export --source mopsi --root data\downloads\mopsi --dims 2 `
                        --out data\trajectories\mopsi
python -m trajio stats  --source mopsi --root data\downloads\mopsi

python src\viz\generate.py spiral2d -n 600 -o data\synthetic\spiral2d.json
python src\viz\plot.py data\synthetic\spiral2d.json -o data\renders\spiral2d.png

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles"
cmake --build build -j
.\build\tests\ssk_tests.exe [name-filter]
.\build\src\pipelines\simplify-sweep\ssk_simplify.exe --in data\trajectories\mopsi --out data\simplified-trajectories
.\build\src\pipelines\frechet-distance\ssk_frechet.exe --in data\trajectories\mopsi --out data\frechet-distances
.\build\src\pipelines\error-match\ssk_match.exe --in data\trajectories\mopsi --eps 0.5 [--ours NAME] [--validate]
python scripts\comparison_stats.py --matched data\error-matched
```

```sh
export PYTHONPATH=$PWD/src

python -m trajio sources | describe <name> | selftest
python -m trajio export --source mopsi --root data/downloads/mopsi --dims 2 \
                        --out data/trajectories/mopsi
python -m trajio stats  --source mopsi --root data/downloads/mopsi

python src/viz/generate.py spiral2d -n 600 -o data/synthetic/spiral2d.json
python src/viz/plot.py data/synthetic/spiral2d.json -o data/renders/spiral2d.png

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tests/ssk_tests [name-filter]
./build/src/pipelines/simplify-sweep/ssk_simplify --in data/trajectories/mopsi --out data/simplified-trajectories
./build/src/pipelines/frechet-distance/ssk_frechet --in data/trajectories/mopsi --out data/frechet-distances
./build/src/pipelines/error-match/ssk_match --in data/trajectories/mopsi --eps 0.5 [--ours NAME] [--validate]
python scripts/comparison_stats.py --matched data/error-matched
```

**Nothing in the project is platform-specific — only the shell is.** No `_WIN32`, no `.exe` in
CMake, every Python file opened with an explicit `encoding="utf-8"`, `.gitattributes` pinning
`eol=lf`, and `source`/`simplified` fields in result documents built with `/` rather than
`std::filesystem::path` so the JSON moves between machines. Keep it that way: when adding a
command to a doc, add both forms. The translation table is at the end of `docs/running.md`.

`trajio selftest` runs every dataset parser against fixtures and needs no downloads — it is
the fast check that ingestion still works.
