# Our algorithm, and running it against the baselines

What our simplifier does, step by step; how `src/algo/simplify/convex_shadow.cpp` implements
it; and how `ssk_match` and `scripts/comparison_stats.py` run it through the protocol in
[comparison.md](comparison.md).

**Source:** Algorithm 1 (`Simplify`) of Siu-Wing Cheng, Haoqiang Huang and Le Jiang,
*Simplification of Trajectory Streams*, SoCG 2025,
[arXiv:2503.23025](https://arxiv.org/abs/2503.23025), Section 2. The polygon routines come
from the earlier prototype in `archive/streaming-simplification`.

**One deliberate deviation:** only the grid points of the *boundary* cells start a segment,
which gives O(1/ε) candidates instead of O(1/ε²). Section 3 shows that both of the paper's
guarantees survive this.

---

## 1. The problem and what is guaranteed

The input is a polyline `τ = v1, v2, …` arriving one vertex at a time, a tolerance `δ > 0`, and
`ε ∈ (0, 1)`. The output is a polyline `σ` with

- `d_F(σ, τ) ≤ (1 + ε)·δ`, the continuous Fréchet distance; and
- `|σ| ≤ 2·κ(τ, δ) − 2`, where `κ(τ, δ)` is the fewest vertices that *any* curve within δ of
  τ can have.

This is the paper's Theorem 3. The output vertices are grid points and points of shadows, not
input vertices, so `ConvexShadow` is a `Simplifier<2>` and not a `SubsetSimplifier<2>`.

```cpp
#include "simplify/convex_shadow.hpp"

ssk::simplify::ConvexShadow cs({{"delta", 0.5}, {"eps", 0.1}});   // eps defaults to 0.1
const ssk::simplify::Curve<2> out = cs.simplify(curve);
```

The constructor throws if `delta <= 0` or if `eps` is outside (0, 1). No timestamps are needed.

## 2. The algorithm, step by step

### Words

- **Stabbing in order.** A segment stabs regions `O1, …, Om` in order when it has points
  `x1, …, xm`, in that order along it, with each `xa` in `Oa`. A segment that stabs the balls
  `B(va, r)` of consecutive vertices in order is within Fréchet distance `r` of that piece of
  the input: match each input edge `va va+1` to the piece `xa xa+1`.
- **`G_v`.** Take the grid of side `h = εδ / (2√2)` that has `v` itself as a grid point, and keep
  the cells that meet the ball of radius `ρ = (1 + ε/2)·δ` around `v`. `conv(G_v)` contains
  `B(v, ρ)` and lies inside `B(v, (1+ε)δ)`, because `ρ` plus one cell diagonal (`εδ/2`) is
  `(1+ε)δ`.
- **Starts.** These are the grid points (cell corners) of the cells of `G_v` that meet the
  circle `|x − v| = ρ`. The paper uses every grid point of `G_v`.
- **Shadow `S[p]`.** For a start `p`, this is the set of points where a segment from `p` can
  end after stabbing `conv(G_v)` of every vertex of the current run in order.
- **`F(S, p) = { y : the segment p→y meets S }`.** Picture `p` as a light source: `F` is `S`
  plus the dark region behind it, an unbounded convex wedge.

### Steps

1. **Open a run** at a vertex `v`. Every start `p` gets the shadow `S[p] = {p}`.
2. **Feed the next vertex `w`.** For every start, set
   `S[p] ← conv(G_w) ∩ F(S[p], p)`. If `p ∈ S[p]`, then `F` is the whole plane and the new
   shadow is all of `conv(G_w)`. This always happens on a run's first step.
3. **If some shadow is non-empty**, the run continues. The segment it would emit right now is
   any start `p` with a non-empty shadow, together with any point `q` of that shadow (the
   paper's buffer `σ̃`).
4. **If every shadow is empty**, the run ended at the previous vertex. Emit the segment
   `p→q` from step 3 and **open a new run at `w`**, the vertex that emptied every shadow. The
   gap between the old `q` and the new start is an edge of the output, matched to the input
   edge between the two vertices.
5. **At the end of the stream**, emit the open run's segment. If the stream ended right after a
   run was opened, emit that vertex itself.

The output is `p1 q1 p2 q2 … pk qk`, possibly ending in a bare vertex. Only the current shadows
are kept, so working storage depends on ε and not on the length of the stream.

## 3. Why it is correct

**Lemma 1.** After the vertices `vs … vj` of a run,
`S[p] = { x ∈ conv(G_vj) : the segment p→x stabs conv(G_vs), …, conv(G_vj) in order }`.

- If `p→x` passes through some `s` in the old shadow, the earlier stabs lie on `[p, s]`, and `x`
  stabs the new region.
- If `p→x` stabs everything in order, its stab point for the previous vertex is in the old
  shadow and lies on the segment.

**The error bound.** Each emitted segment stabs its run's `conv(G_va)` in order, and every
point of `conv(G_v)` is within `(1+ε)δ` of `v`. Match each input edge to the piece of segment
between consecutive stab points. Match each joining edge `q → p'` to its input edge: `q` lies
in `conv(G_vi)` and `p'` is within `(1+ε)δ` of `vi+1`. Everything is then within `(1+ε)δ`.

**The size bound (Lemma 2).** Take an optimal curve `u1 … uκ` within δ, and let `xa` be the point
matched to `va`. If `b` is the last vertex matched onto `u1u2`, then `x1u2` stabs
`B(v1, δ), …, B(vb, δ)` in order.

- **Start on the grid.** Snap the start to a grid point. That moves every stab point by at most
  the snapping distance, so some start survives at least through `vb`.
- **Charge each run.** Each run therefore reaches at least as far as an optimal edge does, and
  each emitted segment can be charged to a distinct optimal edge. So there are at most `κ − 1`
  segments, and `|σ| ≤ 2κ − 2`.
- **Why the restart vertex matters.** This charging works only because a new run opens at the
  vertex that broke the old one, which is matched onto a later optimal edge. Opening it one
  vertex earlier, as the archive prototype did, voids the argument.

**Starting from the boundary only.** Extend `x1u2` *backwards*, past `x1`, until it leaves
`B(v1, ρ)` at a point `y` on the circle.

- **The longer segment still works.** `y→u2` contains `x1→u2`, so it stabs every ball in the
  same order, and `y` stabs the first ball itself.
- **There is a boundary start near `y`.** The cell containing `y` meets the circle, so it is a
  boundary cell, and its nearest corner is within half a cell diagonal, `εδ/4`, of `y`.
- **The snapped segment stays inside the grid regions.** Starting from that corner, the segment
  stabs the balls of radius `δ + εδ/4 ≤ ρ`, all of which lie inside `conv(G_v)`.
- **The start is close enough to the vertex.** The corner is within `ρ + εδ/2 = (1+ε)δ` of `v1`.

The charging argument then goes through unchanged.

## 4. How it is implemented

Everything is in `src/algo/simplify/convex_shadow.cpp`. Shadows are counter-clockwise convex
polygons stored as `std::vector<Vec>`.

### Once per run of `simplify`

- **`boundary_points(h, ρ)`** walks the cells `[i, i+1] × [j, j+1]` in grid units around the
  origin. For each cell it computes the distance to its nearest point (`near`) and to its
  farthest corner (`far`), and keeps the cell when `near ≤ ρ/h < far`, meaning it reaches the
  ball without lying wholly inside it. The four corners of every kept cell are deduplicated and
  scaled by `h`.
- **`hull(starts)` is `conv(G_v)`.** The boundary cells form a closed ring, so the hull of their
  corners contains everything the ring encloses, including every interior cell.
- **Both are offsets.** The grid moves with each vertex, so `G_v`, its hull and its starts are
  the same shapes translated, and they are built once and added to `v` when needed.

### `run()`

1. **Local coordinates.** Every vertex is shifted by the first one. The algorithm is
   translation-invariant, so this only changes rounding. It matters because δ can be a
   millionth of the coordinates themselves (Mopsi is degrees around 62, NGSIM is State Plane
   feet). The output is shifted back at the end.
2. **`start(v, starts)`** makes one candidate `{p = v + offset, shadow = {p}}` per start.
3. **`advance(cands, w, cells, r)`** builds `around = conv(G_w)` once and updates every
   candidate:
   - if `closest(shadow, p)` returns `p` itself, then `p` is in its shadow and the new shadow is
     `around`. `contains` counts a point as inside only when it is clearly inside (by `1e-9`,
     relative to the edge), because a point taken as inside can widen a shadow to all of
     `conv(G_w)` or become an output vertex. A start judged outside when it is really inside is
     harmless: the rays of `cover` then span the whole circle, and the result is `around`
     anyway;
   - otherwise the new shadow is `clip(around, cover(shadow, p, w, r))`;
   - a result with fewer than three vertices is dropped (see below).
4. **If no candidate survives,** `emit_segment` uses the candidate list from **before** this
   vertex, which is the paper's buffer. It emits the start nearest the run's first vertex, then
   the point of that start's shadow nearest the run's last vertex. A new run opens at `w`.
5. **At the end,** emit the open run, or the last vertex alone if it has only just opened a run.
   `emit` skips a point equal to the previous one.

### How `cover` builds `F` exactly where it matters

1. Measure the angle of every shadow vertex from `p`, relative to the direction of the
   shadow's centroid. The extremes are the two tangent directions.
2. Shoot rays from `p` across that range, at most 30° apart, each `1.3×` longer than the far
   side of `D(w, r)`.
3. Take the convex hull of the shadow and the ray tips.

Every ray tip lies in `F`, so the hull lies in `F`. The chords between consecutive tips stay
beyond `D(w, r) ⊇ conv(G_w)`, so after clipping to `around` the result is
`conv(G_w) ∩ F` exactly. `clip` cuts `around` by each edge of the cover in turn
(Sutherland–Hodgman).

### Zero-area shadows

A clip can leave a point or a segment where two polygons merely touch. The cover of such a
shadow would be a whole line, not a ray, which is too large. Dropping these shadows keeps every
approximation on the inside. It cannot cost the size bound, because Lemma 2's witness stabs
balls strictly inside `conv(G_v)`, so its shadow always has area.

### What differs from the paper

| | paper | here |
|---|---|---|
| starts | every grid point of `G_v` | grid points of boundary cells only (Section 3) |
| which `p`, `q` to emit | any | `p` nearest the run's first vertex, and `q` the point of `S[p]` nearest its last vertex |
| zero-area shadows | kept | dropped |
| dimension | any `d` | 2 |
| exact arithmetic | assumed | doubles; every test near a boundary errs inward, and `run()` works with δ less 16 ulps of the largest coordinate. Without that, points ours places exactly on `(1+ε)δ` could round past it when written back in absolute coordinates. A δ no larger than that margin returns the input unchanged |

### What differs from the archive prototype

| | archive `--shadow convex` | here |
|---|---|---|
| restart | at the previous vertex, which voids Lemma 2 | at the vertex that emptied every shadow |
| ball | inscribed m-gon of `B(v, (1+ε)δ)` | `conv(G_v)` |
| grid | anchored at the coordinate origin, starts at cell centres | anchored at each vertex, starts at grid points |
| starts | every cell centre within `(1+ε)δ`, about `8π(1+ε)²/ε²` | boundary cells, O(1/ε) |
| shadow size | capped at 64 vertices, which shrinks shadows | uncapped (Theorem 3 bounds it) |
| duplicate and collinear filter | yes | no |

### Starts per run

| ε | boundary starts | every cell centre (archive) |
|---|---:|---:|
| 0.05 | 920 | ≈ 11 084 |
| 0.1 | 472 | ≈ 3 041 |
| 0.3 | 168 | ≈ 472 |
| 0.5 | 120 | ≈ 226 |

## 5. The pipeline: `ssk_match`

`src/pipelines/error-match/`: the library is `match.{hpp,cpp}` (`ssk_match_lib`), beside a thin
`main` in `error_match.cpp`. It is step 3 of [comparison.md](comparison.md).

**It reads** each baseline operating point that `ssk_frechet` measured,
`data/frechet-distances/<algorithm>/<dataset>/m<rate>/<name>.json`, and takes from it the
measured distance `d*`.

**It runs ours at `delta = d* / (1 + ε)`.** Ours guarantees `(1+ε)·delta`, so it is within
the baseline's error by construction, and the comparison asks the protocol's question: at that
error, how many vertices? Running at `delta = d*` would let ours be up to `(1+ε)` times worse
than the baseline. This is the problem comparison.md raises under *Where ε goes*.

**It measures ours in the same pass**, with the same Fréchet computer and the same relative
tolerance (`1e-9` of the bounding-box diagonal) that `ssk_frechet` used on the baselines. The
protocol reports measured error on both sides, and measuring here saves reading 786 642
results back from disk.

**It writes**, under a root per ε:

```
data/error-matched/eps<E>/<algorithm>/<dataset>/m<rate>/<name>.json   ours, one per baseline run
data/error-matched/eps<E>/summary/<dataset>/<name>.json               every pair's numbers
```

Here `<algorithm>` is the baseline matched against. Each result is a trajectory document, so
`io::read_trajectory` and `viz/plot.py` read it directly. It holds the points (no `t`, since
ours invents its vertices), `algorithm: "convex-shadow"`, `mode: "error-matched"`,
`params {delta, eps}`, `stats {input_size, output_size, compression}`, its own
`frechet {distance, tolerance, bound, within_bound, computer}`, and under `matched` the
baseline document's `algorithm`, `params`, `stats`, `frechet` and `simplified`.

The **summary** holds one entry per pair, with the baseline's `output_size`, `target`,
`algorithm_runs`, `distance` and `tolerance` beside ours. It exists so that the statistics read
one small file per trajectory instead of every result.

| Flag | Default | What |
|---|---|---|
| `--in DIR` | required | a dataset directory of trajectory documents |
| `--eps E` | required | ours alone, so each value is its own run, and a result must say which it used |
| `--distances DIR` | `data/frechet-distances` | the baselines' measurements |
| `--out DIR` | `data/error-matched` | root of the result tree |
| `--tol REL` | `1e-9` | Fréchet bisection tolerance, relative to the input's diagonal, as in `ssk_frechet` |
| `--rates N` | 6 | deepest compression rate to look for |
| `--limit N` | all | stop after N trajectories |
| `--shard I --shards N` | 0 / 1 | process every Nth trajectory, starting at I |
| `--skip-existing` | off | skip a trajectory whose summary exists |

**It checks the guarantee as it goes.** Every run also takes the exact decision
`Frechet::within(input, ours, (1+ε)·delta)`, with no tolerance, which costs about a ninth of the
measurement. It is stored as `frechet.bound` and `frechet.within_bound`, the fields
`viz/plot.py` reports, and the driver prints a count of any run that fails it.

**It is resumable.** The summary is written last, so a trajectory that has one is finished.
After an interruption, rerun with `--skip-existing`.

### Parallel by construction

Every trajectory writes only its own paths, and every ε writes under its own root. Any number of
processes, whether shards of one dataset, different datasets or different ε, run side by side
with no queue, no lock and no shared file. Shards stride through the sorted file list, as in
the other drivers ([pipeline.md](pipeline.md) explains why, and the one way striding fails).

Unlike `ssk_frechet`, **every dataset is worth sharding here**. Ours costs 0.5–7 ms per input
point per run, depending on ε (Section 7), so even Mopsi is tens of core-hours. One ε over the whole
corpus, 8 processes at a time, skipping the `mopsi-sample` subset so it is not counted twice:

```powershell
foreach ($ds in Get-ChildItem data\trajectories -Directory -Exclude mopsi-sample) {
  $in = $ds.FullName
  0..7 | ForEach-Object -Parallel {
    .\build\src\pipelines\error-match\ssk_match.exe --in $using:in --eps 0.3 `
      --shard $_ --shards 8 --skip-existing
  } -ThrottleLimit 8
}
```

```sh
for ds in data/trajectories/*/; do
  [ "$(basename "$ds")" = mopsi-sample ] && continue
  seq 0 7 | xargs -P 8 -I{} ./build/src/pipelines/error-match/ssk_match \
    --in "$ds" --eps 0.3 --shard {} --shards 8 --skip-existing
done
```

`ForEach-Object -Parallel` needs PowerShell 7; on Windows PowerShell 5.1 use `Start-Job` in a
loop instead. Repeat the loop for each ε. A single dataset at a time is the same command
without the loop:

```powershell
.\build\src\pipelines\error-match\ssk_match.exe --in data\trajectories\mopsi --eps 0.3
```

```sh
./build/src/pipelines/error-match/ssk_match --in data/trajectories/mopsi --eps 0.3
```

## 6. The statistics: `scripts/comparison_stats.py`

The script uses the standard library only, so it needs no venv. It reads every
`eps*/summary/<dataset>/<name>.json`, spreading the reading over one process per core
(`--jobs N`), and prints the following.

**Per dataset:**

- **coverage**: trajectories in each (rate, baseline) cell, and our matches beside them. This
  is the per-cell count comparison.md asks for, which shows where DOTS is missing;
- **fixed rate**: each baseline's median measured error at `N/2^m`, which is the headline
  table, plus how often DOTS hit its budget exactly and how many runs its search took;
- **matched**, per ε: for each (baseline, m),
  - the median of our size over theirs;
  - the share of runs where ours is smaller, equal or larger;
  - the median of our measured error over theirs;
  - **Pareto dominance in both directions**: ours smaller and no less accurate, and theirs
    smaller and no less accurate.

**Pooled over all datasets:** the matched table again, ratios only, since errors are in each
dataset's own units.

**Theorem 3, checked**, per ε:

- **Error.** A pair violates the bound when its exact decision `within_bound` is false. The
  table also shows the worst `(distance − tolerance) / bound`, which is how far over the bound a
  run certainly is.
- **Size.** `κ(τ, delta)` is bounded from above by the input and by every curve on disk that is
  certainly within `delta`. Every stored distance is the top of a bisection, so each one is a
  certificate: each baseline counts, and so does ours at any ε. A pair is a certain violation
  when ours exceeds `2κ − 2` under that bound.
- **What a pass means.** κ is only bounded from above, so passing does not prove the size
  bound; failing does break it.

**Notes** carry the caveats comparison.md says must travel with the numbers: DOTS reaches its
budgets by a search, ours does not keep a subsequence, errors are per-dataset units, and SED is
not scored.

Everything above, plus a **rate–distortion curve** per algorithm (median compression and
median error per rate, with ours per ε and per matched baseline), is written as JSON to
`--out`, by default `<matched>/comparison.json`.

```powershell
python scripts\comparison_stats.py --matched data\error-matched
```

```sh
python scripts/comparison_stats.py --matched data/error-matched
```

## 7. Cost of the full comparison

**Measured.** `ssk_match` was run on a random sample of every dataset, about 40 000
point-runs each (a point-run is one input point simplified once), at three values of ε, with 7
processes at a time. The rates therefore already include the contention a parallel run sees.
None of the 2 778 sampled runs exceeded the bound.

| ε | ms per point-run, across the nine datasets |
|---|---|
| 0.1 | 3.6 – 6.8 |
| 0.3 | 0.86 – 1.58 |
| 0.5 | 0.51 – 0.85 |

**The full comparison is 844 377 088 point-runs per ε,** because each of the 786 642 baseline
operating points reruns its whole trajectory (17.9 runs per input point on average). Measuring
ours costs little in comparison: our outputs have 0.29–0.63 of the baselines' Fréchet cells,
so roughly half of the 6.6 core-hours `ssk_frechet` took, about 3.5 core-hours per ε.

| core-hours | ε = 0.1 | ε = 0.3 | ε = 0.5 |
|---|---:|---:|---:|
| `geolife` | 709 | 163 | 107 |
| `mopsi` | 221 | 42 | 34 |
| `ngsim-i-80` | 155 | 33 | 16 |
| `ngsim-us-101` | 136 | 33 | 16 |
| `ngsim-lankershim` | 46 | 9 | 5 |
| `mot-mot20` | 35 | 8 | 4 |
| `ngsim-peachtree` | 28 | 7 | 4 |
| `mot-dancetrack` | 15 | 3 | 2 |
| `mot-mot17` | 11 | 3 | 2 |
| **total** | **1 356** | **300** | **189** |
| **wall clock, 8 processes** | **≈ 7 days** | **≈ 37 h** | **≈ 24 h** |

All three values of ε come to about **1 850 core-hours, roughly 9.6 days on this 8-core
machine**, and ε = 0.1 is three quarters of that. Two ways to make it tractable:

- **Do ε = 0.5 and 0.3 first.** On the full corpus they take about 2.5 days together; read the
  statistics before committing to more.
- **Run ε = 0.1 on `mot-dancetrack` first, not the whole corpus.** That is 10.3 M point-runs,
  ≈ 15 core-hours, ≈ 2 h on 8 processes. Only run it on everything if that says it is worth a
  week. `mopsi-sample` is too small to decide: 40 trajectories, 4 115 points.

**Longest single job:** `geolife-001523`, 92 645 points with 18 runs, takes about 2.6 h on one
core at ε = 0.1 and 35 min at ε = 0.3. That is well under the total, but it is the tail at the
end of whichever GeoLife shard holds it.

**Disk:** 68 bytes per output point, and outputs are about half the baselines' size, so
roughly 4–5 GB per ε. That is an estimate from the sample.

**How far to trust this.**

- **Lengths sampled.** The sample drew trajectories of up to about 1 100 points, so that each
  dataset got several. Per-point cost does not depend on trajectory length, and the Fréchet part
  is scaled by cells, which do; but no long track was timed directly.
- **Cache.** I/O was measured with warm caches.
- **Calibration.** As with the other drivers, take the first dataset's measured rate as the
  calibration.

## 8. Verification so far

- **`ConvexShadow` stays within its bound.** On 12 Mopsi trajectories at ε ∈ {0.1, 0.3} × δ ∈
  {5e-5, 2e-4, 1e-3} degrees, 0 of 72 runs exceeded `(1+ε)δ` under the repo's exact
  `Frechet<2>::within`.
- **It is smaller and faster than the archive.** Against the archive prototype on the same
  inputs, it produced up to 28% fewer points and ran 2.5–4.4× faster.
- **`ssk_match` and the stats script run end to end.** On 5 Mopsi trajectories at ε = 0.3 and
  0.5, there were 0 violations of either half of Theorem 3, the resume path skipped every
  finished trajectory, and both print as shown above. The strict decision at `(1+ε)δ` found 0 violations
  in 2 006 runs: the calibration samples of all nine datasets at ε = 0.5 and 0.3, plus 3
  trajectories each of `mot-dancetrack`, `mopsi` and `ngsim-us-101` at ε = 0.5. `viz/plot.py`
  draws a result document and reports the check.
- **Not done:** there are no unit tests for `ConvexShadow` or `ssk_match` in `tests/`.

## 9. Open concerns

1. **Measured error sits near the bound.** Every run starts on the boundary ring, at least about
   δ from its first vertex, so our measured error lands close to `(1+ε)·delta` even where a
   closer start existed. In the trial above, the median per cell was 0.87–0.995 of the
   baseline's error. The bound and the vertex count are unaffected, and the comparison is still
   fair, because ours never exceeds the baseline's error. A start could be slid along its own segment toward the vertex
   after the fact; the paper does not do this.
2. **Units.** `delta` is in coordinate units, which are degrees for Mopsi and GeoLife. This is
   consistent within the protocol, since every error is measured the same way, but it is not
   metres.
3. **SED is not scored.** Ours has no timestamps, so it cannot be scored under SED alongside
   SQUISH and DOTS.
4. **κ is only bounded from above,** so the size check can refute the bound but never prove it.
5. **Volume.** Each ε writes 786 642 result documents plus 47 916 summaries, roughly 4–5 GB
   (Section 7).
6. **Batch interface only.** `ConvexShadow::simplify` takes the whole curve. The algorithm is
   streaming, but there is no `push`/`finish` yet, so latency and per-stream memory cannot be
   measured.
