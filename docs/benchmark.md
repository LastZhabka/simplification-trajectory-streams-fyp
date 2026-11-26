# Benchmark: the four implementations of ours

How fast each greedy stabber is, on the same inputs, as ε shrinks. All four are Algorithm 1
of *Faster δ-simplification of Trajectory Streams*, or Algorithm 2 which computes the same
thing faster; [convex-shadow-32.md](convex-shadow-32.md) §10–11 describes them:

| `--ours` | short | region | starts | per vertex |
|---|---|---|---|---|
| `greedy-stabbing-shared` | shared | m-gon | M = ⌈2π/ε⌉ on the δ-circle | O(M + m²) = O(ε⁻¹), Algorithm 2 |
| `greedy-stabbing-naive` | naive | m-gon | the same | O(Mm) = O(ε^(−3/2)), Algorithm 1 |
| `greedy-stabbing-naive-paper-boundary` | paper-boundary | `conv(G_v)` | grid points of the cells of `G_v` meeting the circle, O(ε⁻¹) | O(ε^(−5/3)) |
| `greedy-stabbing-naive-paper` | paper | `conv(G_v)` | every grid point of `G_v`, O(ε⁻²) | O(ε^(−8/3)) |

Algorithm 1 costs O(sides of the region) per start and vertex. `conv(G_v)` has Θ(ε^(−2/3)) sides,
the hull of the lattice points in a disc ([convex-shadow-32.md](convex-shadow-32.md) §4), where the
m-gon has Θ(ε^(−1/2)); hence the paper variants' bounds.

The question is the one Algorithm 2 exists to answer: **what does sharing the offsets and the
`Q[p]` between starts buy, and what does the paper's grid cost**, at the ε a user would pick.

## Running it

Steps 1–6 of [running.md](running.md) first: the selection reads `data/trajectories`, and the
δ values are read from `data/frechet-distances`.

```powershell
cmake --build build -j
python scripts\benchmark.py [--eps 0.5,0.25,0.2,0.1,0.05] [--time-limit 300] [--jobs 8] `
                            [--per-dataset 10] [--out data\benchmark\benchmark.json]
```

```sh
cmake --build build -j
python scripts/benchmark.py [--eps 0.5,0.25,0.2,0.1,0.05] [--time-limit 300] [--jobs 8] \
                            [--per-dataset 10] [--out data/benchmark/benchmark.json]
```

About 25 minutes with the defaults on the machine below. `benchmark.json` is rewritten every
1 000 cells, so a stopped run leaves the cells it finished, fewer than its `cells_total`, and
`--resume` goes on from there; the rest is written at the end. Beside `--out` it writes:

| file | what |
|---|---|
| `benchmark.json` | every cell, the selection and the settings |
| `times/<ours>.json` | every time of one implementation, grouped by trajectory with its N: each δ and ε, its status, the best and every repeat, the output size, and `ns_per_vertex`, the best over N; a cell that hit the limit has `ns_per_vertex_at_least`, the limit over N |
| `summary/<ours>.json` | that implementation's distributions per ε, over all cells and per dataset: seconds, ns per input vertex and output size, as count, min, p10, p25, median, p75, p90, p99, max and mean; and how many finished, hit the limit or failed |
| `pairwise.json` | for each pair a, b, per ε, over all cells and per dataset: the distribution of b's time over a's, how many times faster a is. `finished` counts only cells both finished; `with_bounds` also counts those one of them did not, with the limit standing in for its time |

and prints the tables below as Markdown; redirect it to keep them (`data/benchmark/report.md`
for the run below).

**Running only some of them.** `--ours` names the implementations to run, comma-separated, and
`--reuse` a finished `benchmark.json` to take every other one's cells from. The script refuses one
whose selection, ε or time limit differs, so the same trajectories, δ and ε are timed. Its
`benchmark.json` and derived files then cover all seven, the reused cells as they were, and
`settings.reused` names what came from where. **A run stopped part way** goes on from its last
checkpoint with the same command and `--resume`: the cells of the implementations run that
`--out` already holds are kept, after the same checks as `--reuse`, and only the others run.
`settings.resumed` counts the kept cells, and `wall_seconds` adds up the runs. The three
no-decomp stabbers against the run at 1 007 trajectories, which was resumed twice, both times
from the checkpoint at 10 000 cells:

```powershell
python scripts\benchmark.py --per-dataset 111 --time-limit 200 `
    --ours greedy-stabbing-no-decomp,greedy-stabbing-no-decomp-paper-boundary,greedy-stabbing-no-decomp-paper `
    --reuse data\benchmark\per-dataset-111\benchmark.json --out data\benchmark\no-decomp-111\benchmark.json
```

```sh
python scripts/benchmark.py --per-dataset 111 --time-limit 200 \
    --ours greedy-stabbing-no-decomp,greedy-stabbing-no-decomp-paper-boundary,greedy-stabbing-no-decomp-paper \
    --reuse data/benchmark/per-dataset-111/benchmark.json --out data/benchmark/no-decomp-111/benchmark.json
```

## What is measured

A **cell** is one trajectory, one δ, one ε and one implementation. Each is one `ssk_bench`
process (`src/pipelines/benchmark/`), which reads the trajectory, builds the implementation
and times `simplify()` alone with `steady_clock` — not the file read. A run under a second is
repeated, up to 5 times or until the repeats add up to 1 s, and the **best** repeat is
reported: it is the one least disturbed by the rest of the machine. Each cell also reports the
output's vertex count and a 64-bit FNV-1a digest of its coordinates, so outputs meant to agree
to the bit can be compared.

**The time limit** is per cell, 300 s by default, on the whole process. A cell that has not
printed by then is killed and recorded as `time_limit_exceeded`, printed as **TLE**: no time,
no vertex count. The limit is enforced from outside because a thread cannot be killed, and a
deadline check inside the stabbers would put the benchmark into the algorithm. Reading the
file is milliseconds even for the largest trajectory, so the limit is in effect the run's.
At 300 s only the paper stabber on the two largest trajectories at ε = 0.05 hits it; a first
run at 60 s cut off 17 cells.

**Parallelism.** Cells run 8 at a time. The machine is an Intel Core Ultra 9 275HX, 8
performance cores and 16 efficiency cores, no hyperthreading, so concurrency beyond 8 puts
runs on the slower cores. Measured on one naive run at ε = 0.05 run concurrently with copies
of itself: 1 job 1.412 s; 8 jobs 1.406–1.495 s (6% spread); 16 jobs 14%; 24 jobs 20%. At 8,
timings are within 6% of a quiet machine. The cells are shuffled with a fixed seed, so no
implementation always runs beside the slow ones.

## What is selected, and why

**The cost is per vertex and set by ε, not by δ.** A calibration pass before choosing
(mopsi-006752, geolife-001527, ngsim-002303; δ at m2 and m6) found each implementation's time
within a few percent between m2 and m6 at the same ε, and linear in N. So the selection is
about which traces and how long, and is made by rule in `scripts/benchmark.py`, never by hand:

- **A sample of every dataset**: 10 trajectories of each of the nine, drawn at random among
  those with at least 500 points, seeded by the dataset's name. Below 500 points the shared
  stabber takes microseconds, which measures the timer, not the algorithm. Every dataset gets
  the same number, so the four motion types — GPS walking, cycling and driving (GeoLife, Mopsi),
  10 Hz freeway and arterial driving (NGSIM, four sites), pixel-space pedestrian and dancer
  tracks (MOT17, MOT20, DanceTrack) — all weigh in, rather than GeoLife's 18 670 drowning the
  rest. A trajectory whose three δ are not all positive is passed over.
- **The long tail** beside it (`TAIL`): the largest trajectory and the one at the 99th
  percentile of size of GeoLife and Mopsi, and the largest of NGSIM us-101 and lankershim, MOT20
  and DanceTrack. That is where a super-linear cost in ε shows as a time limit:

  | dataset | rule | trajectory | N |
  |---|---|---|---:|
  | geolife | largest | geolife-001523 | 92 645 |
  | geolife | p99 | geolife-001527 | 11 611 |
  | mopsi | largest | mopsi-006752 | 31 741 |
  | mopsi | p99 | mopsi-000446 | 7 871 |
  | ngsim-us-101 | largest | ngsim-002303 | 3 719 |
  | ngsim-lankershim | largest | ngsim-000414 | 2 183 |
  | mot-mot20 | largest | mot-002319 | 3 315 |
  | mot-dancetrack | largest | mot-000104 | 2 402 |

That is **98 trajectories, 265 218 points**. `--per-dataset` changes the sample; at 111 it is
1 007 trajectories and 1.4 M points, about 15 CPU-hours, nearly all of it the paper stabber.

**δ is where the baselines operate.** For each trajectory, δ is the Fréchet distance
Douglas-Peucker-n reached with N/4, N/16 and N/64 vertices (`m2`, `m4`, `m6` in
`data/frechet-distances/douglas-peucker-n/`). That is the range of step 7 of
[running.md](running.md), in each dataset's own units, so no δ is picked by hand and none
means nothing on its trajectory. Three of them, low to high compression, since δ changes the
output and should not change the time.

**ε** is 0.5, 0.25, 0.2, 0.1, 0.05: from coarse to the regime where the gap between O(ε⁻¹)
and O(ε^(−3/2)) and beyond is widest.

That is 98 trajectories × 3 δ × 5 ε × 4 implementations = **5 880 cells**.

## Measured

2026-11-23, the defaults, at commit f7af71e plus the benchmark: 5 880 cells in 1 523 s of wall
clock. **5 875 finished, 5 hit the 300 s limit, 0 errors.** The shared and the naive stabber's
outputs are **identical to the bit in all 1 470 cells** where both finished, as §11 of
[convex-shadow-32.md](convex-shadow-32.md) says they must be.

This run predates the fix to `geometry::closest` of 2026-11-24 (`docs/running.md` step 7). The
fix changes only which point of a shadow is emitted, once per output segment; on these 1 470
cases it changes 3 outputs' coordinates and no vertex count, so the timings stand.

It also predates the shared stabber's exact marking of 2026-11-24
([convex-shadow-32.md](convex-shadow-32.md) §11, step 2), which makes shared faster and leaves
the other three as they were. So shared's columns below are its old times, and its speedups
are 1.1–1.2× low.

#### Shared, before and after the exact marking

The same 1 470 cells, shared only, both builds interleaved in one shuffled queue at 8 jobs.
Each cell gets its best of 3 rounds, each round being `ssk_bench`'s own best of up to 5:

| ε | before | after | summed time | median cell | ns / vertex |
|---|---:|---:|---:|---:|---:|
| 0.5 | 0.854 s | 0.775 s | −9.2% | −6.5% | 1 073 → 975 |
| 0.25 | 1.150 s | 1.014 s | −11.8% | −7.3% | 1 445 → 1 275 |
| 0.2 | 1.342 s | 1.159 s | −13.6% | −9.7% | 1 686 → 1 456 |
| 0.1 | 2.054 s | 1.741 s | −15.2% | −11.7% | 2 581 → 2 188 |
| 0.05 | 3.132 s | 2.608 s | −16.7% | −13.1% | 3 936 → 3 277 |

**Outputs identical to the bit in all 1 470.** A throwaway build that, at every vertex, also
computed every start's sides by brute force counted the following per vertex:

| ε | starts alive | changes | marks before | marks after | builds before | builds after | needed |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 11.1 | 1.33 | 1.84 | 1.33 | 2.43 | 2.02 | 1.99 |
| 0.25 | 21.2 | 2.06 | 2.93 | 2.06 | 3.13 | 2.51 | 2.47 |
| 0.2 | 25.9 | 2.40 | 3.40 | 2.40 | 3.43 | 2.71 | 2.68 |
| 0.1 | 49.8 | 3.66 | 5.23 | 3.66 | 4.70 | 3.59 | 3.54 |
| 0.05 | 98.3 | 5.35 | 7.76 | 5.35 | 6.39 | 4.71 | 4.66 |

How to read it:

- A *change* is a start whose sides differ from the start before's.
- *Needed* counts the builds where the survivor's sides differ from those of the last build,
  plus the first build of each vertex.
- After the change the marks are exactly the changes. Before it, 28–31% of the marks were
  spurious, and 17–26% of the builds were wasted.
- Neither marking ever missed a change.
- What still exceeds *needed*, 0.03–0.05 per vertex, is a side that changes and changes back
  between two survivors.

### Time, summed per ε

Best-of seconds summed over the 294 (trajectory, δ) cells:

| ε | shared | naive | paper-boundary | paper |
|---|---:|---:|---:|---:|
| 0.5 | 0.9 s | 3.5 s | 38.7 s | 74.4 s |
| 0.25 | 1.3 s | 6.6 s | 79.5 s | 252.5 s |
| 0.2 | 1.4 s | 8.6 s | 97.1 s | 369.1 s |
| 0.1 | 2.3 s | 19.1 s | 213.7 s | 1 429.1 s |
| 0.05 | 3.5 s | 41.3 s | 560.6 s | 4 301.8 s over 289, **5 TLE** |

The shared stabber does all 1 470 of its cells in 9.4 s together; the paper stabber needs more
than 300 s for some single cells.

### Pairwise: how many times faster

Median over the 294 cells of each ε of the column's time over the row's; a cell where one hit
the limit counts at the limit. The full distributions, per dataset too, are in
`pairwise.json`.

| ε | shared vs naive | shared vs paper-boundary | shared vs paper | naive vs paper-boundary | naive vs paper | paper-boundary vs paper |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 | 4.3× | 45× | 88× | 10.8× | 21× | 1.9× |
| 0.25 | 6.2× | 73× | 234× | 11.7× | 38× | 3.2× |
| 0.2 | 7.9× | 88× | 341× | 11.2× | 43× | 3.8× |
| 0.1 | 10.6× | 112× | 831× | 10.8× | 77× | 7.0× |
| 0.05 | 15.6× | 205× | ≥ 2 992× | 13.3× | ≥ 186× | ≥ 13.9× |

The smallest speedup of shared over naive in any cell is 1.6× at ε = 0.5 and 3.0× at 0.05;
over paper-boundary 12× and 45×. **Shared is the fastest in every one of the 1 470 cells.**

Down a column it follows the bounds: shared over naive grows about as ε^(−1/2), 3.6× for ε ten
times finer, the m it saves per start; anything over paper grows fastest, since paper's starts
are O(ε⁻²). Naive over paper-boundary stays at 11–13×, both having O(ε⁻¹) starts: paper-boundary
has 7–9× the starts (472 against 63 at ε = 0.1) and a region of 2–2.5× the sides (28 against
11).

### Per dataset

Median speedup of shared over naive:

| dataset | trajectories | ε 0.5 | 0.25 | 0.2 | 0.1 | 0.05 |
|---|---:|---:|---:|---:|---:|---:|
| geolife | 12 | 3.3× | 4.9× | 4.9× | 8.0× | 10.8× |
| mopsi | 12 | 3.5× | 4.9× | 5.9× | 8.7× | 11.8× |
| mot-dancetrack | 11 | 2.3× | 3.0× | 3.0× | 4.3× | 6.0× |
| mot-mot17 | 10 | 4.5× | 9.0× | 10.9× | 16.0× | 25.2× |
| mot-mot20 | 11 | 3.9× | 5.7× | 6.7× | 9.3× | 12.1× |
| ngsim-i-80 | 10 | 6.0× | 10.3× | 11.8× | 19.3× | 28.8× |
| ngsim-lankershim | 11 | 6.4× | 7.4× | 11.1× | 16.8× | 21.8× |
| ngsim-peachtree | 10 | 5.6× | 7.2× | 10.5× | 15.8× | 19.4× |
| ngsim-us-101 | 11 | 6.2× | 9.5× | 11.0× | 18.4× | 29.2× |

It grows with 1/ε on every dataset. It is largest on the smooth NGSIM vehicle tracks and
smallest on DanceTrack; why the trace matters this much has not been looked into. Against
paper-boundary and paper the same order holds, at 26–69× and 48–133× for ε = 0.5 and
77–396× and 1 164–5 518× for ε = 0.05 (`report.md`).

### Time per input vertex

Nanoseconds, summed time over summed N:

| ε | shared | naive | paper-boundary | paper |
|---|---:|---:|---:|---:|
| 0.5 | 1 183 | 4 378 | 48 598 | 93 474 |
| 0.25 | 1 594 | 8 337 | 99 866 | 317 348 |
| 0.2 | 1 797 | 10 846 | 121 978 | 463 925 |
| 0.1 | 2 831 | 23 989 | 268 645 | 1 796 099 |
| 0.05 | 4 419 | 51 945 | 704 593 | 9 470 277 * |

\* over the 289 cells that finished, so it understates.

From ε = 0.5 to 0.05, ten times finer, shared is 3.7× slower, naive 12×, paper-boundary 14×
and paper 101× at least. All grow more slowly than their worst cases — 10×, 32×, 46× and 464× —
which assume every start survives every vertex; on these traces they do not. The order and the
spread between them are the bounds'.

### What the paper's grid buys

Vertices, summed over all 294 cells of each ε (paper's are paper-boundary's within a vertex or
two wherever both finished):

| ε | shared | paper-boundary | shared's extra |
|---|---:|---:|---:|
| 0.5 | 39 555 | 37 715 | +4.9% |
| 0.25 | 42 917 | 41 848 | +2.6% |
| 0.2 | 43 636 | 42 718 | +2.1% |
| 0.1 | 45 140 | 44 702 | +1.0% |
| 0.05 | 46 061 | 45 811 | +0.5% |

**The grid's denser starts save 0.5–4.9% of the vertices, for 45–205× the time.** The saving
shrinks with ε while the cost grows. At a fixed δ a smaller ε tightens the guarantee to
(1 + ε)δ, so every implementation takes *more* vertices as ε shrinks: a smaller ε buys a
tighter error, not a smaller output. For that tighter error, **shared at ε = 0.05 is 11×
faster than paper-boundary at ε = 0.5** — 3.5 s against 38.7 s on the same 294 cells — with a
bound of 1.05δ instead of 1.5δ, at 22% more vertices.

### Where the limit was hit

All five are the paper stabber at ε = 0.05: geolife-001523 (92 645 points) at all three δ, and
mopsi-006752 (31 741 points) at m2 and m6 — at m4 it finished inside 300 s. Shared and naive
never came near the limit.

### Every cell

`report.md` ends with every cell of the eight long-tail trajectories. Every cell of the run is
in `benchmark.json` and, by implementation and trajectory, in `times/<ours>.json`.

## Measured at 1 007 trajectories

2026-11-24, `--per-dataset 111 --time-limit 200`, at commit 584ba5c (the exact marking), into
`data/benchmark/per-dataset-111/`. The selection is 111 per dataset plus the long tail:
**1 007 trajectories, 1 414 193 points, 60 420 cells**. It took 9 278 s of wall clock, 2 h 35 min
at 8 jobs. **60 414 finished, 6 hit the 200 s limit, 0 errors.** All 6 are the paper stabber at
ε = 0.05, on geolife-001523 (92 645 points) and mopsi-006752 (31 741), at all three δ each.
Shared and naive are **identical to the bit in all 15 105 cells**.

Best-of seconds summed over the 3 021 (trajectory, δ) cells of each ε:

| ε | shared | naive | paper-boundary | paper |
|---|---:|---:|---:|---:|
| 0.5 | 4.5 s | 19.1 s | 207.3 s | 402.2 s |
| 0.25 | 5.9 s | 36.7 s | 430.5 s | 1 414.3 s |
| 0.2 | 7.6 s | 47.1 s | 535.5 s | 2 056.0 s |
| 0.1 | 10.1 s | 105.0 s | 1 161.9 s | 8 039.4 s |
| 0.05 | 14.8 s | 222.7 s | 3 026.6 s | 37 677.6 s over 3 015, **6 TLE** |

Median speedup of shared over cells, and the smallest in any cell. A cell that hit the limit
counts at the limit:

| ε | vs naive | vs paper-boundary | vs paper |
|---|---:|---:|---:|
| 0.5 | 4.7× (min 1.2×) | 50× (min 7.5×) | 97× (min 19×) |
| 0.25 | 7.0× (min 1.4×) | 80× (min 13×) | 261× (min 29×) |
| 0.2 | 7.2× (min 1.5×) | 81× (min 10×) | 313× (min 26×) |
| 0.1 | 12.6× (min 2.7×) | 137× (min 17×) | 974× (min 88×) |
| 0.05 | 18.1× (min 3.6×) | 241× (min 26×) | ≥ 3 368× (min 357×) |

Nanoseconds per input vertex, the distribution over cells as p10 / median / p90:

| ε | shared | naive | paper-boundary | paper |
|---|---:|---:|---:|---:|
| 0.5 | 658 / 947 / 1 669 | 3 831 / 4 507 / 5 468 | 40 k / 49 k / 58 k | 79 k / 94 k / 111 k |
| 0.25 | 794 / 1 215 / 2 353 | 7 320 / 8 682 / 10 276 | 83 k / 101 k / 121 k | 274 k / 332 k / 393 k |
| 0.2 | 935 / 1 476 / 3 210 | 9 194 / 11 141 / 13 320 | 103 k / 125 k / 150 k | 399 k / 485 k / 571 k |
| 0.1 | 1 159 / 1 920 / 4 375 | 20 k / 25 k / 30 k | 218 k / 271 k / 328 k | 1.60 M / 1.92 M / 2.28 M |
| 0.05 | 1 598 / 2 791 / 6 660 | 43 k / 53 k / 63 k | 569 k / 713 k / 868 k | 8.2 M / 9.8 M / 11.8 M |

How it compares with the 98-trajectory run above:

- Shared's speedups are 0.9–1.2× those of the 98-trajectory run. They are higher at every ε
  but 0.2, by about what the exact marking gained (1.1–1.2×), though the selection differs
  too.
- Per vertex, shared has the widest spread of the four, 2.5–4.2× from p10 to p90. Its cost per
  vertex depends on how many starts survive, which the trace sets. The naive stabbers spread
  only 1.4–1.5×.
- Shared's extra vertices over paper-boundary are the same: +4.8% at ε = 0.5, +0.5% at 0.05.
- Per dataset, shared over naive is again largest on NGSIM (30× at ε = 0.05 on I-80 and
  US-101) and smallest on DanceTrack (7×).

Every table, per dataset too, is in `report.md` beside the data. Every cell, with
`ns_per_vertex`, is in `times/<ours>.json`.

Why shared spreads so much more is the trace. At ε = 0.05 its median per dataset runs from
1 720 ns on NGSIM I-80 and 1 728 on US-101, the peak of its histogram, up to 6 852 on
DanceTrack, the tail. Naive's runs only from 48 823 (DanceTrack) to 61 252 (MOT17).

δ adds to it, for shared only. Its median at ε = 0.05 is 2 364 ns at m2, 2 956 at m4 and
3 624 at m6 (884, 970, 1 021 at ε = 0.5); the other three stay within 3% across the three δ.
So *What is selected* above, that the time is set by ε and not by δ, holds for them and no
longer for shared.

#### What sets shared's time: the rebuilds

The instrumented build of *Shared, before and after the exact marking* was run once more, on
every (trajectory, δ) of this run at ε = 0.5 and 0.05: 6 042 runs, 75 s. It counts, per vertex,
the starts alive and the rebuilds of `Q[p]`. The counts are in `shared_counts.json` beside the
data. Medians over cells, and a least-squares fit of each cell's ns per vertex to its counts:

| | ε = 0.5 | ε = 0.05 |
|---|---:|---:|
| starts alive per vertex | 11.2 of 13 | 99.9 of 126 |
| rebuilds per vertex (bound 2m) | 1.47 (12) | 2.79 (30) |
| naive ≈ a + b·alive | 218 + 391·alive, R² 0.55 | 1 768 + 509·alive, R² 0.84 |
| shared ≈ a + b·alive + c·rebuilds | 307 + 0.4·alive + 418·rebuilds, R² 0.91 | 853 + 3.3·alive + 569·rebuilds, R² 0.95 |
| shared on alive alone | R² 0.01 | R² 0.01 |

Shared's time is its rebuilds, each about what naive pays per live start, since it is the same
work. A live start costs shared 0.4–3 ns.

- **Per dataset**, at ε = 0.05, the rebuilds per vertex run from 1.16 (US-101) and 1.20
  (I-80) to 9.63 (DanceTrack), in the order of shared's time. The live starts run only from 93.9
  to 125.0, in the order of naive's. Paper-boundary and paper are 13.4–13.7× and 183–192× naive
  on every dataset.
- **Per δ**, the rebuilds grow from 1.88 to 3.18 to 4.39 at m2, m4, m6 (1.21, 1.59, 1.85 at
  ε = 0.5). The live starts stay at 98–101.
- **Segment length does not explain it.** MOT17 has the longest segments, a median of 177 input
  vertices per output segment at ε = 0.05, and among the fewest rebuilds. The rank correlation
  between segment length and shared's time over all cells is 0.15.

Why the trace matters: `ℓ_k` meets the δ-circle around `T_1` iff
`−δ − b < min_j n_k·(T_j − T_1) < δ − b`, with `(1 + ε/2)δ ≤ b ≤ (1 + ε)δ`. So a line can split
the starts only while the run has gone more than `εδ/2`, and nowhere more than `(2 + ε)δ`, in
the direction `−n_k`.

- A straight run soon leaves that range in every direction but those near perpendicular to it.
- A run that turns within a few δ, like a dancer, stays in it in many directions.

That fits the per-dataset order. Why a larger δ adds rebuilds has not been measured.

### Histograms

`scripts/plot_benchmark.py` draws time per input vertex from `times/`, into `plots/` beside
it. It draws all seven implementations, so it takes the `times/` of the direct run below, which
holds this run's four as well:

- One histogram per implementation and ε (`histograms/<ours>_eps<E>.png`, 35).
- `01_overview.png`, all of them as a grid.
- `02_by_eps.png`, each implementation's ε on one axis.

The histograms have a logarithmic x axis and a dashed median. The 57 cells that hit the limit
have only a lower bound, so they are counted in their panel's note rather than binned.

The figures name the implementations as the write-up does: Algorithm 2 (shared), Algorithm 1
(naive), Paper-Boundary and Paper-Grid (paper), and the direct ones Direct-m-gon, Direct-Boundary
and Direct-Grid (no-decomp, no-decomp-paper-boundary, no-decomp-paper). Three line charts follow,
the median against ε from 0.5 on the left to 0.05 on the right, one line per implementation:

- `03_eps_lines_log.png`, log y.
- `04_eps_lines_linear.png`, linear y.
- `05_eps_band.png`, log y, shaded from p10 to p90.

Two more are for the write-up, as PDF and PNG: white, in Computer Modern, no title or footer,
drawn 3.15 in wide to sit two in a row of a 6.5 in text block at their own size:

- `06_writeup_eps`, 05 with the four decomposed lines named on it rather than in a legend, and a
  marker per region and starts for grey print. A direct line is dashed in its twin's colour,
  with its twin's marker hollow and no band, and a key tells dashed from solid. A name sits
  under its line, or over the line and its dashed twin where the next line down is closer than
  half a decade, never between the two.
- `07_writeup_distributions`, 02 as one panel: a row per implementation on one log axis, each
  ε a kernel density over log time, light to dark as ε shrinks, each row scaled to the same
  height. A direct row is dashed, over its twin's.

It needs matplotlib, from the venv.

Their legend gives k, the exponent in time ∝ ε^k fitted over the five ε. Each is set beside
its worst case:

| | shared | naive | paper-boundary | paper | no-decomp | no-decomp-paper-boundary | no-decomp-paper |
|---|---:|---:|---:|---:|---:|---:|---:|
| measured k | −0.47 | −1.08 | −1.16 | −2.01 | −1.13 | −1.44 | −2.30 |
| worst case over these five ε, from the counts | −0.88 | −1.40 | −1.43 | −2.24 | | | |
| worst case | −1 | −3/2 | −5/3 | −8/3 | −3/2, times log | −5/3, times log | −8/3, times log |

Every one is shallower than its bound, as *Time per input vertex* above says: not every
start survives every vertex. A direct stabber's bound is its twin's times the log of the
hull's sort, while its shadow has O(sides of the region) sides. To draw the figures:

```powershell
.\.venv\Scripts\python.exe scripts\plot_benchmark.py --times data\benchmark\no-decomp-111\times
```

```sh
.venv/bin/python scripts/plot_benchmark.py --times data/benchmark/no-decomp-111/times
```

## The direct stabbers at 1 007 trajectories

2026-11-26, the three no-decomp stabbers of [convex-shadow-32.md](convex-shadow-32.md) §12, the
paper's own step by a hull in the dual, at commit ff7019a. Run as `--per-dataset 111
--time-limit 200 --ours greedy-stabbing-no-decomp,greedy-stabbing-no-decomp-paper-boundary,greedy-stabbing-no-decomp-paper
--reuse data/benchmark/per-dataset-111/benchmark.json`, into `data/benchmark/no-decomp-111/`.
That is the same 1 007 trajectories, δ, ε and limit, and the same machine, with the four
decomposed stabbers' cells taken from the run above. Its `times/`, `summary/`, `pairwise.json`
and `report.md` cover all seven.

**45 315 cells, 45 264 finished, 51 hit the 200 s limit, 0 errors.** The run took 6 h 19 min of
wall clock at 8 jobs, summed to its checkpoints. It was stopped twice, once when memory ran low
and once when the laptop restarted, and went on each time with `--resume`. The cells done
since the last checkpoint were run again.

The limit was hit only by the grid variants. no-decomp-paper hit it 3 times at ε = 0.1 and
45 times at ε = 0.05; no-decomp-paper-boundary hit it 3 times at ε = 0.05. Paper hit it 6
times at ε = 0.05 in the run above.

Best-of seconds summed over the 3 021 cells of each ε:

| ε | no-decomp | no-decomp-paper-boundary | no-decomp-paper |
|---|---:|---:|---:|
| 0.5 | 30.2 s | 355.6 s | 677.4 s |
| 0.25 | 59.2 s | 857.4 s | 2 779.0 s |
| 0.2 | 76.8 s | 1 067.3 s | 4 050.4 s |
| 0.1 | 174.4 s | 2 939.5 s | 19 359.8 s over 3 018, **3 TLE** |
| 0.05 | 409.0 s | 9 149.8 s over 3 018, **3 TLE** | 112 551.8 s over 2 976, **45 TLE** |

Nanoseconds per input vertex, p10 / median / p90 over the cells that finished. Beside each,
the median over cells of its time over its naive twin's, counting a cell that hit the limit at
the limit:

| ε | no-decomp | over naive | no-decomp-paper-boundary | over paper-boundary | no-decomp-paper | over paper |
|---|---:|---:|---:|---:|---:|---:|
| 0.5 | 5 485 / 7 215 / 8 646 | 1.63× | 60 k / 83 k / 103 k | 1.76× | 113 k / 159 k / 196 k | 1.73× |
| 0.25 | 10 560 / 14 067 / 17 180 | 1.67× | 138 k / 200 k / 250 k | 2.04× | 443 k / 643 k / 807 k | 1.99× |
| 0.2 | 13 491 / 18 090 / 22 200 | 1.68× | 174 k / 247 k / 308 k | 2.04× | 662 k / 946 k / 1.18 M | 2.01× |
| 0.1 | 29 k / 41 k / 51 k | 1.71× | 449 k / 680 k / 857 k | 2.58× | 3.28 M / 4.85 M / 6.33 M | 2.56× |
| 0.05 | 66 k / 96 k / 118 k | 1.89× | 1.54 M / 2.32 M / 2.93 M | 3.36× | 21.8 M / 31.9 M / 42.2 M | 3.29× |

- **The decomposition alone is worth 1.6–3.4×.** It is the only difference between each pair.
  For the m-gon the gain grows little with finer ε, 1.63× to 1.89×. For `conv(G_v)` it grows
  from 1.7× to 3.3×.
- **The direct stabbers grow faster as ε shrinks.** Their k is −1.13, −1.44 and −2.30 against
  their twins' −1.08, −1.16 and −2.01. From ε = 0.5 to 0.05 the median grows 13×, 28× and 201×
  against 12×, 15× and 104×. The hull's sort is a log factor. What else grows with
  `conv(G_v)`'s side count has not been measured.
- **Shared over no-decomp-paper**, the paper's algorithm as written, is a median 154× at ε = 0.5,
  478× at 0.25, 561× at 0.2, 2 268× at 0.1 and at least 10 062× at 0.05.

**The same outputs, up to rounding**, over the cells where both finished:

| | cells | same vertex count | identical to the bit |
|---|---:|---:|---:|
| no-decomp against naive | 15 105 | 15 105 | 10 171 (67.3%) |
| no-decomp-paper-boundary against paper-boundary | 15 102 | 15 098 | 10 224 (67.7%) |
| no-decomp-paper against paper | 15 057 | 15 053 | 9 769 (64.9%) |

The 4 cells whose count differs are the same 4 for both grid variants, all at m2:

| trajectory | ε | δ | direct | twin |
|---|---:|---:|---:|---:|
| geolife-007395 | 0.5 | 2.34e-5 | 890 | 888 |
| mot-001680 | 0.5 | 0.158 | 481 | 479 |
| mot-000208 | 0.5 | 0.493 | 68 | 66 |
| mopsi-002662 | 0.05 | 2.3e-5 | 306 | 308 |

All four are on `conv(G_v)`, which both grid variants share. `mot-000208` is the run
convex-shadow-32.md §12 traces: the hull's survivors differ from naive's at 115 of its vertices.
The benchmark does not decide the error. §12's 4 481 runs do: none went over the bound, and 4 469 were within
3.4e-7 δ of naive.
