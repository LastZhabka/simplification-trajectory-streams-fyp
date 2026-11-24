# Benchmark: the four implementations of ours

How fast each greedy stabber is, on the same inputs, as ε shrinks. All four are Algorithm 1
of *Faster δ-simplification of Trajectory Streams*, or Algorithm 2 which computes the same
thing faster; [convex-shadow-32.md](convex-shadow-32.md) §10–11 describes them:

| `--ours` | short | region | starts | per vertex |
|---|---|---|---|---|
| `greedy-stabbing-shared` | shared | m-gon | M = ⌈2π/ε⌉ on the δ-circle | O(M + m²) = O(ε⁻¹), Algorithm 2 |
| `greedy-stabbing-naive` | naive | m-gon | the same | O(Mm) = O(ε^(−3/2)), Algorithm 1 |
| `greedy-stabbing-naive-paper-boundary` | paper-boundary | `conv(G_v)` | grid points of the cells of `G_v` meeting the circle, O(ε⁻¹) | O(ε^(−3/2)) |
| `greedy-stabbing-naive-paper` | paper | `conv(G_v)` | every grid point of `G_v`, O(ε⁻²) | O(ε^(−5/2)) |

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

About 25 minutes with the defaults on the machine below. Results are written only at the end.
Beside `--out` it writes:

| file | what |
|---|---|
| `benchmark.json` | every cell, the selection and the settings |
| `times/<ours>.json` | every time of one implementation, grouped by trajectory: each δ and ε, its status, the best and every repeat, the output size |
| `summary/<ours>.json` | that implementation's distributions per ε, over all cells and per dataset: seconds, ns per input vertex and output size, as count, min, p10, p25, median, p75, p90, p99, max and mean; and how many finished, hit the limit or failed |
| `pairwise.json` | for each pair a, b, per ε, over all cells and per dataset: the distribution of b's time over a's, how many times faster a is. `finished` counts only cells both finished; `with_bounds` also counts those one of them did not, with the limit standing in for its time |

and prints the tables below as Markdown; redirect it to keep them (`data/benchmark/report.md`
for the run below).

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
are O(ε⁻²). Naive over paper-boundary stays at 11–13×, both being O(ε^(−3/2)): paper-boundary
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
and paper 101× at least. All grow more slowly than their worst cases — 10×, 32×, 32× and 316× —
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
