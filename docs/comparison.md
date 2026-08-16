# Comparing against the baselines

Why the experiments fix a compression rate, what that buys, and what it costs.

## What we are actually doing

This project builds a streaming simplifier that takes a bound `delta` and returns a
simplification within that Fréchet distance of the input. To claim it is any good, it has to
be put beside the three algorithms everyone else uses — Douglas–Peucker, SQUISH and DOTS,
ported in `src/algo/simplify/` — on 47 916 real trajectories from four datasets.

The experiment that does it, in one paragraph: **take every trajectory, simplify it with each
baseline down to a fixed fraction of its vertices, measure how much Fréchet error that cost,
then run ours at exactly that error and see how many vertices it needs.** Repeat at 1/2, 1/4,
1/8, 1/16, 1/32 and 1/64 of the input, over the whole corpus. If ours consistently needs fewer
vertices for the same error, that is the result; if it does not, that is also the result.

Two artefacts come out of it. A **table** — at each compression rate, the error each algorithm
incurred — and a **curve** per algorithm through (rate, error) space, which is the same data
read the other way round.

The rest of this document is why that experiment is shaped the way it is: why a fixed
compression rate is the axis rather than a fixed error, what has to be reported alongside the
numbers for the comparison to be honest, and — at the end — exactly how each of the three
baselines is driven to a vertex count it was never designed to take.

## Four algorithms, four incomparable knobs

Nothing can be held constant across these algorithms, because no two of them take the same
kind of parameter:

| | knob | unit | what it bounds |
|---|---|---|---|
| ours | `delta` | coordinate units | Fréchet distance to the input |
| Douglas–Peucker (DPn) | `count` | vertices | output size |
| SQUISH | `buffer_size` | vertices | output size |
| DOTS | `lssd_threshold` | units²·points | summed squared SED along a shortcut |

There is no setting at which they are "run the same way". Worse, three of the four units
depend on the coordinate system — a Mopsi document is degrees, NGSIM is State Plane feet, MOT
is pixels — so even one algorithm's knob does not mean the same thing across datasets. See
[simplification.md](simplification.md) for how far that goes: DOTS' upstream default of
`10000.0` is sensible in metres and collapses a trajectory to its endpoints in degrees.

## What is comparable

Whatever produced it, every simplification has two properties that are *measured* rather than
configured:

- **compression rate** — kept vertices over input vertices,
- **error** — Fréchet distance from the simplification to the input.

That pair is the operating point, and it is the only common ground. So the experiment fixes
one of the two and measures the other.

## The protocol

**1. K-simplification.** Run each baseline at a fixed budget of `N / 2^M` vertices, for
M = 1…6. DPn and SQUISH take a budget natively and hit it exactly; DOTS has no budget knob and
needs a search on its threshold (99.3% of budgets exact, worst miss 0.4% under). The mechanics
are in *How each baseline is driven to the budget* below.

**2. Measure.** Compute the Fréchet distance from each simplification back to its input. Call
it δ*. This is where a baseline's own units stop mattering: whatever `lssd_threshold` meant,
the result is now described by a rate and a distance.

**3. Match.** Run our algorithm at `delta = δ*` and compare output sizes. Ours is guaranteed
within δ*, so the question left is the only one worth asking: **at the same error, how much
more can we compress?**

The protocol works because of an asymmetry worth stating plainly. The baselines can be set by
*size* but not by error; ours can be set by *error* but not by size. So each algorithm is
driven by the quantity it can hit exactly, and the comparison is made on the other one. No
search is needed on either side — except for DOTS, which can be set by neither.

### Where ε goes

The matching supplies **δ and nothing else**. If our algorithm also takes an approximation
parameter — trading how close to the optimal δ-simplification it gets against time or space —
that parameter is ours alone: no baseline can supply a value for it, and there is nothing to
match it against. So it is swept, not matched. Fix δ from the baseline, sweep ε over its
range, and the result is a third axis showing what the approximation costs us in output size
at a fixed error. Reporting one ε without saying which is how a result becomes unreproducible.

One consequence worth stating before the numbers exist: if ε loosens the guarantee rather than
the search — that is, if the output is within `(1 + ε)δ` rather than within δ — then matching
`delta = δ*` no longer puts us at the baseline's error, and the comparison has to use the
measured error on both sides or it is not a comparison at all.

## What has to be reported honestly

**Our guarantee is one-sided.** Running at `delta = δ*` bounds our error at δ*; it does not
make it equal to δ*. If we come in well under, "the same Fréchet distance" flatters us and
the size difference is not the whole story. Report our *measured* error next to the baseline's,
not just the δ we asked for.

**Nobody else is optimising Fréchet.** SQUISH and DOTS minimise synchronous Euclidean
distance, which is time-aware and punishes a shortcut that arrives at the right place at the
wrong time; DPn minimises perpendicular distance, and its budget form does not bound even
that. Scoring all four under the continuous Fréchet distance measures them against *our*
objective, which they never agreed to. It is the right primary metric — it is what our
algorithm guarantees — but the SED-based pair should be scored under SED too, or the result
reads as a rigged race.

**Compare like with like on the output.** All three baselines return a subsequence of the
input. If our simplification may place vertices that are not input points, it is solving a
strictly easier size problem, and that belongs in the caption rather than in a footnote.

**Say which numbers came from a search.** Bisecting DOTS to a vertex budget costs it the
online property that is the entire point of the algorithm. A table that mixes it with the
streaming algorithms without saying so is claiming something false.

**And say where a baseline is missing.** DOTS holds 271 190 of the 272 230 operating points it
should, and the shortfall is entirely its own — DPn and SQUISH are complete.

Two causes. `geolife-013552` is excluded outright because its budget search does not finish in
reasonable time, while DPn and SQUISH handle it in under a second. And **1 034 of its results
were removed** after a full-corpus audit found them not to be subsequences of their inputs, a
defect reproduced in the original `DotsSimplifier` — see
`archive/[2026-08-16] Incident - 1034 corpus documents removed.md`.

The shortfall is 0.38% overall but **19.9% in its worst cell**: at m1 on `ngsim-us-101`, DOTS
covers 2 280 of 2 847 trajectories against 2 847 for the other two. Any error averaged over
that cell is averaged over a fifth fewer trajectories for DOTS, and the missing fifth is not
random — it is where the decode misbehaved. A per-cell trajectory count belongs in any table
that puts the three side by side.

## Other protocols, and why this one

Fixing the rate is one of five defensible choices:

| protocol | search needed | gives |
|---|---|---|
| fix compression, compare error | exact budgets: DPn, SQUISH ✓, DOTS ✗ | one column per rate |
| fix error, compare compression | search for *every* baseline | one column per δ |
| **error-matched pairing** (above) | none, except DOTS' budget | direct "same error, smaller output" |
| **rate–distortion curves** | none at all | both readings at once |
| Pareto dominance | none | how often we win outright |

**Rate–distortion curves are the better frame, and cost less.** Let every algorithm sweep its
*own* knob over a natural grid — `count` for DPn, `buffer_size` for SQUISH,
`lssd_threshold` for DOTS, `delta` for ours — and plot each run as a point at (rate, error).
Each algorithm traces a curve. The horizontal gap between two curves answers "same error, who
compresses more"; the vertical gap answers "same size, whose error is lower". Neither question
needs the operating points to line up, so **no algorithm has to be searched at all** — which
deletes the DOTS bisection, the most expensive step in the pipeline by a factor of ten.

The fixed `N / 2^M` grid still earns its place as the *headline* table: a curve is the honest
picture, but "at 8× compression, error X vs Y vs Z vs ours" is the sentence a reader
remembers, and aligned rates are what make that sentence sayable. Do both — they come from
one set of runs, since a swept curve contains the fixed-rate points whenever a baseline can be
driven to them exactly.

Pareto dominance is the fallback that never needs anything to align: count the (trajectory,
operating point) pairs where our result is both smaller and no less accurate. It survives
ties, unreachable budgets and short trajectories, and it is the claim least sensitive to how
the grid was chosen.

## How each baseline is driven to the budget

The target for a trajectory of `N` points at rate M is `k = ceil(N / 2^M)`, and the three
baselines reach it in three different ways.

**Douglas–Peucker — set `count = k`. Exact, one run.** This is the reason the port is
psimpl's DPn rather than its tolerance form: DPn ranks all sub-polylines in a max-heap by
their furthest point and promotes the globally furthest one until `k` keys exist, so the
budget is the loop's termination condition and it comes out exact by construction. Measured
445/445 budgets exact. The tolerance form had to be bisected instead and landed exactly on
only 79.7% of the same budgets, missing by up to 33% on small ones.

**SQUISH — set `buffer_size = k`. Exact, one run.** The knob already is a point budget: the
algorithm keeps a fixed-size buffer and evicts the cheapest point whenever it fills, so the
output is what survives in the buffer. Measured 138/138 exact. One constraint from upstream —
`buffer_size > 4` — so M is capped at the largest value with `N / 2^M > 4`, and short
trajectories drop out of the deeper rates.

**DOTS — search `lssd_threshold`. 6.1 runs per budget.** The only one with no budget knob. Its
threshold bounds error, not size, so the size it produces has to be searched for. Three things
make that affordable:

- *The count is monotone in the threshold*, so the search is a bisection rather than a scan —
  verified over a 561-step log grid per trajectory, 0 non-monotone steps.
- *All six rates are solved as one problem, not six.* Smaller budgets need larger thresholds,
  so the rates are solved from M = 1 upward and each solved threshold becomes the lower bound
  of the next search. One traversal of the threshold range instead of six.
- *Every evaluation is memoised and the bracket is seeded by interpolation* — the counts
  already seen bound the next target and give a secant step, since count falls roughly
  linearly in log(threshold).

Together those take it from 40 runs per budget to 6.1, with identical results on all 2966
budgets measured. DOTS lands exactly on 99.3% of budgets; the rest miss low by at most 0.4%,
because the count only changes at thresholds where a shortcut becomes legal and occasionally
two counts are adjacent across one such step. Record the achieved count next to the target —
it is not always the target.

**Two floors apply to all three.** A budget below 2 is meaningless, since both endpoints are
always kept, and below 5 for SQUISH. When a rate falls under its floor, the trajectory
contributes no row at that rate rather than a clamped one — clamping would silently mix rates
within a column.

**And DOTS stops being online.** Bisecting a threshold means running the algorithm to
completion repeatedly over a stored trajectory, which is exactly what DOTS was designed not to
need. Its numbers in a fixed-rate table describe DOTS-the-simplification-quality, not
DOTS-the-streaming-algorithm.

## Cost

Measured over 600 Mopsi trajectories (778,766 points, 1.65% of the corpus), sweeping M = 1…6:

| | this sample | extrapolated to 47.1M points |
|---|---|---|
| DPn at a budget | 0.4 s | ~22 s |
| SQUISH at a budget | 2.9 s | ~3 min |
| DOTS at a fixed threshold | 6.0 s | ~6 min |
| DOTS searched to a budget | 76 s | ~1.2 h |

The last row is the whole cost of the fixed-rate protocol, and it is the reason the curve
frame is worth having: sweeping DOTS' own threshold is the ~6 min row, not the ~1.2 h one.
That row is already down from ~6 h by the three measures above, and the remaining factor is
DOTS itself, not the search around it. Cheaper still is not needing it.

None of this is worth parallelising around: the sweep is embarrassingly parallel across
trajectories, so the ~1.2 h is a core-hour figure, not a wall-clock one.
