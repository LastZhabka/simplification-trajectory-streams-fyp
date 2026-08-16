# Continuous Fréchet distance

The measure everything in [comparison.md](comparison.md) is scored against. `src/algo/frechet/`
is a port of the ACM SIGSPATIAL GIS Cup 2017 winning implementation, verified against two
independent references.

This document is the **library**: what was ported, how it was verified, and what it costs. The
driver that applies it to a corpus — `ssk_frechet`, its flags, its output layout and how to
shard it — is [pipeline.md](pipeline.md); the end-to-end walkthrough is
[running.md](running.md).

## The interface

`frechet/frechet.hpp`, header-only because it is templated on dimension the way the
simplifiers are:

```cpp
template <std::size_t D> using Point = std::array<double, D>;
template <std::size_t D> using Curve = std::vector<Point<D>>;

template <std::size_t D>
class Frechet {
 public:
  bool   within(const Curve<D>& a, const Curve<D>& b, double r) const;
  double distance(const Curve<D>& a, const Curve<D>& b, double tol) const;
};
```

**`within` is the primitive.** Deciding *is d_F(a, b) ≤ r?* is what the algorithm actually
computes, and it is exact — no tolerance, no approximation. `distance` bisects it, so it costs
a tolerance and about nine decisions. Code that only needs a yes/no answer should ask for one.

The instance holds a scratch buffer for the free-space frontline, so a driver measuring a
corpus reuses one allocation across millions of calls. That is upstream's design, kept.

## What was ported

**Source:** [`mwernerds/frechetrange`](https://github.com/mwernerds/frechetrange), MIT,
`include/frechetrange/detail/dv/frechet_distance.hpp` and `include/frechetrange/distance_sqr.hpp`.

That repository consolidates the three winning entries of the **ACM SIGSPATIAL GIS Cup 2017**,
whose task was range queries under the Fréchet distance. `dv` is the entry of **Fabian Dütsch
and Jan Vahrenhold, first place**; the consolidation is Martin Werner's, with Dütsch credited
for adapting the implementations. We take only the `dv` detail header — not the repository's
unified layer, which pulls in `boost::geometry` and which its own README describes as work in
progress.

Neither ported file needs anything outside the standard library, which is why this could come
in at all: `src/algo` takes no third-party dependencies and the build must work offline.

*Kept:* the algorithm exactly. The free-space diagram traversed one column at a time carrying
only the frontline — the reachable beginnings of the previous column's right segments — so
working memory is O(min(n, m)) rather than O(nm). The `BEGIN_NOT_REACHABLE = 2.0` sentinel and
`begin <= 1.0` reachability test. The corner pre-check. The single-point special case. The
two-sided monotone point-matching filter, run in both directions before the diagram is built,
which rejects most non-matching pairs without touching a cell. `get_reachable_begin`'s four
cases. The line-circle intersection by quadratic formula, including its `a == 0` branch for a
degenerate segment.

*Changed:* upstream is templated on a point type plus a coordinate getter and a distance
functor; this takes `Curve<D>` directly, since every curve here is
`std::vector<std::array<double, D>>`. Upstream's template recursion over dimensions became a
loop — **counting down from `D-1` to `0`, matching the recursion's order**, so the
floating-point sums are identical rather than merely equivalent. Names follow this repository.

*Deviations, both deliberate:*

- Upstream leaves the per-cell `begin`/`end` scalars **uninitialised** when the line-circle
  solve finds no intersection: it assigns `begin` and returns, leaving `end` at whatever it
  held. The value is never read in a decision — a segment whose begin is unreachable is
  skipped — but it is still uninitialised. Ours initialises them. No result changes.
- Upstream's `is_bounded_by` requires non-empty input as a precondition; ours returns instead
  of reading out of bounds.

*Dropped:* `USE_POSITIVE_FILTER`, a greedy pre-filter from Baldus & Bringmann that upstream
leaves compiled out by default, and the spatial index and grid built around the decision
procedure for range queries — we need distances, not neighbour queries.

*Not ported:* a value routine. Upstream has none; `distance` is ours, and it is a plain
bisection.

## Verified against three implementations

Not just against the code it came from. **32 000 decisions over 4 000 random curve pairs**, at
radii spanning 0.25× to 4× the true distance so both answers are exercised, including 0.999×
and 1.001× where the decision flips:

| Reference | What it is | Differences |
|---|---|---|
| `dv` | the upstream header this was ported from | **0** |
| `bb` | Baldus & Bringmann, GIS Cup **second place** — a different algorithm | **0** |
| prototype | the earlier project's own Alt–Godau, in `archive/` | **0** |

The second row is the one that matters: `bb` is an independent implementation by different
authors, so agreement is evidence about the *answer*, not merely about the transcription. The
third is the implementation this project used before and would otherwise be silently replacing.

`tests/test_frechet.cpp`, 9 cases, covers what can be computed by hand — identical curves at
0, parallel lines at their separation, a polyline against its own endpoints at 0, a spike
against its chord at the spike height, a segment against itself reversed at its full length
(the matching must stay monotone) — plus 3D, degenerate curves of one point and of repeated
points, decision and value agreeing either side of the answer, and a simplification never
lying further from its input than the largest step it skips.

## Cost

Measured on real trajectory–simplification pairs from the corpus, MinGW g++ 15.1 `-O2`:

| | ns per free-space cell |
|---|---|
| one decision | **3.76** |
| the value, bisected to 1e-9 | **34.45** |

The value costs **9.2×** the decision, not the 30–50 a naive bisection count suggests, because
the filters reject most radii long before a full traversal. A free-space diagram is `n × m`
cells, so scoring the whole corpus is **0.69 trillion cells, about 6.6 core-hours** — exact,
with no per-shortcut approximation needed.

### Where that time goes

`n × m` means cost is quadratic in trajectory length at a fixed compression rate, so it lands
almost entirely on GeoLife:

| Dataset | Gcells | one core |
|---|---:|---:|
| `geolife` | 537.0 | 308 min |
| `mopsi` | 89.6 | 51 min |
| `ngsim-us-101` | 29.2 | 17 min |
| `ngsim-i-80` | 26.6 | 15 min |
| `mot-mot20` | 4.6 | 3 min |
| `ngsim-lankershim` | 2.6 | 1 min |
| `mot-dancetrack` | 1.7 | 1 min |
| `mot-mot17` | 1.2 | 1 min |
| `ngsim-peachtree` | 1.1 | 1 min |
| **total** | **693.6** | **6.6 h** |

**GeoLife is 77% of the work**, so splitting by dataset alone achieves nothing — it finishes in
five hours gated by one core. The split has to be *within* GeoLife, which is what
`ssk_frechet --shard` is for.

### It is bounded, unlike DOTS

The important difference from the simplification sweep: this cost is **predictable in advance**.
Cells are `n × m` exactly and each costs a known 34.45 ns, so the table above is arithmetic
rather than extrapolation. The worst single trajectory in the corpus, `geolife-001523` at
92 645 points, comes to 25.3 Gcells — **15 minutes on one core** for all its rates and
algorithms together.

Nothing here can stall. There is no search, no frontier and no data-dependent branching, so
unlike DOTS — where cost depended on a trajectory's *geometry* rather than its length, and one
64 000-point track consumed 90+ CPU-minutes without finishing — the bill is known before the
job starts. Working memory is O(min(n, m)) for the frontline, a few hundred kilobytes whatever
the input.

### The bottleneck is I/O, not arithmetic

The 6.6 core-hours are compute. The job also reads **every one of the ~786 000 result documents
and their inputs**, which is the access pattern that made the corpus audit take about 2.5 hours
of wall clock across nine parallel processes while running at 5% CPU.

So a full corpus pass should be planned as **roughly three hours of wall clock**, with the ~50
minutes of eight-core compute hidden inside it — not the 50 minutes alone. That figure is
extrapolated from the audit rather than measured for this pipeline, so treat it as an
expectation to check on the first dataset, not a promise.
