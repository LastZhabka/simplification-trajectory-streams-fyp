# ConvexShadow32: a regular m-gon instead of conv(G_v)

> **§1–9 describe implementations that are gone** (UPDATES 65): `ConvexShadow32`,
> `ConvexShadow32SideNorms` and `ConvexShadowInvertLin`, as of commit 902a3c5 (the last kept in
> `archive/removed-simplify/`). The m-gon argument of §1–5 is what the greedy stabbers of
> §10–11 use.

`ConvexShadow32` (`src/algo/simplify/convex_shadow_32.{hpp,cpp}`, `--ours convex-shadow-32`)
is `ConvexShadow` with one change: the region around each input vertex `v` is not `conv(G_v)`,
the hull of a grid, but `P_v`, a **regular m-gon** around `v`. This page says what the region
has to satisfy, why a regular m-gon satisfies it with **m = O(ε^(−1/2))** sides, why no convex
polygon can do with fewer, and what that is against `conv(G_v)`'s corners. The algorithm itself
is in [algorithm.md](algorithm.md); everything there except the region still holds.

---

## 1. What the region has to satisfy

Write `r = (1+ε)δ` and `ρ = (1+ε/2)δ`. The region `P` around a vertex `v` must lie between two
circles around `v`:

```
B(v, ρ)  ⊆  P  ⊆  B(v, r)
```

- **`P ⊆ B(v, r)` gives the error bound.** Every emitted segment stabs the regions of its
  run's vertices in order, so every point matched to `v` is within `r` of it
  ([algorithm.md](algorithm.md) §3, *The error bound*).
- **`B(v, ρ) ⊆ P` gives the size bound.** Lemma 2's witness, snapped to a start, stabs the
  balls of radius `δ + εδ/4` around the run's vertices. Those balls must lie inside the regions
  for the witness's shadow to stay non-empty, and `δ + εδ/4 ≤ ρ`. (The guarantee needs only
  `δ + εδ/4`; `ρ` is what `conv(G_v)` contains, so the two regions are held to the same check,
  `region_covers_ball`.)

`conv(G_v)` lies between exactly these circles: it contains `B(v, ρ)` by construction, and its
corners are within `ρ` plus a cell diagonal, `εδ/2`, of `v`. Any convex polygon between them
does the same job. Nothing about Lemma 1 or Theorem 3 uses the grid except through these two
inclusions and the starts.

**The room between the circles is a ring of width `r − ρ = εδ/2`**, and how many sides a
polygon needs to fit in it is the whole question.

## 2. A regular m-gon fits with m = O(ε^(−1/2))

Take the regular m-gon with circumradius `r` around `v`. Its vertices are on the outer circle,
so `P ⊆ B(v, r)`. Its edges are chords of that circle, each subtending an angle `2π/m` at `v`,
and the point of an edge nearest `v` is its midpoint, at distance `r·cos(π/m)`: that is the
inradius. So

```
B(v, ρ) ⊆ P   ⇔   r·cos(π/m) ≥ ρ   ⇔   cos(π/m) ≥ ρ/r = (1 + ε/2)/(1 + ε) = 1 − t,

                                           where t = (r − ρ)/r = ε / (2(1 + ε)).
```

**The bound.** `1 − cos x ≤ x²/2` for every `x`, so the condition holds as soon as
`(π/m)²/2 ≤ t`, that is

```
m ≥ π / √(2t) = π·√((1 + ε)/ε)   ≤   π·√(2/ε)   for ε < 1.
```

The fewest sides that work, `m* = ⌈π / arccos(ρ/r)⌉`, is therefore at most
`⌈π·√((1+ε)/ε)⌉ = O(ε^(−1/2))`.

**Why the exponent is ½: the sagitta.** An edge dips below the outer circle by its sagitta,
`r − r·cos(π/m) = r·(1 − cos(π/m)) ≈ r·π²/(2m²)`. The dip falls with the *square* of the number
of sides, and it must be at most the ring's width, `εδ/2`, which is linear in ε. So
`π²/(2m²) ≈ ε/2`, and `m ≈ π/√ε`. A polygon with straight edges approximates a circle to second
order; that is the whole reason.

## 3. No convex polygon does with fewer

Let `Q` be any convex polygon with `B(v, ρ) ⊆ Q ⊆ B(v, r)`, and take one of its edges, on a
line `ℓ`.

- `Q` is on one side of `ℓ` and contains `B(v, ρ)`, so `ℓ` is at distance `d ≥ ρ` from `v`.
- Both ends of the edge are within `r` of `v`. Seen from `v`, a point of `ℓ` within `r` is at
  most `arccos(d/r) ≤ arccos(ρ/r)` from the foot of the perpendicular.
- So each edge spans an angle of at most `2·arccos(ρ/r)` at `v`. `v` is inside `Q`, so the
  angles of all the edges add up to `2π`.

Hence every such polygon has at least `π / arccos(ρ/r)` edges, which is `m*`: **the regular
m-gon is optimal among all convex polygons between the two circles**, `conv(G_v)` included.

The count is also Θ(ε^(−1/2)) from below. `1 − cos x ≥ 2x²/π²` on `[0, π]`, so
`arccos(1 − t) ≤ π·√(t/2)`, and

```
2·√((1 + ε)/ε)   ≤   m*   ≤   ⌈π·√((1 + ε)/ε)⌉.
```

## 4. Against conv(G_v)

`conv(G_v)` is the hull of the corners of the grid cells, side `h = εδ/(2√2)`, that meet a
disk of radius `ρ`; in grid units that disk has radius `ρ/h = Θ(1/ε)`. The convex hull of the
integer points in a disk of radius `n` has Θ(n^(2/3)) vertices (Balog & Bárány, 1991), and the
measured counts below grow at that rate. So `conv(G_v)` has Θ(ε^(−2/3)) corners where the
m-gon has Θ(ε^(−1/2)) sides; the grid's corners are forced by the lattice, not by the ring.

| ε | starts (both) | `conv(G_v)` corners | m* | `2√((1+ε)/ε)` | `⌈π√((1+ε)/ε)⌉` | m-gon inradius / δ | needed: `1+ε/2` |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.9 | 72 | 8 | 5 | 2.9 | 5 | 1.537 | 1.450 |
| 0.5 | 120 | 12 | 6 | 3.5 | 6 | 1.299 | 1.250 |
| 0.3 | 168 | 16 | 7 | 4.2 | 7 | 1.171 | 1.150 |
| 0.1 | 472 | 28 | 11 | 6.6 | 11 | 1.055 | 1.050 |
| 0.05 | 920 | 44 | 15 | 9.2 | 15 | 1.027 | 1.025 |
| 0.01 | 4 552 | 136 | 32 | 20.1 | 32 | 1.005 | 1.005 |
| 0.005 | 9 080 | 252 | 45 | 28.4 | 45 | 1.003 | 1.003 |

From ε = 0.5 to 0.005, a factor of 100, the grid hull's corners grow 21× (exponent 0.66) and the
m-gon's sides 7.5× (exponent 0.44, tending to ½ as `1 + ε → 1`). In every row `m*` equals the
upper bound `⌈π√((1+ε)/ε)⌉`.

**Why the count matters.** Every shadow is the region of the vertex just read, cut by lines
through the previous shadow. The region's edges are where every shadow starts, so fewer of
them means shadows with fewer vertices, and every step's `behind` and every validation check
costs less.

## 5. What does not change

- **The starts.** They are still the grid points of the boundary cells, `O(1/ε)` of them, with
  the same argument behind them ([algorithm.md](algorithm.md) §3, *Starting from the boundary
  only*). Some are farther from `v` than the m-gon's inradius and lie outside `P_v`. That is
  harmless: a start is matched to its run's first vertex by distance, which is at most `r`, and
  the run's first step gives every start the whole region of the next vertex.
- **Both guarantees**, `d_F ≤ (1+ε)δ` and `|σ| ≤ 2κ − 2`, by the same proofs with `P_v` for
  `conv(G_v)`.
- **Everything else in the code.** `convex_shadow_32.cpp` is `convex_shadow.cpp` with the
  region replaced; the restart rule, the emitted points and the 16-ulp margin are the same,
  and `behind` cuts by the same lines. It finds the two tangent vertices more simply (§6).

## 6. In the code

`ConvexShadow32::sides(ε)` computes `m*` from `cos(π/m) ≥ (1+ε/2)/(1+ε)`, adjusting up or down
from `⌈π/arccos⌉` so that rounding in `acos` cannot leave it a side off. It depends on ε alone.
The m-gon's vertices are `r·(cos 2πi/m, sin 2πi/m)` with `r` taken 4 ulps short, since every
test near the error bound errs inward and `cos` and `sin` round.

`--validate` checks it with `region_convex`, `region_covers_ball` and `region_inside_disk`, the
checks `conv(G_v)` gets too; `grid_hull` does not apply and is not counted.

**The tangents come from the facing test.** `faces` is the sign of `cross(b − a, p − a)`, the
dot product of `p − a` with the edge's outward normal, widened by the same 1e-9 tolerance as in
`convex_shadow.cpp`. The tangents from `p` touch the shadow where a facing edge meets one facing
away, so the one pass that cuts by the facing edges also finds them: the vertex where facing
starts is the left tangent's, the vertex where it stops is the right one's. `convex_shadow.cpp`
orders the vertices by angle around `p` instead, split at the direction of the shadow's centre
([algorithm.md](algorithm.md) §4, *Finding the tangents near a corner*). Here no centre and no
ordering are needed. From a `p` an ulp off a corner, both edges at that corner are within the
tolerance and face `p`, so the tangents go to their far ends and never along the noise
direction from `p` to the corner. When every edge faces `p`, the shadow is in line with it,
there is no vertex where facing changes, and the cuts by its edges leave no more than the
shadow, which is in `F`.

Checked against the angular version (2026-11-17) on the first 150 Mopsi and 40 GeoLife
trajectories, ε ∈ {0.5, 0.3, 0.1}, δ at 0.0005, 0.003 and 0.02 of each trajectory's extent,
validated: 1 701 runs, equal vertex counts in all, 1 620 identical to the bit and the rest
within 6e-9 δ, 0 over the bound by the exact decision, 0 failures in 155 163 318 checks.

## 7. Measured

Against `convex-shadow` (2026-11-16), on the same baseline operating points: 1 260 matched
runs on samples of Mopsi, GeoLife and DanceTrack, over five (dataset, ε) cells, all run with
`--validate`.

| | `convex-shadow` | `convex-shadow-32` |
|---|---:|---:|
| over the bound, by the exact decision | 0 | 0 |
| over `2κ − 2`, wherever a curve on disk bounds κ | 0 | 0 |
| validation failures | 0 in 121 528 577 | 0 in 120 041 121 |
| output vertices, all 1 260 runs | 33 040 | 33 312 (+0.8%) |
| time, 40 Mopsi trajectories, ε = 0.5 | 18.2–18.8 s | 12.3–12.6 s (1.5× faster) |
| time, 12 GeoLife trajectories, ε = 0.5 | 37.7–37.8 s | 24.8–25.4 s (1.5×) |
| time, 6 Mopsi trajectories, ε = 0.1 | 18.4–18.7 s | 8.1–8.2 s (2.3×) |

Times are of `ssk_match`, one process at a time, unvalidated, two runs each; they include
measuring every output's Fréchet distance, so the algorithm alone gained more. The gain grows as
ε falls, as §4 predicts: 6 sides against 12 corners at ε = 0.5, 11 against 28 at ε = 0.1.

**Sizes, run by run:** equal in 1 043, the m-gon larger in 161 and smaller in 56. Per ε, the
mean of its size over `convex-shadow`'s is 1.023 at 0.5, 1.004 at 0.3 and 1.018 at 0.1, and the
median is 1 everywhere. The largest ratio is 2, on runs of a few vertices.

**Why a few more vertices.** Both regions hold the same inner ball and stay in the same disk,
so both guarantees are the same, but the m-gon is the smaller of the two: it has 90% of
`conv(G_v)`'s area at ε = 0.5, 96% at 0.3 and 97% at 0.1. A smaller region empties a shadow a
little sooner, so a run ends a little earlier. The size bound is a bound on the worst case,
and neither region comes near it here.

## 8. ConvexShadow32SideNorms: a shadow as m offsets and two rays

`ConvexShadow32SideNorms` (`src/algo/simplify/convex_shadow_32_sidenorms.{hpp,cpp}`,
`--ours convex-shadow-32-sidenorms`) is the same algorithm on the same m-gon, with the same
starts, restart rule and emitted points. It keeps each shadow as what bounds it: a candidate is
`p`, the m offsets, the two rays that cut its shadow and the two tangents aimed from it for the
next step. No polygon is stored.

**Why m offsets and two rays are enough.** Every region is the m-gon moved to a vertex, so it is
`nₖ·x ≤ nₖ·v + oₖ` with the same m outward unit normals `nₖ`. `F(S, p)` adds only the two
tangents from `p`, and every tangent passes through `p`. So every shadow is

```
S = { x : nₖ·x ≤ cₖ for every k }  ∩  W,      W the wedge at p between two rays, or the plane
```

**One step**, for a start `p` and the next vertex `v`:

1. `F` keeps offset `k` iff `nₖ·p ≥ cₖ`. If `p` is beyond or on a side, going away from `p` only
   goes deeper inside it: for `y = p + λ(s − p)`, `λ ≥ 1`, `nₖ·y ≤ nₖ·s ≤ cₖ`. It drops the rest.
2. Each new offset is the minimum of the kept one and the region's, `nₖ·v + oₖ`.
3. The previous tangents become the rays, and the corners are built: the region cut by every
   kept offset stricter than the region's, then by the two rays. Fewer than three corners means
   the shadow is empty.
4. If `p` is beyond a side, new tangents are aimed through the corners furthest right and left
   as seen from `p`. The corners are then dropped. Emission and the validator rebuild them the
   same way from the m-gon of the vertex the offsets were last taken at.

**Why the rays and the tangents are kept apart.** Both bound the same set, and the first version
rebuilt a shadow from its tangents. A tangent can run almost along an edge of the shadow, when
`p` is nearly in line with that edge. Clipping by it again is then ill-conditioned and moved the
corner there by about 3e-9 δ. In 2 of 51.7 M checks the next shadow had a corner in the lost
sliver (`shadow_behind`, ε = 0.1). Rebuilt from the rays that cut it, a shadow is the same
polygon to the bit.

**Two different questions about a side.** Whether `p` is beyond side `k` decides whether `F` keeps
its offset; it changes once, from inside to beyond, and never back. Whether side `k` is an edge
of `S` decides the corners; it changes from step to step, because a stricter offset can slide a
side past its neighbours' corners, and a ray cuts sides away. So the corners are not where
neighbouring side lines meet, and each step has to compute them. On 46.5 M candidate steps of the
sample below (ε = 0.5 and 0.1), an m-gon side was gone from the offsets alone in 27% and once
the rays cut in 59%. Each step makes 2.58 clip passes on average, almost all of them the two
rays.

**Once outside, `p` stays outside, exactly.** A side `p` is beyond is kept, its offset only
falls, and `nₖ·p` is the same double at every step. So "is `p` in its shadow" compares `nₖ·p`
with `cₖ`, against normals that are constants rather than differences of shadow vertices. No
tolerance, no exact-equality case and no final hull are needed. `ConvexShadow32` needs all three
because it reads each edge's direction off the polygon.

**What reads the corners** is the tangents and the emitted end point. The tangent vertices
are ordered on each side of a reference direction, as in `convex_shadow.cpp`; the reference is
the inward normal of the side `p` is furthest beyond, which the whole shadow is ahead of. The end
point's `closest` skips an edge a few ulps long. A ray passing within an ulp of a corner leaves
such an edge, its direction is noise, and believing it put a stationary track's last vertex
outside its own shadow: in 2 of the runs below, before the skip, the end point moved by the
inradius.

**Measured** (2026-11-17) on the same 1 701 validated runs as §6:

| | `convex-shadow-32` | `convex-shadow-32-sidenorms` |
|---|---:|---:|
| output vertices | 63 355 | 63 355, equal in every run |
| outputs identical to the bit | | 1 580; the rest within 2.3e-9 δ |
| over the bound, by the exact decision | 0 | 0 |
| validation failures | 0 in 155 163 318 | 0 in 155 163 318 |
| time at ε = 0.5, unvalidated, two runs | 10.7–12.1 s | 6.6–7.6 s (1.6×) |
| time at ε = 0.1 | 61.5–62.2 s | 34.0–35.1 s (1.8×) |
| lines of `.cpp` | 336 | 299 |

Storing no polygon changed no output: the same 1 701 runs are identical to the bit with the
version that kept one, and run in the same time within noise (29.4–31.8 s against 30.5–34.5 s
at ε = 0.1, three runs each).

## 9. ConvexShadowInvertLin: shared minima and per-start cones, O(ε⁻¹) per vertex

`ConvexShadowInvertLin` (`src/algo/simplify/convex_shadow_invert_lin.{hpp,cpp}`,
`--ours convex-shadow-invert-lin`) implements Algorithm 1 of the manuscript *Planar Streaming
Simplification in O(ε⁻¹) Space and Update Time* (16 November 2026). It inverts §8: the m
offsets are shared by every start, and each start keeps only a cone.

**What it stores**, for the open run:

| | size |
|---|---|
| the run's first vertex `v_s` and the M start offsets on the circle of radius `ρ` | O(M) |
| `minima`: per normal, the smallest offset of the run's m-gons so far, `Hᵢ` | m |
| per start: alive, and a cone at `p`, the whole plane or two rays | O(M) |
| the buffered segment `(p, q)` | O(1) |

No offsets per start and no polygon per start: a side `p` is strictly beyond always has offset
`Hᵢ` (Lemma 3 of the manuscript), so one shared vector serves every start.

**What it computes for a vertex `v`**, with incoming offsets `aᵢ = nᵢ·v + h`:

1. **Blocks.** For each of the 2m thresholds `Hᵢ` and `aᵢ`, the starts strictly beyond it form
   one arc of the circle; two binary searches, on the two monotone halves of `nᵢ·p`, find its
   ends. The at most 4m ends cut the starts into at most 4m + 1 blocks, on which every sign is
   the same.
2. **Per block**, at its first start: `bᵢ = min(aᵢ, Hᵢ)` where `Hᵢ < nᵢ·p`, `aᵢ` elsewhere, and
   `B` their intersection, by one deque pass over the normals in order, O(m). The corner where
   the sides facing `p` start and the one where they stop are the tangent contacts, the same for
   the whole block.
3. **Per start**, O(1): the cone is cut to the rays through those two corners; an empty or
   zero-width cone kills the start. A start inside `B` keeps its cone.
4. **The first survivor**: `B` cut by its cone is its shadow; `q` is its point nearest `v`, and
   `(p, q)` becomes the buffer.
5. Survivors: `Hᵢ ← min(Hᵢ, aᵢ)`. None: the buffer is emitted and a run starts at `v`.

That is O(m log M + m log m + m² + M) per vertex, O(ε⁻¹) with M = Θ(ε⁻¹) and m = Θ(ε^(−1/2)).

**What differs from the other m-gon implementations.** The starts are the manuscript's:
M = ⌈4π(1 + ε/2)/ε⌉ points on the circle of radius `ρ` around `v_s`, 32 at ε = 0.5 and 132 at
ε = 0.1, against 120 and 472 boundary grid points. Comparisons with `p` are strict. All starts
are at the same distance from `v_s`, so the first survivor is emitted. `--validate` checks the
region and the segments only; there are no per-start shadows to check.

**Measured** (2026-11-17), on the 1 701 validated runs of §6:

- Against a brute-force reference that runs the manuscript's recurrence with an explicit
  polygon per start (the same starts and strict rules; it passes 0 failures in 44 257 693
  checks): equal vertex counts in every run, 1 597 identical to the bit and the rest within
  7.6e-10 δ.
- 0 over the bound by the exact decision; 0 failures in 37 187 checks.
- 63 819 vertices against `convex-shadow-32`'s 63 355 (+0.7%; larger in 191 runs, smaller in
  1), the price of the sparser starts.
- Unvalidated, two runs each:

| | `convex-shadow-32` | `convex-shadow-32-sidenorms` | `convex-shadow-invert-lin` |
|---|---:|---:|---:|
| ε = 0.5 | 9.7–9.9 s | 5.6 s | 0.61–0.66 s |
| ε = 0.1 | 57.7–58.2 s | 29.3–30.2 s | 1.32–1.35 s |

## 10. GreedyStabbingNaive: the naive planar algorithm, O(ε^(−3/2)) per vertex

`GreedyStabbingNaive` (`src/algo/simplify/greedy_stabbing.{hpp,cpp}`, `--ours
greedy-stabbing-naive`; `greedy-stabbing` before entry 62)
implements Algorithm 1 of the write-up *Faster δ-simplification of Trajectory Streams*, its
naive planar algorithm. The region is §6's m-gon; each start keeps its own offsets, as in §8,
and a cone instead of a polygon, as in §9.
It is the first `GreedyStabber` (`src/algo/simplify/greedy_stabber.hpp`): the loop over the
vertices, the restarts and the last segment are the base's, and the steps below its overrides.
Its geometry is `src/algo/geometry/`, and it uses nothing of the other implementations: the
region is a `ConvexPolygon` and the side count its own, the same m as `ConvexShadow32::sides`.

**What it stores**, for the open run:

| | size |
|---|---|
| per start `p`: the m offsets `d_k`, the minimum of `b + n_k·v` over the run's vertices after its first | O(Mm) |
| per start: `R(S[p], p)`, the whole plane or two rays from `p` | O(M) |

**What it computes for a vertex `v`**, for each start, O(m):

1. `d_k ← min(d_k, b + n_k·v)`, with `b` the m-gon's apothem.
2. `Q[p]`: offset `d_k` on the sides `p` is strictly beyond (`V[p]`), `b + n_k·v` on the rest,
   by one deque pass over the normals in angular order.
3. `R(Q[p], p)`: the whole plane if `V[p]` is empty; otherwise the rays through the two corners
   where the sides in `V[p]` start and stop; nothing if `Q[p]` is empty or has no area.
4. `R ← R ∩ R(Q[p], p)`. An empty or zero-width cone kills the start.

With no survivor, the first start that survived the step before is emitted, with the point of
`S[p] = Q[p] ∩ R` nearest the run's last vertex, and a run starts at `v`. That is O(Mm) per
vertex, O(ε^(−3/2)) with M = Θ(ε⁻¹) and m = Θ(ε^(−1/2)).

**What differs from the other m-gon implementations.**

- **The starts**: M = ⌈2π/ε⌉ points on the circle of radius `δ` around the run's first vertex.
  That is 13, 21 and 63 at ε = 0.5, 0.3 and 0.1, against §9's 32, 49 and 132 on the circle of
  radius `ρ`. Every point of the `δ`-circle is within `εδ/2` of a start, and the greedy may
  start on the circle instead of in the disk, so the witness stabs balls of radius
  `(1 + ε/2)δ = ρ`, which the m-gon contains.
- **The normals** are at angles `2πk/m`, not §6's `π(2k+1)/m`. Any rotation will do.
- **The paper's region and starts**, the same algorithm otherwise. `GreedyStabbingNaivePaper`
  (`--ours greedy-stabbing-naive-paper`) takes both from `G_v`, the cells of side `εδ/(2√2)`,
  with `v` at a corner, that meet the disk of radius `ρ`. The region is `conv(G_v)`, its hull,
  taken in cells, where grid points in a line are exactly in a line, then scaled and turned
  into its edges' half-planes. The starts are every grid point of `G_v`.
  `GreedyStabbingNaivePaperBoundary` (`--ours greedy-stabbing-naive-paper-boundary`) has the same
  region and keeps only the corners of the cells that meet the circle, the starts of §6–8. All of
  it is built in `greedy_stabbing.cpp`, not taken from those implementations, and the starts are
  ordered nearest `v` first, so the start emitted is the nearest that survived:

  | ε | m-gon sides | `conv(G_v)` sides | circle of radius `δ` | every grid point | boundary cells |
  |---|---:|---:|---:|---:|---:|
  | 0.5 | 6 | 12 | 13 | 225 | 120 |
  | 0.3 | 7 | 16 | 21 | 461 | 168 |
  | 0.1 | 11 | 28 | 63 | 3 009 | 472 |
- **`--validate` checks nothing yet**: the validator was taken out to keep it small, and a
  validated run reports 0 checks. The numbers below were measured with it in (entry 58).

**Measured** (2026-11-23), on the first 150 Mopsi and 40 GeoLife trajectories, ε ∈ {0.5, 0.3,
0.1}, with δ at 0.0005, 0.003 and 0.02 of each trajectory's extent: 1 710 runs.

- Against a brute-force reference with an explicit polygon per start, stepped as
  `S ← A_v ∩ F(S, p)` from that polygon's own edges and tangents (no offsets, no `V`, no cones;
  the same starts and emission): equal vertex counts in every run, 1 614 identical to the bit and
  the rest within 4.3e-8 δ. The reference at first disagreed in one run. The fault was its own:
  a cut through a corner left an edge a few ulps long, whose direction is noise. Once it skips
  such edges, the two agree. `greedy_stabbing.cpp` never uses a shadow's edges.
- 0 over the bound by the exact decision; 0 failures in 28 348 997 checks of the region, every
  start's shadow and the segments.
- Vertices and time, unvalidated, one run each, all four in the same session:

| | `convex-shadow-32` | `convex-shadow-32-sidenorms` | `convex-shadow-invert-lin` | `greedy-stabbing-naive` |
|---|---:|---:|---:|---:|
| ε = 0.5 | 21 194, 13.9 s | 21 194, 8.3 s | 21 416, 0.81 s | 21 774, 1.89 s |
| ε = 0.3 | 22 352, 24.9 s | 22 352, 11.6 s | 22 503, 1.03 s | 22 709, 3.11 s |
| ε = 0.1 | 23 792, 84.5 s | 23 792, 43.5 s | 23 862, 1.69 s | 23 960, 10.3 s |

That is +2.7%, +1.6% and +0.7% vertices against `convex-shadow-32`, the price of the sparser
starts, at 3.7–4.4× the speed of `convex-shadow-32-sidenorms`.

## 11. GreedyStabbingShared: Algorithm 2, O(ε⁻¹) per vertex

`GreedyStabbingShared` (`src/algo/simplify/greedy_stabbing.{hpp,cpp}`, `--ours
greedy-stabbing-shared`) implements Algorithm 2 of the same write-up: §10 with the computations
shared between starts. Same region, starts and emission; another `GreedyStabber`.

**What it stores**, for the open run:

| | size |
|---|---|
| the run's first vertex `c` | O(1) |
| the m offsets `d_k`, shared by every start | m |
| per start: its index on the circle, `p`, and `R(S[p], p)` | O(M) |

**What it computes for a vertex `v`**:

1. `d'_k = min(d_k, b + n_k·v)`.
2. **Marks.** Each line `n_k·x = d'_k` meets the circle of radius `δ` around `c` in at most two
   points, found by angle: `cos(θ − α_k) = (d'_k − n_k·c)/δ`. Only there can the set of sides a
   start is beyond, `V[p]`, change. The starts are evenly spaced, so where to look is found by
   dividing the crossing's angle by the spacing, with no sort. What marks a start is the test
   `Q[p]` is built with, `n_k·p > d'_k`, asked of the starts a spacing either side of the
   crossing: a start is marked where its answer differs from the start before's. That is at
   most one per crossing, and none for a line whose two crossings fall between the same two
   starts. The angle and the dot product round differently, so the angle only says where; the
   test decides, strict at a start exactly on the line as `Q[p]` is. Until 2026-11-24 each
   crossing marked every start within a quarter spacing of it and the next instead: 1.4 marks
   per change on the benchmark's cells, 28–31% of them spurious, 17–26% of builds wasted
   ([benchmark.md](benchmark.md)).
3. **Per start**, in order: after a mark, or at the first start, `Q[p]` and its two corners are
   built from `p`, O(m). Until the next mark they are the same, as they depend on `p` only
   through `V[p]`, and each start only forms its two rays from them and cuts its cone, O(1).
4. With survivors, `d ← d'`; without, the start that survived last emits as in §10, with the
   `d` it reached.

That is O(m) crossings, O(M) marks and starts, and at most 2m + 1 builds of O(m): O(M + m²)
per vertex, O(ε⁻¹). The text's sort of the crossings, O(m log m), is not needed.

**Against §10**: every start sees the same `Q[p]` and corners as there, so the output is meant
to be the same to the bit, and is.

**Measured** (2026-11-23), all four stabbers on the 1 710 runs of §10 (150 Mopsi and 40 GeoLife
trajectories, ε ∈ {0.5, 0.3, 0.1}, δ at 0.0005, 0.003 and 0.02 of each extent), no validator:

- **`greedy-stabbing-shared` is `greedy-stabbing-naive` to the bit** in all 1 710 runs.
- **0 over the bound** by the exact decision, for all four.
- **Against the brute-force reference of §10**, generalised to take the region as corners
  (the m-gon, or the hull of the grid taken directly) and the starts as a list, and stepping an
  explicit polygon per start: equal vertex counts in every run for all three kinds of Algorithm 1
  (`greedy-stabbing-naive-paper` at ε 0.5 and 0.3 only, 1 140 runs; the reference is slow with
  3 009 starts). The naive one is identical to the bit in 1 614 runs, the rest within 4.3e-8 δ.
  The two paper ones are identical in 1 622 and 1 076 runs, but up to 5 δ apart in a few:
- **A start on the region's boundary.** When a vertex repeats the one before, `G_v`'s grid
  points on the hull of the previous region lie exactly on the new region's edges. Whether such
  a start is beyond a side is then decided by rounding: ours may find it an ulp beyond, give it a
  cone of width about π, or kill it, where it is in the region and should see the whole plane.
  The reference, with its own predicates, keeps it. At a repeated vertex at (123.456, −78.9),
  δ = 0.7: ε = 0.5 keeps all 225 starts but 9 lose the whole plane; ε = 0.3 kills 1 of 461 and
  restricts 7 more; ε = 0.1 kills 3 of 3 009 and restricts 11. Every run more than 1e-6 δ apart
  is on a trajectory with exact repeats (in `mopsi-000111` the second point is the first, and
  the outputs differ from the first vertex). This errs inward, so the bound holds; and the
  starts the size bound relies on are well inside the region, so it holds too. The circle's
  starts are at δ, inside the m-gon, and never meet it.
- **Vertices and time**, one run each; at ε 0.3 and 0.1 beside one or two other single-threaded
  jobs on 8 cores:

| ε | `greedy-stabbing-naive` | `greedy-stabbing-shared` | `greedy-stabbing-naive-paper-boundary` | `greedy-stabbing-naive-paper` |
|---|---:|---:|---:|---:|
| 0.5 | 21 774, 2.02 s | 21 774, 0.62 s | 21 003, 24.1 s | 20 994, 47.1 s |
| 0.3 | 22 709, 3.94 s | 22 709, 0.80 s | 22 389, 37.9 s | 22 382, 114 s |
| 0.1 | 23 960, 12.3 s | 23 960, 1.67 s | 23 770, 176 s | 23 764, 1 052 s |

The sharing is 3.3×, 4.9× and 7.4× at ε 0.5, 0.3 and 0.1. The paper's `conv(G_v)` and grid
starts take 0.8–3.6% fewer vertices than the m-gon and circle, for 12–85× the naive stabber's
time and 39–630× the shared one's.
