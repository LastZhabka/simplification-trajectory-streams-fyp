# ConvexShadow32: a regular m-gon instead of conv(G_v)

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
