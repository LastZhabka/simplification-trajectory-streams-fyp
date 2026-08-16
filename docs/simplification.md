# Simplification algorithms

Three baselines to compare a Fréchet-bounded streaming simplifier against, ported from their
reference implementations behind one interface: **Douglas–Peucker**, **SQUISH** and **DOTS**.

Everything lives in `src/algo/simplify/` and builds into `ssk` with no new dependencies.

## The interface

`simplify/simplifier.hpp`:

```cpp
template <std::size_t D> using Point = std::array<double, D>;
template <std::size_t D> using Curve = std::vector<Point<D>>;

struct Context { std::vector<double> t; };     // per-point input that is not geometry

class Params {                                  // name -> value, algorithm-specific
  double get(const std::string& key) const;                 // throws when absent
  double get(const std::string& key, double fallback) const;
  void set(const std::string& key, double value);
};

template <std::size_t D>
class Simplifier {                              // every algorithm
 public:
  virtual std::string name() const = 0;
  virtual bool needs_time() const { return false; }

  Curve<D> simplify(const Curve<D>&, const Context& = {}) const;

 protected:
  virtual Curve<D> run(const Curve<D>&, const Context&) const = 0;
};

template <std::size_t D>
class SubsetSimplifier : public Simplifier<D> {  // output is a subsequence of the input
 public:
  std::vector<std::size_t> indices(const Curve<D>&, const Context& = {}) const;

 protected:
  Curve<D> run(const Curve<D>&, const Context&) const final;   // materialises the indices
  virtual std::vector<std::size_t> run_indices(const Curve<D>&, const Context&) const = 0;
};
```

Four decisions worth stating:

**Dimension is a template parameter.** An algorithm that only works in the plane derives from
`Simplifier<2>`; one that works in any dimension stays a template. That is not decoration —
Douglas–Peucker is `DouglasPeucker<D>` and genuinely runs in 3D, while DOTS and SQUISH are
2D because their reference implementations are, and saying so in the type is better than
failing at run time.

**`simplify` returns points, because that is all every algorithm can promise.** A min-# or
optimal-placement method may emit vertices that appear nowhere in the input, so points are
the only universal output. All three algorithms ported here happen to pick a *subsequence*,
which is worth more — positions keep the link back to the source record and its timestamp,
and make the correspondence an error measure needs explicit — so they derive from
`SubsetSimplifier` and additionally offer `indices`. Code that needs points takes a
`Simplifier`; code that needs provenance asks for a `SubsetSimplifier` and the compiler
checks it. Adding an algorithm that invents points means deriving from `Simplifier` and
overriding `run`, and nothing else changes.

**Hyper-parameters go through `Params`**, a `string -> double` map, so a driver can build any
algorithm from a config file without knowing which one it is. Each algorithm reads what it
needs in its constructor and a missing required key is an error, not a silent default.

**Extra inputs go through `Context`.** Today it holds `t`, because the SED-based algorithms
need a clock and nothing else needs anything. `needs_time()` lets the base class reject a
missing or mismatched clock once, in one place, rather than in each algorithm.

Two bridges connect this to the on-disk format:

```cpp
template <std::size_t D> Curve<D> curve_of(const io::Document&);   // throws on a dim mismatch
Context context_of(const io::Document&);
```

## Usage

```cpp
#include "simplify/douglas_peucker.hpp"

const auto doc = ssk::io::read_trajectory_file("data/trajectories/mopsi/mopsi-000001.json");
const auto curve = ssk::simplify::curve_of<2>(doc);

ssk::simplify::DouglasPeucker<2> dp({{"tol", 1e-4}});
const auto kept = dp.indices(curve);              // positions into `curve`
const auto out  = dp.simplify(curve);             // or the points themselves
```

| Algorithm | Header | Parameters | Needs `t` | Dimensions |
|---|---|---|---|---|
| Douglas–Peucker | `simplify/douglas_peucker.hpp` | `tol` | no | any |
| SQUISH | `simplify/squish.hpp` | `buffer_size` (> 4) | yes | 2 |
| DOTS | `simplify/dots.hpp` | `lssd_threshold`, `k` = 2.0, `max_vk_size` = 1e6 | yes | 2 |

## Sources, and what changed

Both upstreams were cloned and ported by hand. Neither is vendored — the ported files are
ours, the algorithms are theirs.

### Douglas–Peucker — [psimpl](https://psimpl.sourceforge.net/) v7 (MPL 1.1)

From `simplify_douglas_peucker` and `DPHelper` in `psimpl.h`.

*Kept:* the part most reimplementations omit — psimpl runs a **radial-distance pass at the
same tolerance first** and applies DP to that reduced polyline. Also the point-to-**segment**
distance (not point-to-line), all comparisons in squared distance, the explicit LIFO stack
instead of recursion, and the "key index 0 means none found" sentinel.

**The radial pass costs the error bound**, and that matters here more than it does to psimpl.
DP alone keeps every point within `tol` of its simplification; the radial pass first drops
points within `tol` of the *previous surviving vertex*, and those two errors add. Measured
over 200 random walks of 300 points: 1.88 `tol` with the pass, exactly 1.00 `tol` without it.
So `tol` is not a Hausdorff — let alone a Fréchet — bound on this baseline, and comparing it
against a δ-bounded simplifier means measuring the error rather than assuming it.

*Changed:* psimpl takes a flat iterator range of interleaved coordinates and writes points;
this takes `Curve<D>` and returns indices, so the radial pass carries surviving positions
rather than copying coordinates. `util::scoped_array` became `std::vector`. One numerical
difference: psimpl computes the projection fraction in `float` deliberately, to dodge
integer-type division for integer polylines; we are always `double`, so it stays `double`.
That is marginally more accurate and can move a borderline point.

*Dropped:* the other seven algorithms in psimpl, the `_n` (point-count) variant, and the
error-estimation helpers.

*Why psimpl and not the dots repo's own DP:* `DouglasPeuckerBatchSimplifier::simplifyByIndex(
x, y, out, minimumDistance)` is fixed at 2D and needs Qt. psimpl is header-only, dependency
free and n-dimensional, so it is the only one of the three baselines that can be compared
against a k-d simplifier. psimpl is also pure geometry — it has no notion of a trajectory,
so unlike the dots repo it applies no projection, normalisation or time handling of its own.

### SQUISH — [caoweiquan322/dots](https://github.com/caoweiquan322/dots), `SquishBatchSimplifier`

*Kept:* the structure exactly. An ordered map of surviving points and a priority set keyed
`(SED, index)`; priorities as **integers scaled by 1000**, so ties break on index as upstream;
the first point pinned at `1 << 30` so it can never be evicted; on eviction the removed
point's error is added to *both* neighbours, skipping the predecessor when it is index 0 and
re-keying the successor only if it is actually in the buffer; the `size >= buffer_size` fill
test; the `buffer_size > 4` precondition.

*Changed:* `QMap` → `std::map`, `QMap<QPair<int,int>, bool>` → `std::set<std::pair<int,int>>`
(the mapped `bool` was never read, so the set is the honest structure), upstream's `itr--` /
`itr += 2` → `std::prev` / `std::next`, `qRound` → `std::lround`. Output comes straight off
the ordered map, which upstream sorted explicitly.

### DOTS — [caoweiquan322/dots](https://github.com/caoweiquan322/dots), `DotsSimplifier`

*Kept:* the DAG search in full — both the streaming and the `finished` branch, `getLSSD`'s
O(1) closed form over its eight prefix sums, `needUpdateVK` / `updateVK` / `minimizeISSED` /
`viterbiDecode`, and the numerical guards in `feedData`. The **online interface is kept**
(`feed`, `read_index`, `finish`) because being online is the entire point of the algorithm;
`run()` drives it exactly as upstream's `batchDotsByIndex` does, reading at most one index
per fed point and draining after `finish()`.

*Changed:* Qt containers → standard ones. `vK`/`vL` were `QVector<double>` holding integer
indices and are now `std::vector<int>`; likewise `double minParent`. Upstream reaches its
data through `pX`, `pXSum`, … pointers that exist only to support cascade mode, so with
cascade dropped those became direct member access.

**The one behavioural deviation:** `run()` translates the curve to its own first point in
space and time before feeding. Upstream's guards reject a first timestamp above a year in
seconds or a coordinate above 2 000 km — and our documents carry epoch-millisecond timestamps
around 1.2e12, so DOTS would refuse to run on them outright. Upstream never hits this because
`Helper` normalises during parsing (see *Upstream works in projected metres* below); since
`trajio` deliberately does not, the equivalent step has to happen somewhere, and the adapter
is the least surprising place. LSSD is a sum of squared distances and invariant under
translation, so no result changes. Note upstream centres `x` and `y` on their **mean** and
rebases `t` on its **first value**, where we use the first point for all three — again, both
are translations. Calling `feed()` directly still enforces the original guards, and there is a
test pinning that.

*Dropped:* cascade mode entirely (`feedIndex`, `cascadeRoot`, `batchDotsCascade*`), which is
a distinct variant from the paper rather than part of DOTS proper, and the `getAverageSED` /
`getMaxLSSD` diagnostics. All are straightforward to restore from the clone.

*Dropped from both dots-repo ports:* the Qt GUI (`mainwindow`, `qcustomplot`), the `.pro`
build files, `DotsException`, `Helper`, and `AlgorithmComparison`.

## Parameters in the original implementations

No knob was renamed or invented. Each `Params` key is the upstream parameter, with upstream's
default where it had one:

| Ours | Upstream signature | Type | Upstream default |
|---|---|---|---|
| `tol` | `simplify_douglas_peucker(first, last, tol, result)` | `value_type` (double here) | none — required |
| `buffer_size` | `SquishBatchSimplifier::simplifyByIndex(x, y, t, out, bufferSize)` | `int` | none — required, and `bufferSize > 4` |
| `lssd_threshold` | `DotsSimplifier::setParameters(lssdTh, k, maxVkSize)` | `double` | `10000.0` in the constructor |
| `k` | same call, 2nd argument | `double` | `2.0` |
| `max_vk_size` | same call, 3rd argument | `int` | `1e6` |

So DOTS is the only one of the three that had more than one knob upstream, and it had exactly
the three we expose. `batchDots` passes only `lssdThreshold` and lets `k` and `maxVkSize`
default, which is what our `run()` does too.

Two upstream knobs are deliberately not carried over. psimpl has a second entry point,
`simplify_douglas_peucker_n(first, last, count, result)`, which takes a **point count**
instead of a tolerance — the min-# form of the same algorithm; we ported the `tol` form only.
And DOTS' cascade mode adds `thStart` and `thStep` (`batchDotsCascadeByIndexOptions`), which
went with the rest of cascade mode.

### What DOTS' two secondary knobs actually do

Both control the DAG frontier rather than the error, and both trade search quality against
work done. Neither appears in `batchDots`, so upstream's own batch driver leaves them at
their defaults, and so does ours unless you set them.

**`k` = 2.0 — when to give up on a candidate predecessor.** It sets
`lssdUpperBound = lssdTh * k`. Each new point is tried against the frontier; a candidate
whose LSSD merely exceeds `lssdTh` is skipped for *this* point, but one that exceeds
`lssdTh * k` is marked *terminated* and dropped from the frontier permanently, on the
reasoning that it is now so far behind it will never be viable again:

```cpp
if (distance < lssdTh)            { /* link to it */ }
else if (distance > lssdUpperBound) { terminated[j] = true; ++numTerminated; }
```

Raising `k` keeps candidates alive longer — a more thorough search, more work per point, and
potentially fewer output points. Lowering it prunes sooner and is faster but can discard a
predecessor that would have led to a better path. It never affects correctness of the error
bound, only how close to optimal the result gets.

**`max_vk_size` = 1e6 — how large a frontier layer may grow before it is forced closed.**
`needUpdateVK()` is true when every frontier candidate has terminated *or* when the next
layer reaches `maxVkSize`:

```cpp
return (vK.count() == numTerminated || vL.count() >= maxVkSize);
```

Closing a layer runs `minimizeISSED`, swaps the layers, and — in the streaming path — runs
`viterbiDecode`, which is what actually emits output. So `max_vk_size` is a bound on memory
and on **output latency**: at the default of a million it effectively never fires, and layers
close only when the frontier dies out naturally. Setting it low forces frequent layer swaps,
which makes points come out sooner at the cost of a less optimal path.

### Upstream works in projected metres, not degrees

This is the part that makes those defaults make sense, and it is easy to miss. `Helper`
preprocesses every trajectory before any algorithm sees it — here in `parseMOPSI`:

```cpp
double timestamp = (double)QDateTime::fromString(...).toTime_t();   // epoch SECONDS
mercatorProject(longitude, latitude, x, y);                         // degrees -> metres
Helper::normalizeData(x, true);                                     // centre on the mean
Helper::normalizeData(y, true);
Helper::normalizeData(t, false);                                    // rebase on the first value
```

So upstream's inputs are mean-centred Mercator **metres** and seconds from zero. A `tol` of
10.0 is ten metres, and DOTS' guards — first `|t|` under a year in seconds, `|x|`, `|y|`
under 2 000 km — are sized for exactly that.

**We do none of that projection**, because `trajio` copies coordinates verbatim by design: a
Mopsi document is degrees, NGSIM is State Plane feet, MOT is pixels. Nothing in this port
projects or rescales. The only preprocessing is the translation `Dots::run` applies to get
past those guards, and even that differs slightly from upstream — upstream centres `x` and
`y` on their **mean** and rebases `t` on its **first value**; we translate all three to the
first point. Both are translations and LSSD is invariant under either.

The consequence: **upstream's default parameter values are meaningless on our data**, and a
tolerance in degrees is not a tolerance in metres. Either project the coordinates before
calling these, or sweep the parameter against the data — which is what the next section
measures.

## Parameter scales here

Measured on `mopsi-000003` (130 points, raw lon/lat degrees):

| | meaning | useful range on degrees |
|---|---|---|
| `tol` | perpendicular distance, in coordinate units | `1e-4` → 22 points, `1e-3` → 4 |
| `buffer_size` | a point budget, not an error bound | exact: 10 → 10, 50 → 50 |
| `lssd_threshold` | **squared** SED summed over points, in units²·points | `1e-12` → 126, `1e-9` → 84, `1e-7` → 21 |

`lssd_threshold` is quadratic in the coordinate unit, so it moves fastest of the three when
units change: upstream's default of `10000.0` in metres corresponds to nothing usable in
degrees, where `1e-3` already collapses a trajectory to its endpoints.

SQUISH cannot be given an error bound at all — upstream implements SQUISH, not SQUISH-E, so
there is no `mu`. Its knob is a point budget and nothing else.

## Verified against the originals

Each port was run side by side with the code it came from. Both upstreams were compiled
outside their build systems — psimpl is header-only, and `DotsSimplifier` /
`SquishBatchSimplifier` were built against a ~150-line shim for the Qt API they use
(`QVector`, `QMap`, `QPair`, `qSqrt`/`qRound`, `Q_OBJECT` as nothing, `DotsException::raise`
as a throw). `QVector<bool>` is a `std::deque<bool>` in the shim, because Qt's holds real
bools. The harness then compared outputs point for point over 40 real Mopsi trajectories and
8 synthetic ones (random walks at three scales, spirals, a staircase, a curve full of
repeated points), sweeping 8 tolerances, 7 buffer sizes and 8 LSSD thresholds, plus DOTS
through its online interface at every combination of `k` in {1.05, 2, 8} and `max_vk_size` in
{1, 2, 5, 1e6}, comparing not just the indices but *when* each was emitted, and psimpl in 3D.

**3807 comparisons, one mismatch**, and it is the documented one: on a unit staircase at
`tol = 0.5` every candidate point is exactly equidistant from the chord, so the key is
decided purely by how the projection fraction rounds. Re-running our DP with psimpl's `float`
fraction reproduces psimpl exactly, on that case and everywhere else — the port is otherwise
identical, and only inputs with exact distance ties can diverge.

Three checks go past matching upstream, to what the papers claim:

- `getLSSD`'s eight-prefix-sum closed form agrees with the brute-force sum of squared
  synchronous distances to 4e-9 relative over 200 random curves, so the O(1) form really is
  the quantity DOTS claims it is.
- Every shortcut DOTS emits is under its threshold — 0 violations over 60 runs on 20 Mopsi
  trajectories — and the output is **within 1.11× the true minimum** point count achievable
  under the same constraint, computed exhaustively over all legal shortcuts. That is the
  paper's near-optimality claim, measured.
- DP's error is the 1.88 `tol` above, not `tol`. It is the one place where fidelity to the
  reference costs a guarantee.

## Tests

`tests/test_simplify.cpp`, 19 cases. Every algorithm is checked against the shared contract —
endpoints kept, output a strictly increasing subsequence of valid indices, monotone in its
parameter — plus what is specific to each: DP flattening a straight line to two points and
running in 3D, SQUISH respecting its budget and rejecting `buffer_size <= 4`, DOTS accepting
epoch-ms through `indices` while `feed` still rejects it, and both SED algorithms refusing to
run without a clock.

Two cases cover the interface itself: a toy `Midpoints` algorithm that emits vertices found
nowhere in its input, confirming a non-subset algorithm needs only `run`/`simplify`; and a
check that for all three real algorithms `simplify` returns exactly the points `indices`
names.

One case is the regression net for the section above: *all reproduce the upstream reference
output* holds the indices psimpl, `SquishBatchSimplifier` and `DotsSimplifier` produce on a
fixed 40-point curve, at two parameter values each. Its coordinates are sixteenths so they
are exact doubles, and it fails the moment a port drifts.

Run: `./build/tests/ssk_tests simplify`
