# References

Everything this project builds on: the papers, the code it ports, and the datasets it runs on.
Each entry says what was taken and where the details are.

## The algorithm this project implements

**Cheng, Huang & Jiang.** *Simplification of Trajectory Streams.* SoCG 2025.
The streaming δ-simplification under the continuous Fréchet distance that this project
implements, with the grid-cell redesign agreed at the 2026-08-08 supervisor meeting. The
earlier prototype of it is in `archive/`.

## Ported code

Each of these was cloned, built outside its own build system, and run side by side with our
port. None is vendored — the ported files are ours, the algorithms are theirs. Deviations are
recorded per algorithm.

| What | Source | Licence | Ported into | Documented in |
|---|---|---|---|---|
| Douglas–Peucker (DPn) | [psimpl](https://psimpl.sourceforge.net/) v7 | MPL 1.1 | `simplify/douglas_peucker.hpp` | [simplification.md](simplification.md) |
| SQUISH | [caoweiquan322/dots](https://github.com/caoweiquan322/dots) | see repo | `simplify/squish.{hpp,cpp}` | [simplification.md](simplification.md) |
| DOTS | [caoweiquan322/dots](https://github.com/caoweiquan322/dots) | see repo | `simplify/dots.{hpp,cpp}` | [simplification.md](simplification.md) |
| Continuous Fréchet distance | [mwernerds/frechetrange](https://github.com/mwernerds/frechetrange), `detail/dv` | MIT | `frechet/frechet.hpp` | [frechet-distance.md](frechet-distance.md) |

### The papers behind them

**Douglas & Peucker.** *Algorithms for the reduction of the number of points required to
represent a digitized line or its caricature.* Cartographica 10(2), 1973. The `psimpl` port is
of the **point-count variant (DPn)**, which takes a vertex budget rather than a tolerance.

**Muckell, Hwang, Patil, Lawson, Ping & Ravi.** *SQUISH: an online approach for GPS trajectory
compression.* COM.Geo 2011. The dots repository implements plain SQUISH, not SQUISH-E, so
there is no error parameter — its knob is a point budget.

**Cao & Li.** *DOTS: An online and near-optimal trajectory simplification algorithm.* Journal
of Systems and Software 126, 2017. Note the defect found in the reference implementation,
recorded in `archive/[2026-08-16] Incident - DOTS emits non-monotone indices.md`.

**Alt & Godau.** *Computing the Fréchet distance between two polygonal curves.* International
Journal of Computational Geometry & Applications 5, 1995. The free-space diagram, which is
what every Fréchet implementation here computes.

**Dütsch & Vahrenhold.** *A filter-and-refinement algorithm for range queries based on the
Fréchet distance.* SIGSPATIAL 2017 — the **GIS Cup 2017 winning entry**, and the code we port.

**Baldus & Bringmann.** *A fast implementation of near neighbors queries for Fréchet
distance.* SIGSPATIAL 2017 — the GIS Cup runner-up, used as an **independent cross-check**
rather than ported.

**Werner, M.** [`frechetrange`](https://github.com/mwernerds/frechetrange) consolidates the
three GIS Cup entries into one dependency-free form, with Fabian Dütsch credited for adapting
the implementations. Our port is of that consolidation's `dv` header.

## Datasets

Sizes, counts and everything about parsing are in [datasets.md](datasets.md); this is the
attribution.

**GeoLife GPS Trajectories** — Microsoft Research Asia, 2007–2012, mostly Beijing. Research
use licence bundled in the archive. Cite:

- Zheng, Zhang, Xie & Ma. *Mining interesting locations and travel sequences from GPS
  trajectories.* WWW 2009.
- Zheng, Li, Chen, Xie & Ma. *Understanding mobility based on GPS data.* UbiComp 2008.
- Zheng, Xie & Ma. *GeoLife: A collaborative social networking service among user, location
  and trajectory.* IEEE Data Engineering Bulletin 33(2), 2010.

**Mopsi Routes 2014** — SIPU group, University of Eastern Finland, 2008–2014, around Joensuu.
Citation required by the provider:

- Mariescu-Istodor & Fränti. *Grid-based method for GPS route analysis for retrieval.*
  ACM Transactions on Spatial Algorithms and Systems 3(3), 2017.

**NGSIM** — Next Generation Simulation vehicle trajectories, US Department of Transportation /
FHWA, public domain, via `data.transportation.gov`. Four sites: us-101, i-80, lankershim,
peachtree. Note the epoch discrepancy documented in [datasets.md](datasets.md) — peachtree's
timestamps are not on the Unix epoch.

**MOT video object tracks** — MOT17, MOT20 and DanceTrack bounding-box annotations, taken from
the [Trajectory Simplify Benchmark](https://github.com/foreverYoungGitHub/trajectory-simplify-benchmark)
distribution because `motchallenge.net` was unreachable. Cite the underlying tracking dataset
used:

- Milan, Leal-Taixé, Reid, Roth & Schindler. *MOT16: A benchmark for multi-object tracking.*
  arXiv:1603.00831, 2016 (MOT17 shares its sequences).
- Dendorfer et al. *MOT20: A benchmark for multi object tracking in crowded scenes.*
  arXiv:2003.09003, 2020.
- Sun et al. *DanceTrack: Multi-object tracking in uniform appearance and diverse motion.*
  CVPR 2022.

## Tools

Python `matplotlib` and `numpy` for `src/viz/` only — pinned in `requirements.txt`. Nothing
under `src/trajio/` or `src/algo/` depends on anything outside the standard library, which is
what keeps ingestion portable and the C++ build offline.
