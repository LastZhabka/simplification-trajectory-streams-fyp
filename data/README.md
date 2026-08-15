# data/

Everything the project reads or writes that is not source. All four subdirectories are
git-ignored — the downloads are gigabytes and mostly redistributable only under their own
licences, and the rest is regenerable from them.

```
downloads/       the datasets exactly as fetched, in their own formats   -- see its README
trajectories/    those datasets converted to our format, by trajio
synthetic/       generated test curves, by src/viz/generate.py
renders/         images, by src/viz/plot.py
```

`downloads/` is the only one that cannot be reproduced by running something locally, and
the only one with content worth reading before use. Everything else is output.

## Our format

`trajectories/` and `synthetic/` hold the same thing, and it is the only trajectory format
in this repository: one JSON document per curve.

```json
{
  "dim": 2,
  "name": "geolife/000/20081023025304",
  "points": [
    [116.318417, 39.984702],
    [116.31845, 39.984683]
  ]
}
```

- `dim` is 2 or 3 and every point has exactly that many coordinates.
- `name` is provenance — for a converted dataset, the record the trajectory came from; for
  a generated curve, the generator. Optional.
- Coordinates converted from a dataset are the **source text**, placed straight into the
  JSON number slot, so no float round-trip happens between the dataset file and ours.

Read it with `ssk::io::read_trajectory_file` (C++, `src/algo/io/trajectory.hpp`) or
`trajectory_io.load` (Python, `src/viz/`).

A *result* document — what a simplification run writes — is the same document plus
`algorithm`, `mode`, `params`, `stats`, `input_points` (the original curve, with `points`
now the simplified one) and, when checked, `frechet`. `src/viz/plot.py` renders either.

## Refilling it

From the repository root, with the virtual environment active and `src/` on `PYTHONPATH` —
see the setup section of [`../README.md`](../README.md):

```powershell
$env:PYTHONPATH = "$PWD\src"

# downloads/ -> trajectories/
python -m trajio export --source mopsi --root data\downloads\mopsi --dims 2 `
                        --out data\trajectories\mopsi

# synthetic/ and renders/
python src\viz\generate.py spiral2d -n 600 -o data\synthetic\spiral2d.json
python src\viz\plot.py data\synthetic\spiral2d.json -o data\renders\spiral2d.png
```

Download commands for `downloads/` are in [`downloads/README.md`](downloads/README.md) and
[`../docs/datasets.md`](../docs/datasets.md).

Each export directory also gets an `index.json`:

```json
{"dataset": "mopsi", "count": 6779,
 "trajectories": [{"file": "mopsi-000001.json", "name": "...", "points": 412, "dim": 2}]}
```

`name` there is what traces a file back to the record it came from, which is what makes a
result reproducible without re-running the export.
