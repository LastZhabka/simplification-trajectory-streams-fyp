"""    python -m trajio sources
    python -m trajio describe geolife
    python -m trajio stats    --source mopsi --root data/downloads/mopsi
    python -m trajio export   --source mopsi --root data/downloads/mopsi --dims 2 \
                              --out data/trajectories/mopsi
    python -m trajio selftest
"""

from __future__ import annotations

import argparse
import sys
from typing import Any

from .core import get_source, source_names
from .export import export
from .stats import summarise


def _source_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--source", required=True, choices=source_names())
    p.add_argument("--root", required=True, help="directory (or file) holding the raw data")
    p.add_argument("--opt", action="append", default=[], metavar="KEY=VALUE",
                   help="source-specific option; repeatable. See `describe <source>`")
    p.add_argument("--limit", type=int, default=None, help="stop after N trajectories")


def _build(args: argparse.Namespace) -> Any:
    cls = get_source(args.source)
    options: dict[str, Any] = {}
    for pair in args.opt:
        key, sep, value = pair.partition("=")
        if not sep:
            raise SystemExit(f"--opt expects KEY=VALUE, got {pair!r}")
        spec = cls.OPTIONS.get(key.strip())
        if spec is None:
            known = ", ".join(sorted(cls.OPTIONS)) or "(none)"
            raise SystemExit(f"{args.source}: unknown option {key!r}; known: {known}")
        options[key.strip()] = spec.type(value)
    return cls(args.root, **options)


def cmd_sources(_args: argparse.Namespace) -> int:
    for name in source_names():
        cls = get_source(name)
        print(f"{name:9s} dims={','.join(map(str, cls.dims)):4s} {cls.summary}")
    return 0


def cmd_describe(args: argparse.Namespace) -> int:
    cls = get_source(args.name)
    print(f"{cls.name}: {cls.summary}")
    print(f"  dims  : {', '.join(map(str, cls.dims))}")
    print(f"  axes  : {cls.axes}")
    for key, spec in cls.OPTIONS.items():
        print(f"  --opt {key}={spec.default!r}  {spec.help}")
    print("\n" + "\n".join("  " + ln for ln in (cls.__doc__ or "").strip().splitlines()))
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    print(summarise(_build(args), limit=args.limit))
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    source = _build(args)
    index = export(source, args.out, dims=args.dims, limit=args.limit)
    n3 = sum(1 for row in index if row[3] == 3)
    print(f"wrote {len(index)} files to {args.out}  ({sum(r[2] for r in index):,} points)")
    if args.dims == 3:
        print(f"  3D: {n3}   2D (no altitude on some points): {len(index) - n3}")
    return 0 if index else 1


def cmd_selftest(_args: argparse.Namespace) -> int:
    from ._selftest import run

    return run()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="trajio", description="trajectory datasets -> JSON")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("sources", help="list datasets").set_defaults(func=cmd_sources)

    d = sub.add_parser("describe", help="show one dataset's options")
    d.add_argument("name", choices=source_names())
    d.set_defaults(func=cmd_describe)

    s = sub.add_parser("stats", help="count trajectories and points")
    _source_args(s)
    s.set_defaults(func=cmd_stats)

    e = sub.add_parser("export", help="write one JSON file per trajectory")
    _source_args(e)
    e.add_argument("--out", required=True, help="output directory")
    e.add_argument("--dims", type=int, default=2, choices=[2, 3])
    e.set_defaults(func=cmd_export)

    sub.add_parser("selftest", help="run every parser against fixtures").set_defaults(
        func=cmd_selftest
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
