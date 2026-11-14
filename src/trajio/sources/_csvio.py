"""Streaming rows out of one CSV, or several."""

from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path

__all__ = ["iter_csv_files", "read_dict_rows", "normalise_header"]


def iter_csv_files(root: Path, pattern: str) -> list[Path]:
    """A single file, or every match of `pattern` under a directory, in name order."""
    if root.is_file():
        return [root]
    files = sorted(p for p in root.glob(pattern) if p.is_file())
    if not files:
        raise FileNotFoundError(f"no files matching {pattern!r} under {root}")
    return files


def normalise_header(names: list[str]) -> dict[str, int]:
    """Case- and space-insensitive header index. NGSIM alone ships several spellings."""
    return {name.strip().lower(): i for i, name in enumerate(names)}


def read_dict_rows(paths: list[Path]) -> Iterator[tuple[Path, dict[str, int], list[str]]]:
    """Stream `(file, header_index, row)` across several CSVs, skipping each header line."""
    for path in paths:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh)
            header = normalise_header(next(reader, []))
            for row in reader:
                if row:
                    yield path, header, row
