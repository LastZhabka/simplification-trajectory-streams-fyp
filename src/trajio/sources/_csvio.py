"""Reading a possibly-compressed, possibly-sharded CSV without unpacking it first.

The big datasets ship as `.zip` (AIS daily, Porto `train.csv.zip`), `.gz`, or `.csv.zst`
(newer AIS), and each is tens of GB when expanded. Everything here streams.
"""

from __future__ import annotations

import csv
import gzip
import io
import sys
import zipfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO

__all__ = ["iter_csv_files", "open_csv_text", "read_dict_rows", "normalise_header"]

# POLYLINE strings in the Porto file are long; the csv default limit is stingy.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


@contextmanager
def open_csv_text(path: Path, *, member: str | None = None) -> Iterator[IO[str]]:
    """Yield a text handle for `.csv`, `.csv.gz`, `.zip` (one CSV member) or `.csv.zst`."""
    suffix = path.suffix.lower()
    if suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = [n for n in zf.namelist() if n.lower().endswith(".csv")]
            if member is not None:
                names = [n for n in zf.namelist() if n == member]
            if not names:
                raise FileNotFoundError(f"no CSV member in {path}")
            with zf.open(names[0]) as binary:
                yield io.TextIOWrapper(binary, encoding="utf-8", errors="replace", newline="")
    elif suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
            yield fh
    elif suffix == ".zst":
        try:  # stdlib gained zstd in 3.14; before that this needs the pypi package
            from compression import zstd  # type: ignore[import-not-found]

            with zstd.open(path, "rt", encoding="utf-8", newline="") as fh:  # type: ignore[attr-defined]
                yield fh
        except ImportError:
            try:
                import zstandard  # type: ignore[import-not-found]
            except ImportError:
                raise RuntimeError(
                    f"{path.name} is zstd-compressed; install `zstandard` (pip install "
                    f"zstandard), use Python >= 3.14, or decompress with `zstd -d` first"
                ) from None
            with open(path, "rb") as raw:
                reader = zstandard.ZstdDecompressor().stream_reader(raw)
                yield io.TextIOWrapper(reader, encoding="utf-8", errors="replace", newline="")
    else:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
            yield fh


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


def read_dict_rows(
    paths: list[Path], *, member: str | None = None
) -> Iterator[tuple[Path, dict[str, int], list[str]]]:
    """Stream `(file, header_index, row)` across several CSVs, skipping each header line."""
    for path in paths:
        with open_csv_text(path, member=member) as fh:
            reader = csv.reader(fh)
            try:
                header = normalise_header(next(reader))
            except StopIteration:
                continue
            for row in reader:
                if row:
                    yield path, header, row
