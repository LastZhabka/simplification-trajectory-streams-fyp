"""Turning a flat record stream into per-object groups.

Four of the six datasets ship as one enormous CSV in which the rows of different objects are
interleaved, so "read the file, group by id" is the expensive step, not the parsing. Which
strategy is right depends entirely on the file's row order, and getting it wrong either
shreds trajectories into fragments or exhausts memory:

    contiguous   rows of an object are consecutive        O(1) memory   (pre-sorted files)
    buffered     interleaved, but the object count fits   O(objects)    (NGSIM, one AIS day
                                                                         subset, Porto)
    external     interleaved and too big for memory       O(chunk)      (a full AIS day/year)

Each source picks a sensible default and exposes `grouping=` so it can be overridden.
"""

from __future__ import annotations

import heapq
import os
import pickle
import tempfile
from collections.abc import Callable, Iterable, Iterator
from typing import Any, TypeVar

__all__ = ["group_contiguous", "group_buffered", "group_external", "group_by"]

R = TypeVar("R")


def group_contiguous(
    records: Iterable[R], key: Callable[[R], str]
) -> Iterator[tuple[str, list[R]]]:
    """Flush whenever the key changes. Constant memory; wrong on interleaved input.

    Correct for per-file datasets (GeoLife, Mopsi, MOT) and for CSVs you have already sorted
    by object id.
    """
    current: str | None = None
    batch: list[R] = []
    for rec in records:
        k = key(rec)
        if k != current:
            if batch:
                yield current, batch  # type: ignore[misc]
            current, batch = k, []
        batch.append(rec)
    if batch:
        yield current, batch  # type: ignore[misc]


def group_buffered(
    records: Iterable[R],
    key: Callable[[R], str],
    *,
    idle_rows: int | None = None,
) -> Iterator[tuple[str, list[R]]]:
    """Accumulate every object in a dict; yield groups at end of stream.

    `idle_rows` bounds memory on a file that is *mostly* grouped: an object whose key has not
    been seen for that many rows is assumed finished and flushed early. Leave it `None` when
    correctness matters more than memory -- an object that reappears after being flushed
    comes out as two trajectories, which is a silent data change, not an error.
    """
    groups: dict[str, list[R]] = {}
    last_seen: dict[str, int] = {}
    for i, rec in enumerate(records):
        k = key(rec)
        groups.setdefault(k, []).append(rec)
        if idle_rows is None:
            continue
        last_seen[k] = i
        if i % max(1, idle_rows // 4) == 0:
            cutoff = i - idle_rows
            stale = [g for g, seen in last_seen.items() if seen < cutoff]
            for g in stale:
                yield g, groups.pop(g)
                del last_seen[g]
    yield from groups.items()


def group_external(
    records: Iterable[R],
    key: Callable[[R], str],
    *,
    sort_key: Callable[[R], Any] | None = None,
    chunk_rows: int = 1_000_000,
    tmpdir: str | None = None,
) -> Iterator[tuple[str, list[R]]]:
    """External merge sort on `(key, sort_key)`, then a contiguous group.

    Memory is `chunk_rows` records; disk is roughly the size of the pickled stream. This is
    the path for a full AIS day (~10 M rows) or anything larger, and it is why the reader
    never assumes it can hold a dataset in RAM.
    """
    sort_key = sort_key or (lambda r: 0)
    runs: list[str] = []
    tmp = tempfile.mkdtemp(prefix="trajio-sort-", dir=tmpdir)

    def spill(buf: list[R]) -> None:
        buf.sort(key=lambda r: (key(r), sort_key(r)))
        path = os.path.join(tmp, f"run-{len(runs):05d}.pkl")
        with open(path, "wb") as fh:
            for rec in buf:
                pickle.dump((key(rec), sort_key(rec), rec), fh, protocol=pickle.HIGHEST_PROTOCOL)
        runs.append(path)

    def read_run(path: str) -> Iterator[tuple[str, Any, R]]:
        with open(path, "rb") as fh:
            while True:
                try:
                    yield pickle.load(fh)
                except EOFError:
                    return

    try:
        buf: list[R] = []
        for rec in records:
            buf.append(rec)
            if len(buf) >= chunk_rows:
                spill(buf)
                buf = []
        if buf:
            spill(buf)

        merged = heapq.merge(*(read_run(p) for p in runs), key=lambda item: (item[0], item[1]))
        current: str | None = None
        batch: list[R] = []
        for k, _sk, rec in merged:
            if k != current:
                if batch:
                    yield current, batch  # type: ignore[misc]
                current, batch = k, []
            batch.append(rec)
        if batch:
            yield current, batch  # type: ignore[misc]
    finally:
        for path in runs:
            try:
                os.remove(path)
            except OSError:
                pass
        try:
            os.rmdir(tmp)
        except OSError:
            pass


def group_by(
    strategy: str,
    records: Iterable[R],
    key: Callable[[R], str],
    *,
    sort_key: Callable[[R], Any] | None = None,
    idle_rows: int | None = None,
    chunk_rows: int = 1_000_000,
) -> Iterator[tuple[str, list[R]]]:
    """Dispatch on `"contiguous"` / `"buffered"` / `"external"`."""
    if strategy == "contiguous":
        return group_contiguous(records, key)
    if strategy == "buffered":
        return group_buffered(records, key, idle_rows=idle_rows)
    if strategy == "external":
        return group_external(records, key, sort_key=sort_key, chunk_rows=chunk_rows)
    raise ValueError(f"unknown grouping strategy {strategy!r}")
