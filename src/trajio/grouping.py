"""Turning a flat record stream into per-object groups.

NGSIM ships as one enormous CSV in which the rows of different vehicles are interleaved, so
"read the file, group by id" is the expensive step, not the parsing. Which strategy is right
depends entirely on the file's row order, and getting it wrong either shreds trajectories
into fragments or exhausts memory:

    contiguous   rows of an object are consecutive        O(1) memory   (pre-sorted files)
    buffered     interleaved, but the object count fits   O(objects)
    external     interleaved and too big for memory       O(chunk)
"""

from __future__ import annotations

import heapq
import itertools
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
    """A new group whenever the key changes. Constant memory; wrong on interleaved input."""
    for k, batch in itertools.groupby(records, key):
        yield k, list(batch)


def group_buffered(
    records: Iterable[R], key: Callable[[R], str]
) -> Iterator[tuple[str, list[R]]]:
    """Accumulate every object in a dict; yield groups at end of stream."""
    groups: dict[str, list[R]] = {}
    for rec in records:
        groups.setdefault(key(rec), []).append(rec)
    yield from groups.items()


def group_external(
    records: Iterable[R],
    key: Callable[[R], str],
    *,
    sort_key: Callable[[R], Any] | None = None,
    chunk_rows: int = 1_000_000,
) -> Iterator[tuple[str, list[R]]]:
    """External merge sort on `(key, sort_key)`, then a contiguous group.

    Memory is `chunk_rows` records; disk is roughly the size of the pickled stream.
    """
    sort_key = sort_key or (lambda r: 0)

    def spill(buf: list[R], path: str) -> None:
        buf.sort(key=lambda r: (key(r), sort_key(r)))
        with open(path, "wb") as fh:
            for rec in buf:
                pickle.dump((key(rec), sort_key(rec), rec), fh, protocol=pickle.HIGHEST_PROTOCOL)

    def read_run(path: str) -> Iterator[tuple[str, Any, R]]:
        with open(path, "rb") as fh:
            while True:
                try:
                    yield pickle.load(fh)
                except EOFError:
                    return

    # A consumer that stops early leaves runs open, and Windows will not delete an open file.
    with tempfile.TemporaryDirectory(prefix="trajio-sort-", ignore_cleanup_errors=True) as tmp:
        runs: list[str] = []
        chunks = iter(records)
        while buf := list(itertools.islice(chunks, chunk_rows)):
            runs.append(os.path.join(tmp, f"run-{len(runs):05d}.pkl"))
            spill(buf, runs[-1])

        merged = heapq.merge(*(read_run(p) for p in runs), key=lambda item: (item[0], item[1]))
        for k, batch in itertools.groupby(merged, lambda item: item[0]):
            yield k, [rec for _k, _sk, rec in batch]


def group_by(
    strategy: str,
    records: Iterable[R],
    key: Callable[[R], str],
    *,
    sort_key: Callable[[R], Any] | None = None,
    chunk_rows: int = 1_000_000,
) -> Iterator[tuple[str, list[R]]]:
    """Dispatch on `"contiguous"` / `"buffered"` / `"external"`."""
    if strategy == "contiguous":
        return group_contiguous(records, key)
    if strategy == "buffered":
        return group_buffered(records, key)
    if strategy == "external":
        return group_external(records, key, sort_key=sort_key, chunk_rows=chunk_rows)
    raise ValueError(f"unknown grouping strategy {strategy!r}")
