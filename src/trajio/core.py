"""Core types.

A dataset source yields `Track`s. A track is a list of points, each point a tuple of
coordinate **strings copied verbatim from the source file** -- no parsing, no rounding, no
unit conversion. `z` is `None` where the dataset has no altitude reading.

The writer decides 2D or 3D per track: with `--dims 3`, a track whose points all carry a `z`
is written `"dim": 3`; any missing `z` and it is written `"dim": 2`.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any, Callable, NamedTuple

__all__ = [
    "Point",
    "Track",
    "TrajectorySource",
    "Option",
    "register_source",
    "get_source",
    "source_names",
    "fmt",
]

#: (x, y, z) as source text; z is None when unavailable.
Point = tuple[str, str, str | None]


@dataclass(slots=True)
class Track:
    """One trajectory: an id for provenance, and its points in file order."""

    id: str
    points: list[Point] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.points)

    def has_z(self) -> bool:
        return bool(self.points) and all(p[2] is not None for p in self.points)


def fmt(value: float) -> str:
    """Format a *computed* number (only MOT box centres) without trailing zeros."""
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text or "0"


class Option(NamedTuple):
    """A source option, exposed on the command line as `--opt name=value`."""

    type: Callable[[str], Any]
    default: Any
    help: str


class TrajectorySource(ABC):
    """One dataset's on-disk format."""

    #: registry key, also the CLI `--source` value and the output filename prefix
    name: str = ""
    #: one line for `trajio sources`
    summary: str = ""
    #: dimensionalities this dataset can supply: (2,) or (2, 3)
    dims: tuple[int, ...] = (2,)
    #: what x and y are, for the docs
    axes: str = "x, y"
    #: constructor options beyond `root`
    OPTIONS: dict[str, Option] = {}

    def __init__(self, root: str) -> None:
        self.root = root

    @abstractmethod
    def tracks(self) -> Iterator[Track]:
        """Yield one `Track` per trajectory, points in file order."""

    def __iter__(self) -> Iterator[Track]:
        return self.tracks()

    def __repr__(self) -> str:
        return f"{type(self).__name__}(root={self.root!r})"


_REGISTRY: dict[str, type[TrajectorySource]] = {}


def register_source(cls: type[TrajectorySource]) -> type[TrajectorySource]:
    if not cls.name:
        raise ValueError(f"{cls.__name__} must set a `name`")
    _REGISTRY[cls.name] = cls
    return cls


def get_source(name: str) -> type[TrajectorySource]:
    try:
        return _REGISTRY[name]
    except KeyError:
        raise KeyError(f"unknown source {name!r}; known: {', '.join(source_names())}") from None


def source_names() -> list[str]:
    return sorted(_REGISTRY)
