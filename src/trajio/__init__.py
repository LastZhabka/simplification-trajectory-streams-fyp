"""trajio -- read public trajectory datasets, write one JSON trajectory per file.

Values are copied from the source files verbatim: no projection, no unit conversion, no
filtering, no splitting. Standard library only.

    from trajio import open_source, export

    src = open_source("mopsi", "data/downloads/mopsi")
    export(src, "data/trajectories/mopsi", dims=2)
"""

from __future__ import annotations

from typing import Any

from .core import (
    Option,
    Point,
    Track,
    TrajectorySource,
    get_source,
    register_source,
    source_names,
)
from .export import export
from .stats import summarise

from . import sources as _sources  # noqa: F401  -- import registers the sources

__all__ = [
    "Point",
    "Track",
    "TrajectorySource",
    "Option",
    "open_source",
    "export",
    "summarise",
    "source_names",
    "get_source",
    "register_source",
]

__version__ = "0.2.0"


def open_source(name: str, root: str, **options: Any) -> TrajectorySource:
    """`open_source("geolife", "raw/geolife")`."""
    return get_source(name)(root, **options)
