"""Dataset sources. Importing this package registers every one of them."""

from __future__ import annotations

from .geolife import GeoLife
from .mopsi import Mopsi
from .mot import MOTTracks
from .ngsim import NGSIM

__all__ = ["GeoLife", "Mopsi", "NGSIM", "MOTTracks"]
