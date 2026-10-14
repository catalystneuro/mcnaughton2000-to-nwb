"""Parse McNaughton lab .RMA binary rate map files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


RMA_FILE_SIZE = 32768
MAP_SHAPE = (64, 64)
N_PIXELS = MAP_SHAPE[0] * MAP_SHAPE[1]  # 4096

# The original analysis writes this value in bins the animal never visited.
UNVISITED = -999.0


@dataclass
class RmaFile:
    """Parsed contents of a single .RMA file."""

    path: Path
    rate_map: np.ndarray  # (64, 64) float32, NaN in unvisited bins
    occupancy_map: np.ndarray  # (64, 64) int32, 60 Hz video samples per bin


def read_rma_file(path: Path) -> RmaFile:
    """Read and parse a single .RMA file.

    Format: 32768 bytes total
      - First 4096 float32 big-endian → rate map 64×64
      - Next 4096 int32 big-endian → occupancy map 64×64
    """
    data = path.read_bytes()
    if len(data) != RMA_FILE_SIZE:
        raise ValueError(f"{path}: expected {RMA_FILE_SIZE} bytes, got {len(data)}")

    rate_map = (
        np.frombuffer(data, dtype=">f4", count=N_PIXELS, offset=0)
        .reshape(MAP_SHAPE)
        .astype(np.float32)
    )
    rate_map[rate_map == UNVISITED] = np.nan
    occupancy_map = (
        np.frombuffer(data, dtype=">i4", count=N_PIXELS, offset=N_PIXELS * 4)
        .reshape(MAP_SHAPE)
        .astype(np.int32)
    )

    return RmaFile(path=path, rate_map=rate_map, occupancy_map=occupancy_map)
