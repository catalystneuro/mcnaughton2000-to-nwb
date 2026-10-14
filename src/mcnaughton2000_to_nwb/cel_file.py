"""Parse McNaughton lab Xclust ASCII .CEL spike files."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# The tracker writes (0, 255) when it loses the animal.
LOST_POSITION = (0.0, 255.0)


@dataclass
class CelFile:
    """Parsed contents of a single .CEL file."""

    path: Path
    fields: list[str]
    header: dict[str, str]
    data: pd.DataFrame

    # Derived
    cluster: int | None = field(init=False)
    start_time_sec: float = field(init=False)
    end_time_sec: float = field(init=False)
    has_position: bool = field(init=False)

    def __post_init__(self):
        self.cluster = self._parse_cluster()
        self.start_time_sec = self._parse_time_str(self.header.get("Start time", ""))
        self.end_time_sec = self._parse_time_str(self.header.get("End time", ""))
        self.has_position = "pos_x" in self.fields and "pos_y" in self.fields

    @property
    def spike_times(self) -> np.ndarray:
        """Spike times in seconds."""
        return pd.to_numeric(self.data["time"], errors="coerce").to_numpy(dtype=np.float64)

    @property
    def position(self) -> np.ndarray | None:
        """Position (pos_x, pos_y) in pixels at each spike, NaN where the tracker lost the animal."""
        if not self.has_position:
            return None
        xy = np.column_stack([
            pd.to_numeric(self.data["pos_x"], errors="coerce").to_numpy(dtype=np.float64),
            pd.to_numeric(self.data["pos_y"], errors="coerce").to_numpy(dtype=np.float64),
        ])
        lost = (xy[:, 0] == LOST_POSITION[0]) & (xy[:, 1] == LOST_POSITION[1])
        xy[lost] = np.nan
        return xy

    def _parse_cluster(self) -> int | None:
        raw = self.header.get("Cluster", "")
        if not raw:
            return None
        m = re.search(r"(\d+)", raw)
        return int(m.group(1)) if m else None

    @staticmethod
    def _parse_time_str(time_str: str) -> float:
        """Parse 'H:MM:SS' or 'MM:SS' to seconds."""
        if not time_str:
            return float("nan")
        parts = time_str.strip().split(":")
        try:
            parts = [int(p) for p in parts]
        except ValueError:
            return float("nan")
        if len(parts) == 3:
            return parts[0] * 3600 + parts[1] * 60 + parts[2]
        if len(parts) == 2:
            return parts[0] * 60 + parts[1]
        return float("nan")


def read_cel_file(path: Path) -> CelFile:
    """Read and parse a single .CEL file."""
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    fields: list[str] | None = None
    header_kv: dict[str, str] = {}
    endheader_idx: int | None = None

    for i, line in enumerate(lines[:800]):
        s = line.strip()

        if s.startswith("%") and ":" in s and not re.match(r"^%?\s*fields\s*:", s, flags=re.IGNORECASE):
            kv = s.lstrip("%").strip()
            k, v = kv.split(":", 1)
            header_kv[k.strip()] = v.strip()

        if re.match(r"^%?\s*fields\s*:", s, flags=re.IGNORECASE):
            rhs = s.split(":", 1)[1].strip()
            fields = rhs.split()

        if s == "%%ENDHEADER":
            endheader_idx = i
            break

    if fields is None:
        raise ValueError(f"No Fields line found in {path}")
    if endheader_idx is None:
        raise ValueError(f"No %%ENDHEADER found in {path}")

    data_start = endheader_idx + 1

    df = pd.read_csv(
        path,
        sep=r"\s+",
        engine="python",
        skiprows=data_start,
        names=fields,
        encoding_errors="ignore",
    ).dropna(how="all")

    return CelFile(path=path, fields=fields, header=header_kv, data=df)
