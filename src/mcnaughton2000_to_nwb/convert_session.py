"""Convert one subject-session of the McNaughton Neurolab dataset to NWB.

Each subject folder (e.g. FD4RAT1) represents one continuous recording session
with interleaved task epochs (BL=Baseline, ES=Escher Staircase, MC=Magic Carpet).
The CEL files contain sorted spike data from Xclust with timestamps in seconds
relative to recording start. Flight sessions also include position (pos_x, pos_y)
at each spike time.

The archived files carry Windows 8.3 short names; the original names
(BL2.cell5, ES.cell1.rmap), which give the epoch and the cluster, are recovered
from each tetrode folder's WS_FTP.LOG (see short_names.py).
"""

from __future__ import annotations

import hashlib
import math
import uuid
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from hdmf.common import DynamicTable, DynamicTableRegion, VectorData
from pynwb import NWBFile, NWBHDF5IO
from pynwb.behavior import Position
from pynwb.ecephys import ElectrodeGroup
from pynwb.file import Subject
from pynwb.misc import Units

from .cel_file import CelFile, read_cel_file
from .metadata import (
    EXPERIMENT_DESCRIPTION,
    INSTITUTION,
    LAB,
    RELATED_PUBLICATIONS,
    SEX,
    SESSION_TYPES,
    SHARED_RECORDING,
    SPECIES,
    STRAIN,
    get_session_metadata,
)
from .rma_file import read_rma_file
from .short_names import parse_long_name, read_long_names

TASK_SESSION_TYPES = ("ES", "MC")


@dataclass
class SpikeFile:
    """One original spike file: one cluster on one tetrode within one epoch."""

    tetrode: str
    long_name: str
    epoch_name: str
    cel: CelFile

    @property
    def in_window(self) -> np.ndarray:
        """Rows whose spike falls within the file's epoch window.

        In 48 of the original files the first row is a spike from before the
        window, which also appears in the cluster's file for the previous epoch.
        """
        t = self.cel.spike_times
        return (t >= self.cel.start_time_sec) & (t <= self.cel.end_time_sec)


def _tetrode_dirs(root: Path) -> list[Path]:
    return sorted(
        [p for p in root.iterdir() if p.is_dir() and p.name.upper().startswith("TT")],
        key=lambda p: int(p.name[2:]),
    )


def _read_spike_files(raw_root: Path, stub_test: bool = False) -> list[SpikeFile]:
    spike_files = []
    for tt_dir in _tetrode_dirs(raw_root):
        long_names = read_long_names(tt_dir)
        cel_paths = sorted(p for p in tt_dir.iterdir() if p.is_file() and p.suffix.upper() == ".CEL")
        if stub_test:
            cel_paths = cel_paths[:2]
        for cel_path in cel_paths:
            long_name = long_names[cel_path.name]
            epoch_name, cell = parse_long_name(long_name)
            cel = read_cel_file(cel_path)
            if cel.cluster != cell:
                raise ValueError(f"{cel_path}: cluster {cel.cluster} in header, cell {cell} in name {long_name}")
            spike_files.append(SpikeFile(tt_dir.name, long_name, epoch_name, cel))
    return spike_files


def convert_session(
    base_dir: Path,
    subject_name: str,
    output_dir: Path,
    stub_test: bool = False,
) -> Path:
    """Convert one subject-session to NWB.

    Parameters
    ----------
    base_dir : Path
        Root of the McNaughton_Neurolab directory containing
        1_RAW_(original_files)/ and ANALYZED_(original_files)/.
    subject_name : str
        Subject folder name (e.g. "FD4RAT1", "PREFLI~1").
    output_dir : Path
        Directory for the output .nwb file.
    stub_test : bool
        If True, only process a small subset of data for quick testing.

    Returns
    -------
    Path
        Path to the written NWB file.
    """
    raw_root = base_dir / "1_RAW_(original_files)" / subject_name
    analyzed_root = base_dir / "ANALYZED_(original_files)" / subject_name

    if not raw_root.exists():
        raise FileNotFoundError(f"RAW directory not found: {raw_root}")

    meta = get_session_metadata(subject_name)

    # --- Create NWB file ---
    nwbfile = NWBFile(
        session_description=meta.session_description,
        identifier=str(uuid.uuid4()),
        session_start_time=meta.session_date,
        experiment_description=EXPERIMENT_DESCRIPTION,
        institution=INSTITUTION,
        lab=LAB,
        related_publications=RELATED_PUBLICATIONS,
        experimenter=["Knierim, James J.", "McNaughton, Bruce L.", "Poe, Gina R."],
        keywords=[
            "hippocampus",
            "place cells",
            "tetrode",
            "spaceflight",
            "microgravity",
            "spatial navigation",
            "Neurolab",
            "STS-90",
        ],
        subject=Subject(
            subject_id=meta.rat_id,
            species=SPECIES,
            strain=STRAIN,
            sex=SEX,
            description=(
                f"Adult male Fischer 344 rat ({meta.rat_id}) from Taconic Biosciences, implanted "
                "28 to 30 days before launch. The age is not recorded in the source files."
            ),
        ),
    )

    # --- Discover tetrodes and create electrode groups ---
    tt_dirs = _tetrode_dirs(raw_root)

    device = nwbfile.create_device(
        name="tetrode_array",
        description="Multi-electrode tetrode recording array chronically implanted in hippocampal area CA1",
    )

    electrode_groups: dict[str, ElectrodeGroup] = {}
    electrode_tt_indices: dict[str, list[int]] = {}
    for tt_dir in tt_dirs:
        tt_name = tt_dir.name
        group = nwbfile.create_electrode_group(
            name=tt_name,
            description=f"Tetrode {tt_name}",
            location="CA1 field of hippocampus",
            device=device,
        )
        electrode_groups[tt_name] = group
        electrode_tt_indices[tt_name] = list(range(4 * len(electrode_tt_indices), 4 * len(electrode_tt_indices) + 4))

        # Add 4 electrodes per tetrode
        for ch in range(4):
            nwbfile.add_electrode(
                group=group,
                location="CA1 field of hippocampus",
            )

    spike_files = _read_spike_files(raw_root, stub_test=stub_test)

    epoch_windows = _add_epochs(nwbfile, spike_files)
    unit_index = _add_units(nwbfile, spike_files, electrode_tt_indices)

    # --- Add position data (flight sessions only) ---
    if any(f.cel.has_position for f in spike_files):
        task_windows = [w for name, w in epoch_windows.items() if name.rstrip("0123456789") in TASK_SESSION_TYPES]
        partner = SHARED_RECORDING.get(subject_name)
        position_files = spike_files
        if partner is not None:
            position_files = spike_files + _read_spike_files(base_dir / "1_RAW_(original_files)" / partner)
        _add_position_data(nwbfile, position_files, task_windows, meta.rat_id, partner)

    # --- Add rate maps ---
    if analyzed_root.exists():
        _add_rate_maps(nwbfile, analyzed_root, unit_index, stub_test=stub_test)

    # --- Write ---
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{subject_name}.nwb"
    with NWBHDF5IO(str(out_path), "w") as io:
        io.write(nwbfile)

    return out_path


def _add_epochs(nwbfile: NWBFile, spike_files: list[SpikeFile]) -> dict[str, tuple[float, float]]:
    """Add one epoch per original epoch name (BL1, ES, ES2, ...) and return their windows."""
    windows: dict[str, tuple[float, float]] = {}
    for f in spike_files:
        start, stop = f.cel.start_time_sec, f.cel.end_time_sec
        if f.epoch_name in windows:
            start = min(start, windows[f.epoch_name][0])
            stop = max(stop, windows[f.epoch_name][1])
        windows[f.epoch_name] = (start, stop)

    epochs = nwbfile.create_time_intervals(
        name="epochs",
        description=(
            "Task epochs within the recording session, one per epoch name in the original file names. "
            "Where the original files cut an epoch with different windows on different tetrodes, the "
            "epoch spans all of them; each unit's obs_intervals give the windows used for that unit."
        ),
    )
    epochs.add_column("epoch_name", "Epoch name from the original file names (e.g. BL2, ES, ES1)")
    epochs.add_column("session_type", "Task type: BL (Baseline), ES (Escher Staircase), or MC (Magic Carpet)")
    epochs.add_column("session_type_description", "Description of the task type")

    for name, (start, stop) in sorted(windows.items(), key=lambda item: item[1]):
        session_type = name.rstrip("0123456789")
        epochs.add_interval(
            start_time=float(start),
            stop_time=float(stop),
            epoch_name=name,
            session_type=session_type,
            session_type_description=SESSION_TYPES[session_type],
        )
    return windows


def _time_resolution(spike_times: np.ndarray) -> float:
    """Resolution of times written with six significant digits, as in the CEL files."""
    t_max = float(np.max(spike_times))
    return 10.0 ** (math.floor(math.log10(t_max)) - 5)


def _add_units(
    nwbfile: NWBFile,
    spike_files: list[SpikeFile],
    electrode_tt_indices: dict[str, list[int]],
) -> dict[tuple[str, int], int]:
    """Add one unit per (tetrode, cluster) and return the row index of each."""
    by_unit: dict[tuple[str, int], list[SpikeFile]] = defaultdict(list)
    for f in spike_files:
        by_unit[(f.tetrode, f.cel.cluster)].append(f)

    resolution = _time_resolution(np.concatenate([f.cel.spike_times[f.in_window] for f in spike_files]))

    nwbfile.units = Units(
        name="units",
        description=(
            "Sorted units. Each tetrode was cluster-cut once in Xclust over the whole recording (one "
            "TTn.pos.parms source file per tetrode), and the original spike files hold the spikes of one "
            "cluster within one epoch. A unit merges the files of one cluster on one tetrode, so it is the "
            "same sorted cluster across epochs. Spike times were written with six significant digits, so "
            "they are resolved to 1 ms below 1000 s, 10 ms from 1000 to 9999 s, and 0.1 s from 10000 s; "
            "spikes in one unit can share a time. In 48 original files the first row is a spike from before "
            "the file's epoch window; spikes outside a file's window are left out."
        ),
        electrode_table=nwbfile.electrodes,
        resolution=resolution,
    )
    nwbfile.add_unit_column("tetrode", "Tetrode name (e.g. TT0)")
    nwbfile.add_unit_column("cluster_id", "Cluster ID from Xclust sorting")
    nwbfile.add_unit_column(
        "source_files",
        "Original spike files merged into this unit, as tetrode/long name (e.g. TT1/BL2.cell5)",
        index=True,
    )

    unit_index: dict[tuple[str, int], int] = {}
    for key in sorted(by_unit, key=lambda k: (int(k[0][2:]), k[1])):
        files = sorted(by_unit[key], key=lambda f: f.cel.start_time_sec)
        times = np.sort(np.concatenate([f.cel.spike_times[f.in_window] for f in files]))
        tetrode, cluster_id = key
        unit_index[key] = len(unit_index)
        nwbfile.add_unit(
            spike_times=times,
            obs_intervals=[[f.cel.start_time_sec, f.cel.end_time_sec] for f in files],
            tetrode=tetrode,
            cluster_id=cluster_id,
            source_files=[f"{f.tetrode}/{f.long_name}" for f in files],
            electrodes=electrode_tt_indices[tetrode],
        )
    return unit_index


def _add_position_data(
    nwbfile: NWBFile,
    spike_files: list[SpikeFile],
    task_windows: list[tuple[float, float]],
    rat_id: str,
    partner: str | None,
) -> None:
    """Add the video tracker position within this rat's ES and MC epochs.

    The source holds position only as a column of the spike files, so the
    trace is sampled at spike times. Samples from all spike files of the
    recording are merged; for a shared flight recording that includes the
    other rat's files, which carry the same tracker.
    """
    times_list, xy_list = [], []
    for f in spike_files:
        xy = f.cel.position
        if xy is None:
            continue
        times_list.append(f.cel.spike_times[f.in_window])
        xy_list.append(xy[f.in_window])
    times = np.concatenate(times_list)
    xy = np.concatenate(xy_list)

    in_task = np.zeros(len(times), dtype=bool)
    for start, stop in task_windows:
        in_task |= (times >= start) & (times <= stop)
    times, xy = times[in_task], xy[in_task]

    times, first = np.unique(times, return_index=True)
    xy = xy[first]

    description = (
        "Position from the video tracker, available in the source only at spike times (the pos_x and "
        "pos_y columns of the original spike files), restricted to this rat's ES and MC epochs. "
        "Samples where the tracker lost the animal, written as (0, 255) in the source, are NaN. "
        "Occupancy cannot be computed from these samples; the rate_maps table holds the original "
        "analysis's occupancy, counted from the 60 Hz tracker."
    )
    if partner is not None:
        description += (
            f" This recording was shared with another rat (folder {partner}), and one tracker served both; "
            "the samples here merge both rats' spike files. In each window one rat ran the task while the "
            f"other rested in the sleeping pouch, so the trace in {rat_id}'s ES and MC epochs is taken to be "
            f"{rat_id}'s position, as in the original rate maps. The source files do not name the tracked animal."
        )

    behavior_module = nwbfile.create_processing_module(
        name="behavior",
        description="Behavioral data (position tracking)",
    )
    position = Position(name="position")
    position.create_spatial_series(
        name="spatial_series",
        data=xy,
        timestamps=times,
        reference_frame=(
            "Pixel coordinates (0 to 255) of the video tracker image, as written in the original spike files. "
            "The camera view and the pixel size on the apparatus are not recorded in the source files."
        ),
        unit="pixels",
        description=description,
    )
    behavior_module.add(position)


def _add_rate_maps(
    nwbfile: NWBFile,
    analyzed_root: Path,
    unit_index: dict[tuple[str, int], int],
    stub_test: bool = False,
) -> None:
    """Read .RMA files and store them as a table linked to the units table."""
    rows: list[dict] = []
    for tt_dir in _tetrode_dirs(analyzed_root):
        long_names = read_long_names(tt_dir)
        rma_paths = sorted(tt_dir.glob("*.RMA"))
        if stub_test:
            rma_paths = rma_paths[:2]
        for rma_path in rma_paths:
            long_name = long_names[rma_path.name]
            epoch_name, cluster_id = parse_long_name(long_name)
            key = (tt_dir.name, cluster_id)
            if key not in unit_index:
                print(f"  WARNING: Skipping {tt_dir.name}/{long_name}: no spike file for this cluster")
                continue
            rma = read_rma_file(rma_path)
            rows.append(dict(
                unit=unit_index[key],
                tetrode=tt_dir.name,
                cluster_id=cluster_id,
                epoch_name=epoch_name,
                source_file=long_name,
                archived_file_name=rma_path.name,
                digest=hashlib.sha256(rma_path.read_bytes()).hexdigest(),
                rate_map=rma.rate_map,
                occupancy_map=rma.occupancy_map,
            ))

    # Some files are byte-identical copies of a map stored under another cell
    # or tetrode; which cell such a map belongs to cannot be told, so all
    # copies are left out.
    digest_count = defaultdict(int)
    for r in rows:
        digest_count[r["digest"]] += 1
    for r in rows:
        if digest_count[r["digest"]] > 1:
            print(f"  WARNING: Skipping {r['tetrode']}/{r['source_file']}: identical to another map in this session")
    rows = [r for r in rows if digest_count[r["digest"]] == 1]

    if not rows:
        return

    rate_map_table = DynamicTable(
        name="rate_maps",
        description=(
            "Spatial firing rate maps and occupancy maps from the original analysis, one row per original "
            ".rmap file, linked to the unit whose cluster the file name gives. Maps are 64 x 64 bins over "
            "the tracker's 256 x 256 pixel image: bin [i, j] covers pos_y from 4i to 4i + 3 and pos_x from "
            "4j to 4j + 3. Maps that were byte-identical copies of another map in the same session were "
            "left out."
        ),
        columns=[
            DynamicTableRegion(
                name="unit",
                description="Unit in the units table (same tetrode and cluster)",
                data=[r["unit"] for r in rows],
                table=nwbfile.units,
            ),
            VectorData(name="tetrode", description="Tetrode name (e.g. TT0)", data=[r["tetrode"] for r in rows]),
            VectorData(name="cluster_id", description="Cluster ID from Xclust sorting", data=[r["cluster_id"] for r in rows]),
            VectorData(name="epoch_name", description="Epoch name from the original file name (e.g. ES, ES2, MC)", data=[r["epoch_name"] for r in rows]),
            VectorData(name="source_file", description="Original file name (e.g. ES.cell5.rmap)", data=[r["source_file"] for r in rows]),
            VectorData(
                name="archived_file_name",
                description="Windows 8.3 short name under which the file is archived in NASA OSDR OSD-968",
                data=[r["archived_file_name"] for r in rows],
            ),
            VectorData(
                name="rate_map",
                description=(
                    "64x64 spatial firing rate map from the original analysis. The original documentation "
                    "describes it as spikes in each bin divided by the bin's occupancy; the maps are smoothed "
                    "and the kernel is not recorded in the source files. Unvisited bins are NaN (-999 in the "
                    "source)."
                ),
                data=[r["rate_map"] for r in rows],
            ),
            VectorData(
                name="occupancy_map",
                description=(
                    "64x64 occupancy map: the number of 60 Hz video samples in each bin. Samples where the "
                    "tracker lost the animal (pixel 0, 255) are counted in bin [63, 0]."
                ),
                data=[r["occupancy_map"] for r in rows],
            ),
        ],
    )

    ecephys_module = nwbfile.create_processing_module(
        name="ecephys",
        description="Processed electrophysiology data including spatial firing rate maps",
    )
    ecephys_module.add(rate_map_table)
