"""Convert one subject-session of the McNaughton Neurolab dataset to NWB.

Each subject folder (e.g. FD4RAT1) represents one continuous recording session
with interleaved task epochs (BL=Baseline, ES=Escher Staircase, MC=Magic Carpet).
The CEL files contain sorted spike data from Xclust with timestamps in seconds
relative to recording start. Flight sessions also include position (pos_x, pos_y)
at each spike time.
"""

from __future__ import annotations

import uuid
from collections import defaultdict
from pathlib import Path

import numpy as np
from pynwb import NWBFile, NWBHDF5IO, ProcessingModule
from pynwb.behavior import Position, SpatialSeries
from pynwb.ecephys import ElectrodeGroup
from pynwb.epoch import TimeIntervals
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
    SPECIES,
    STRAIN,
    get_session_metadata,
)
from .rma_file import read_rma_file


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
            age="/",
            description=f"Adult male Fischer 344 rat ({meta.rat_id})",
        ),
    )

    # --- Discover tetrodes and create electrode groups ---
    tt_dirs = sorted(
        [p for p in raw_root.iterdir() if p.is_dir() and p.name.upper().startswith("TT")],
        key=lambda p: int(p.name[2:]),
    )

    device = nwbfile.create_device(
        name="tetrode_array",
        description="Multi-electrode tetrode recording array chronically implanted in hippocampal area CA1",
    )

    electrode_groups: dict[str, ElectrodeGroup] = {}
    for tt_dir in tt_dirs:
        tt_name = tt_dir.name
        group = nwbfile.create_electrode_group(
            name=tt_name,
            description=f"Tetrode {tt_name}",
            location="CA1 field of hippocampus",
            device=device,
        )
        electrode_groups[tt_name] = group

        # Add 4 electrodes per tetrode
        for ch in range(4):
            nwbfile.add_electrode(
                group=group,
                location="CA1 field of hippocampus",
            )

    # --- Parse all CEL files ---
    all_cel_files: list[tuple[str, CelFile]] = []  # (tetrode_name, cel_file)
    for tt_dir in tt_dirs:
        tt_name = tt_dir.name
        cel_paths = sorted(
            p for p in tt_dir.iterdir()
            if p.is_file() and p.suffix.upper() in (".CEL", ".CELL")
        )
        if stub_test:
            cel_paths = cel_paths[:2]
        for cel_path in cel_paths:
            try:
                cel = read_cel_file(cel_path)
                all_cel_files.append((tt_name, cel))
            except Exception as e:
                print(f"  WARNING: Failed to parse {cel_path}: {e}")

    # --- Build epochs from unique (start_time, end_time, session_type) ---
    epoch_set: dict[tuple[float, float, str], None] = {}
    for tt_name, cel in all_cel_files:
        key = (cel.start_time_sec, cel.end_time_sec, cel.session_type)
        if not np.isnan(key[0]) and not np.isnan(key[1]):
            epoch_set[key] = None

    epochs = nwbfile.create_time_intervals(
        name="epochs",
        description="Task epochs within the recording session",
    )
    epochs.add_column("session_type", "Task type: BL (Baseline), ES (Escher Staircase), or MC (Magic Carpet)")
    epochs.add_column("session_type_description", "Description of the task type")

    for start, stop, stype in sorted(epoch_set.keys()):
        epochs.add_interval(
            start_time=float(start),
            stop_time=float(stop),
            session_type=stype,
            session_type_description=SESSION_TYPES.get(stype, "Unknown task type"),
        )

    # --- Add units (merge CEL files by tetrode + cluster) ---
    nwbfile.add_unit_column("tetrode", "Tetrode name (e.g. TT0)")
    nwbfile.add_unit_column("cluster_id", "Cluster ID from Xclust sorting")

    # Build electrode index: tetrode name → electrode row indices
    electrode_tt_indices: dict[str, list[int]] = {}
    idx = 0
    for tt_dir in tt_dirs:
        tt_name = tt_dir.name
        electrode_tt_indices[tt_name] = list(range(idx, idx + 4))
        idx += 4

    # Group CEL files by (tetrode, cluster_id) to merge the same neuron across epochs
    unit_spikes: dict[tuple[str, int], list[np.ndarray]] = defaultdict(list)
    for tt_name, cel in all_cel_files:
        cluster_id = cel.cluster if cel.cluster is not None else -1
        spike_times = cel.spike_times
        mask = ~np.isnan(spike_times)
        spike_times = spike_times[mask]
        if len(spike_times) > 0:
            unit_spikes[(tt_name, cluster_id)].append(spike_times)

    for tt_name, cluster_id in sorted(unit_spikes.keys()):
        merged = np.sort(np.concatenate(unit_spikes[(tt_name, cluster_id)]))
        nwbfile.add_unit(
            spike_times=merged,
            tetrode=tt_name,
            cluster_id=cluster_id,
            electrodes=electrode_tt_indices[tt_name],
        )

    # --- Add position data (flight sessions only) ---
    has_position = any(cel.has_position for _, cel in all_cel_files)
    if has_position:
        _add_position_data(nwbfile, all_cel_files)

    # --- Add rate maps ---
    if analyzed_root.exists():
        _add_rate_maps(nwbfile, analyzed_root, stub_test=stub_test)

    # --- Write ---
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{subject_name}.nwb"
    with NWBHDF5IO(str(out_path), "w") as io:
        io.write(nwbfile)

    return out_path


def _add_position_data(
    nwbfile: NWBFile,
    all_cel_files: list[tuple[str, CelFile]],
) -> None:
    """Extract position data from CEL files and add as SpatialSeries.

    Position (pos_x, pos_y) is recorded at each spike time. We merge across
    all cells, sort by time, and deduplicate to reconstruct the position trace
    as sampled by the video tracker.
    """
    time_pos_pairs: list[tuple[float, float, float]] = []

    for _, cel in all_cel_files:
        if not cel.has_position:
            continue
        times = cel.spike_times
        x = cel.pos_x
        y = cel.pos_y
        if x is None or y is None:
            continue

        mask = ~(np.isnan(times) | np.isnan(x) | np.isnan(y))
        for t, px, py in zip(times[mask], x[mask], y[mask]):
            time_pos_pairs.append((t, px, py))

    if not time_pos_pairs:
        return

    # Sort by time and deduplicate
    time_pos_pairs.sort(key=lambda tup: tup[0])
    times_arr = np.array([t for t, _, _ in time_pos_pairs])
    x_arr = np.array([px for _, px, _ in time_pos_pairs])
    y_arr = np.array([py for _, _, py in time_pos_pairs])

    # Deduplicate: keep first occurrence of each unique time
    _, unique_idx = np.unique(times_arr, return_index=True)
    times_arr = times_arr[unique_idx]
    x_arr = x_arr[unique_idx]
    y_arr = y_arr[unique_idx]

    behavior_module = nwbfile.create_processing_module(
        name="behavior",
        description="Behavioral data (position tracking)",
    )

    position = Position(name="position")
    position.create_spatial_series(
        name="spatial_series",
        data=np.column_stack([x_arr, y_arr]),
        timestamps=times_arr,
        reference_frame="Video tracker pixel coordinates",
        unit="pixels",
        description="Animal position from video tracking, sampled at spike occurrence times",
    )
    behavior_module.add(position)


def _add_rate_maps(
    nwbfile: NWBFile,
    analyzed_root: Path,
    stub_test: bool = False,
) -> None:
    """Read .RMA files and store as a DynamicTable with 2D array columns."""
    from hdmf.common import DynamicTable, VectorData

    tt_dirs = sorted(
        [p for p in analyzed_root.iterdir() if p.is_dir() and p.name.upper().startswith("TT")],
        key=lambda p: int(p.name[2:]),
    )

    # Collect all RMA data
    rows: list[dict] = []
    for tt_dir in tt_dirs:
        rma_paths = sorted(tt_dir.glob("*.RMA"))
        if stub_test:
            rma_paths = rma_paths[:2]
        for rma_path in rma_paths:
            try:
                rma = read_rma_file(rma_path)
                rows.append(dict(
                    tetrode=tt_dir.name,
                    source_file=rma_path.name,
                    session_type=rma.session_type,
                    cell_number=rma.cell_number if rma.cell_number is not None else -1,
                    rate_map=rma.rate_map,
                    occupancy_map=rma.occupancy_map,
                ))
            except Exception as e:
                print(f"  WARNING: Failed to parse {rma_path}: {e}")

    if not rows:
        return

    # Build the DynamicTable
    n = len(rows)
    rate_map_table = DynamicTable(
        name="rate_maps",
        description=(
            "Spatial firing rate maps and occupancy maps from legacy analysis. "
            "Each row corresponds to one .RMA file. Cell-specific maps (cell_number != -1) "
            "can be matched to the cluster_id in the units table "
            "for the same tetrode and session_type."
        ),
        columns=[
            VectorData(name="tetrode", description="Tetrode name (e.g. TT0)", data=[r["tetrode"] for r in rows]),
            VectorData(name="source_file", description="Original .RMA filename", data=[r["source_file"] for r in rows]),
            VectorData(name="session_type", description="Task type (ES/MC/BL)", data=[r["session_type"] for r in rows]),
            VectorData(name="cell_number", description="Cell number from CELL~N filename; -1 if not a per-cell map", data=[r["cell_number"] for r in rows]),
            VectorData(name="rate_map", description="64x64 spatial firing rate map (Hz)", data=[r["rate_map"] for r in rows]),
            VectorData(name="occupancy_map", description="64x64 spatial occupancy map (bin visit counts)", data=[r["occupancy_map"] for r in rows]),
        ],
    )

    ecephys_module = nwbfile.create_processing_module(
        name="ecephys",
        description="Processed electrophysiology data including spatial firing rate maps",
    )
    ecephys_module.add(rate_map_table)
