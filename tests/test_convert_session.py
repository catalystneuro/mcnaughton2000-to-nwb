"""Convert a small synthetic copy of a shared flight recording (Rat1 and Rat2)."""

from pathlib import Path

import numpy as np
import pytest
from pynwb import NWBHDF5IO

from mcnaughton2000_to_nwb.convert_session import convert_session
from mcnaughton2000_to_nwb.short_names import assign_short_names

FIELDS = "id t_px t_py t_pa t_pb t_mx t_my t_ma t_mb t_maxwd t_maxht time pos_x pos_y head_dir".split()


def _cel_text(cluster: int, start: str, end: str, spikes: list[tuple[float, float, float]]) -> str:
    header = [
        "%%BEGINHEADER",
        "% Program: \tXclust",
        "% Directory: \t/data/SHUTTLE/e100-04.20.98-09:57/TT1",
        "% Source file:\tTT1.pos.parms",
        "% Fields: \t" + "\t".join(FIELDS),
        f"% Cluster:\t{cluster}",
        f"% Start time:\t{start}",
        f"% End time:\t{end}",
        "%%ENDHEADER",
    ]
    rows = [
        "\t".join(str(v) for v in [i, 300, 700, 300, 400, 1, 1, 1, 1, 7, 1000, t, x, y, -50])
        for i, (t, x, y) in enumerate(spikes)
    ]
    return "\n".join(header + rows) + "\n"


def _rma_bytes(rate: np.ndarray, occupancy: np.ndarray) -> bytes:
    return rate.astype(">f4").tobytes() + occupancy.astype(">i4").tobytes()


def _write_folder(folder: Path, files: dict[str, bytes | str]) -> None:
    """Write files under their 8.3 short names, with a WS_FTP.LOG of their long names."""
    folder.mkdir(parents=True)
    long_to_short = {v: k for k, v in assign_short_names(list(files)).items()}
    log_lines = []
    for long_name, content in files.items():
        path = folder / long_to_short[long_name]
        if isinstance(content, str):
            path.write_text(content)
        else:
            path.write_bytes(content)
        log_lines.append(f"101.12.08 16:20 A D:\\jim's\\NASA\\{long_name} <-- nba18 /data/jim/NASA {long_name}")
    (folder / "WS_FTP.LOG").write_text("\n".join(log_lines) + "\n")


@pytest.fixture(scope="module")
def nwb_dir(tmp_path_factory) -> Path:
    base = tmp_path_factory.mktemp("McNaughton_Neurolab")
    raw = base / "1_RAW_(original_files)"
    # Rat1 rests (BL1) from 10 to 20 s and runs the staircase (ES) from 30 to 40 s;
    # Rat2 does the reverse. Both spike files carry the one tracker's position.
    _write_folder(raw / "FD4RAT1" / "TT1", {
        # The first row lies before the window, as in 48 of the original files.
        "BL1.cell1": _cel_text(1, "0:10", "0:20", [(9.5, 5, 5), (12, 1, 1), (15, 2, 2)]),
        "ES.cell1": _cel_text(1, "0:30", "0:40", [(31, 10, 20), (33, 0, 255), (35, 30, 40)]),
        "ES.cell2": _cel_text(2, "0:30", "0:40", [(32, 50, 60)]),
        "ES.cell3": _cel_text(3, "0:30", "0:40", [(36, 52, 62)]),
    })
    _write_folder(raw / "FD4RAT2" / "TT12", {
        "ES.cell1": _cel_text(1, "0:10", "0:20", [(11, 100, 100), (14, 110, 110)]),
        "BL1.cell1": _cel_text(1, "0:30", "0:40", [(34, 70, 80)]),
    })

    rate = np.full((64, 64), -999.0)
    rate[5, 2] = 4.0
    occupancy = np.zeros((64, 64))
    occupancy[5, 2] = 60
    duplicate = _rma_bytes(np.full((64, 64), 1.0), occupancy)
    _write_folder(base / "ANALYZED_(original_files)" / "FD4RAT1" / "TT1", {
        "ES.cell1.rmap": _rma_bytes(rate, occupancy),
        "ES.cell2.rmap": duplicate,
        "ES.cell3.rmap": duplicate,
    })

    out = base / "nwb"
    convert_session(base, "FD4RAT1", out)
    convert_session(base, "FD4RAT2", out)
    return out


def test_units(nwb_dir):
    with NWBHDF5IO(str(nwb_dir / "FD4RAT1.nwb"), "r") as io:
        units = io.read().units.to_dataframe()
    unit = units.iloc[0]
    assert (unit.tetrode, unit.cluster_id) == ("TT1", 1)
    np.testing.assert_array_equal(unit.spike_times, [12, 15, 31, 33, 35])
    assert list(unit.source_files) == ["TT1/BL1.cell1", "TT1/ES.cell1"]
    np.testing.assert_array_equal(unit.obs_intervals, [[10, 20], [30, 40]])


def test_epochs(nwb_dir):
    with NWBHDF5IO(str(nwb_dir / "FD4RAT1.nwb"), "r") as io:
        epochs = io.read().epochs.to_dataframe()
    assert list(epochs.epoch_name) == ["BL1", "ES"]
    assert list(epochs.session_type) == ["BL", "ES"]


def test_position_is_restricted_to_own_task_epochs(nwb_dir):
    with NWBHDF5IO(str(nwb_dir / "FD4RAT1.nwb"), "r") as io:
        series = io.read().processing["behavior"]["position"]["spatial_series"]
        timestamps, data = series.timestamps[:], series.data[:]
    # Samples from both rats' files within Rat1's ES epoch only.
    np.testing.assert_array_equal(timestamps, [31, 32, 33, 34, 35, 36])
    assert np.isnan(data[2]).all()  # lost position (0, 255)
    np.testing.assert_array_equal(data[3], [70, 80])  # from Rat2's spike file

    with NWBHDF5IO(str(nwb_dir / "FD4RAT2.nwb"), "r") as io:
        timestamps = io.read().processing["behavior"]["position"]["spatial_series"].timestamps[:]
    np.testing.assert_array_equal(timestamps, [11, 12, 14, 15])


def test_rate_maps_link_to_units(nwb_dir):
    with NWBHDF5IO(str(nwb_dir / "FD4RAT1.nwb"), "r") as io:
        nwbfile = io.read()
        table = nwbfile.processing["ecephys"]["rate_maps"]
        # The two identical maps are left out.
        assert list(table["source_file"][:]) == ["ES.cell1.rmap"]
        unit = table["unit"][0]
        assert (unit.tetrode.iloc[0], unit.cluster_id.iloc[0]) == ("TT1", 1)
        rate = table["rate_map"][0]
    assert rate[5, 2] == 4.0
    assert np.isnan(rate[0, 0])
