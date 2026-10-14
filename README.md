# mcnaughton2000-to-nwb

Convert the McNaughton Neurolab (STS-90) hippocampal place cell dataset to [NWB](https://www.nwb.org/).

## Dataset

This converts data from the Neurolab Space Shuttle mission (April 1998), in which hippocampal CA1 place cells were recorded from three freely moving rats in microgravity. The rats navigated a 3D "Escher Staircase" track and a flat "Magic Carpet" track while spike data was collected via chronically implanted multi-tetrode arrays.

**Publication:** Knierim, J. J., McNaughton, B. L. & Poe, G. R. Three-dimensional spatial selectivity of hippocampal neurons during space flight. *Nature Neuroscience* **3**, 209-210 (2000). [doi:10.1038/72910](https://doi.org/10.1038/72910)

**Source data:** [NASA OSDR OSD-968](https://osdr.nasa.gov/bio/repo/data/studies/OSD-968)

**NWB dataset:** [DANDI:001754](https://dandiarchive.org/dandiset/001754)

## Installation

```bash
pip install -e .
```

## Usage

Convert all 8 sessions (3 rats, preflight + flight days):

```bash
python -m mcnaughton2000_to_nwb.convert_all_sessions \
    --base /path/to/McNaughton_Neurolab \
    --output ./output
```

Convert a single subject:

```bash
python -m mcnaughton2000_to_nwb.convert_all_sessions \
    --base /path/to/McNaughton_Neurolab \
    --output ./output \
    --subjects FD9RAT1
```

The `--base` directory should contain `1_RAW_(original_files)/` and `ANALYZED_(original_files)/` subdirectories.

## Output

Each NWB file contains:

- **Units**: sorted units with spike times, tetrode, cluster ID, the original spike files merged into each unit, and the window of each file as `obs_intervals`
- **Electrodes**: tetrode electrode groups mapped to hippocampal CA1
- **Epochs**: one interval per original epoch name (`BL1`, `ES`, `ES2`, `MC`, ...), with its task type: Escher Staircase, Magic Carpet, or Baseline (rest in the Rodent Sleeping Pouch)
- **Position**: video tracker pixel coordinates at spike times, within the rat's ES and MC epochs (flight sessions only)
- **Rate maps**: the original 64x64 firing rate and occupancy maps, each linked to its unit

## Notes on the Source Files

The archived files carry Windows 8.3 short names (`ESCELL~1.RMA`, `BL950D~1.CEL`). The original names (`ES.cell1.rmap`, `BL2.cell5`) give the epoch and the cluster. They are recovered from the `WS_FTP.LOG` in each tetrode folder by replaying the Windows short-name assignment (`short_names.py`), which reproduces the archived names in every folder.

On each flight day Rat 1 and Rat 2 were recorded in one acquisition session with one video tracker, and took turns: one rat ran the task while the other rested in the sleeping pouch. Position in each rat's file is therefore limited to that rat's ES and MC epochs, and merges the tracker samples from both rats' spike files. The source holds position only at spike times, and none for the preflight sessions.

Spike times in the source are written with six significant digits (10 ms resolution from 1000 s, 0.1 s from 10000 s). Rate maps that are byte-identical copies of another map in the same session are left out.

## Tests

```bash
pip install -e ".[test]"
pytest tests
```
