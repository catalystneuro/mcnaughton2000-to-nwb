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

- **Units** — spike-sorted neurons with spike times, tetrode assignment, and cluster ID
- **Electrodes** — tetrode electrode groups mapped to hippocampal CA1
- **Epochs** — task intervals (Escher Staircase, Magic Carpet, Baseline)
- **Position** — video-tracked x/y coordinates (flight sessions only)
- **Rate maps** — pre-computed 64x64 spatial firing rate and occupancy maps from the original analysis
