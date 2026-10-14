# McNaughton Neurolab NWB Conversion Notes

## Dataset
- **Publication**: Knierim, McNaughton & Poe (2000) "Three-dimensional spatial selectivity of hippocampal neurons during space flight." *Nature Neuroscience* 3(3):209-212. DOI: 10.1038/72910
- **Mission**: Neurolab STS-90 Space Shuttle mission, April–May 1998
- **Source data**: NASA OSDR legacy dataset, originally in Xclust .CEL spike files and binary .RMA rate maps

## Subjects
- 3 rats (Fischer 344, adult males) with tetrode arrays in hippocampal area CA1
- Rat 1 ↔ PREFLI~1, FD4RAT1, FD9RAT1
- Rat 2 ↔ PREFLI~2, FD4RAT2, FD9RAT2
- Rat 3 ↔ PREFLI~3, FD4RAT3 (no FD9 — technical issues)

## Recording Sessions
| Subject | Date | Time | Context |
|---------|------|------|---------|
| PREFLI~2 | 1998-04-13 | 16:37 UTC | Preflight, Rat 2 |
| PREFLI~1 | 1998-04-14 | 12:53 UTC | Preflight, Rat 1 |
| PREFLI~3 | 1998-04-14 | 13:49 UTC | Preflight, Rat 3 |
| FD4RAT1+2 | 1998-04-20 | 09:57 UTC | Flight Day 4, Rats 1+2 in one session, taking turns |
| FD4RAT3 | 1998-04-20 | 15:28 UTC | Flight Day 4, Rat 3 separate |
| FD9RAT1+2 | 1998-04-25 | 12:45 UTC | Flight Day 9, Rats 1+2 in one session, taking turns |

## Data Structure

### CEL Files (Raw Sorted Spikes)
- ASCII format from Xclust spike sorting software
- Fields: id, t_px, t_py, t_pa, t_pb (spike waveform peak amplitudes on 4 channels), t_mx, t_my, t_ma, t_mb (integrated waveform areas), t_maxwd, t_maxht (waveform shape), time (seconds), pos_x, pos_y (pixels), head_dir
- **Time units are SECONDS** (not 10µs ticks as the original ChatGPT script assumed)
- Position data (pos_x, pos_y) available in flight sessions only; missing in preflight
- head_dir is always -50 (sentinel for "no data") — not stored in NWB
- Each file represents one sorted cluster from one task epoch on one tetrode
- Each tetrode was cluster-cut once (one TTn.pos.parms source file), so a cluster ID names the same cluster in every epoch
- Times are written with six significant digits: 10 ms resolution from 1000 s, 0.1 s from 10000 s
- In 48 files the first row is a spike from before the header's Start time; it also appears in the previous epoch's file. Spikes outside the header window are dropped.
- The tracker writes (0, 255) when it loses the animal; stored as NaN

### File Names
- All archived copies (local and OSD-968) use Windows 8.3 short names. The original long names (BL2.cell5, ES.cell1.rmap) are listed in each folder's WS_FTP.LOG in the order the files were written, and replaying the Windows short-name assignment over that list reproduces the archived names in all 133 folders (one stray "BL1~1 copy.csv" aside). `short_names.py` does this.
- The cell number in each long spike-file name equals the Cluster header in all 982 spike files.

### RMA Files (Rate Maps)
- Binary format: 32,768 bytes
- First 4096 float32 big-endian → 64×64 firing rate map (Hz)
- Next 4096 int32 big-endian → 64×64 occupancy map (counts)
- One map per cluster and epoch (ES.cellN.rmap, MC.cellN.rmap); the ~N in short names is not a cell number
- Unvisited bins are -999 (stored as NaN). Maps are smoothed with an unrecorded kernel.
- 64x64 bins over the 256x256 pixel tracker image, row = pos_y // 4, column = pos_x // 4; occupancy counts 60 Hz video samples
- 14 maps in six groups are byte-identical copies of each other (e.g. FD4RAT1 TT1/ES.cell1 and TT11/ES.cell1) and are left out

### Session Types
- **BL**: Baseline — rest in the Rodent Sleeping Pouch (original README)
- **ES**: Escher Staircase — 3D track with 90° yaw/pitch turns
- **MC**: Magic Carpet — flat 2D track

## NWB File Structure (per subject-session)
- **Units table**: spike_times (seconds), obs_intervals, tetrode, cluster_id, source_files, electrodes
- **Epochs**: TimeIntervals with start/end times, epoch_name and session_type for BL/ES/MC periods
- **Position**: SpatialSeries (pos_x, pos_y in pixels) within the rat's ES and MC epochs — flight sessions only
- **Rate maps**: DynamicTable in processing/ecephys with a unit link and 64×64 rate_map and occupancy_map columns

## Issues Found in Original (ChatGPT-generated) NWB Files
1. **Time units wrong**: multiplied by 1e-5 when times were already in seconds
2. **No position data**: pos_x/pos_y columns were discarded
3. **No epochs**: task periods not marked
4. **Placeholder dates**: used 2001-01-01 instead of actual 1998 dates
5. **Rate maps as TimeSeries**: misused TimeSeries with fake timestamps
6. **Missing metadata**: no species, strain, experimenter, DOI, experiment description
7. **No electrode groups**: tetrodes not properly modeled

## Open Questions
- Exact age of rats at time of recording (not in the source files or OSD-968)
- Camera view and pixel size on the apparatus. OSD-968 says two synchronized cameras recorded position; the spike files hold one x, y pair.
- Whose position the shared flight tracker followed in each window. By the protocol it is the rat in an ES or MC epoch, but a reproduction attempt (catalystneuro/reproduction-registry, 001754-knierim-2000) found Rat2's units not spatially selective with it in any window.
- Whether Rat 1↔PREFLI~1 mapping is correct (tetrode numbering differs between preflight and flight for Rat 2)
