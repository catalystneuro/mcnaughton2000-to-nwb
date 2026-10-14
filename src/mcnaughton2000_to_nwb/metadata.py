"""Session and subject metadata for the McNaughton Neurolab dataset.

The Neurolab STS-90 mission flew April 17 – May 3, 1998. Three rats with
hippocampal tetrode implants were recorded preflight (ground) and in-flight
on Flight Day 4 (April 20) and Flight Day 9 (April 25).

Rat 1 and Rat 2 were recorded in the same acquisition session on each flight
day (different tetrode banks) and took turns on the task apparatus. Rat 3 was
recorded in a separate session.

Reference: Knierim, McNaughton & Poe (2000) Nature Neuroscience 3(3):209-212
DOI: 10.1038/72910
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


RELATED_PUBLICATIONS = ["doi:10.1038/72910"]

EXPERIMENT_DESCRIPTION = (
    "Hippocampal place cell recordings from rats aboard the Space Shuttle "
    "during the Neurolab STS-90 mission (April–May 1998). Three rats with "
    "chronically implanted tetrode arrays in hippocampal area CA1 were trained "
    "to traverse a three-dimensional track ('Escher staircase') and a flat "
    "two-dimensional track ('Magic Carpet') for medial forebrain bundle "
    "stimulation reward. Recordings were made preflight on the ground and "
    "in-flight on Flight Day 4 and Flight Day 9. In baseline (BL) epochs the "
    "rat rested quietly in the Rodent Sleeping Pouch. Rat 1 and Rat 2 were "
    "recorded in the same acquisition session on each flight day and took "
    "turns: while one rat ran the Escher staircase or Magic Carpet, the other "
    "was in a baseline epoch. The source files hold sorted spikes and the "
    "original rate maps; no EEG or other continuous voltage was archived."
)

INSTITUTION = "University of Arizona"
LAB = "McNaughton Lab"

SPECIES = "Rattus norvegicus"
STRAIN = "Fischer 344"
SEX = "M"

SESSION_TYPES = {
    "BL": "Baseline — the rat rested quietly in the Rodent Sleeping Pouch",
    "ES": "Escher Staircase — three-dimensional track with 90° yaw and pitch turns",
    "MC": "Magic Carpet — flat two-dimensional track",
}


@dataclass
class SessionMetadata:
    """Metadata for a single subject-session (one NWB file)."""

    subject_folder: str  # e.g. "FD4RAT1"
    rat_id: str  # e.g. "Rat1"
    session_date: datetime
    session_description: str
    recording_directory: str  # original path from CEL header


# Map each subject folder to its rat identity and session context.
# Recording directory dates (from CEL headers) provide the actual session datetimes.
SUBJECT_SESSION_MAP: dict[str, dict] = {
    "FD4RAT1": dict(
        rat_id="Rat1",
        session_date=datetime(1998, 4, 20, 9, 57, tzinfo=timezone.utc),
        session_description="Flight Day 4 recording — Rat 1. Recorded in the same acquisition session as Rat 2; the two rats took turns on the task apparatus.",
    ),
    "FD4RAT2": dict(
        rat_id="Rat2",
        session_date=datetime(1998, 4, 20, 9, 57, tzinfo=timezone.utc),
        session_description="Flight Day 4 recording — Rat 2. Recorded in the same acquisition session as Rat 1; the two rats took turns on the task apparatus.",
    ),
    "FD4RAT3": dict(
        rat_id="Rat3",
        session_date=datetime(1998, 4, 20, 15, 28, tzinfo=timezone.utc),
        session_description="Flight Day 4 recording — Rat 3. Separate recording session; partial data recovered due to technical issues.",
    ),
    "FD9RAT1": dict(
        rat_id="Rat1",
        session_date=datetime(1998, 4, 25, 12, 45, tzinfo=timezone.utc),
        session_description="Flight Day 9 recording — Rat 1. Recorded in the same acquisition session as Rat 2; the two rats took turns on the task apparatus.",
    ),
    "FD9RAT2": dict(
        rat_id="Rat2",
        session_date=datetime(1998, 4, 25, 12, 45, tzinfo=timezone.utc),
        session_description="Flight Day 9 recording — Rat 2. Recorded in the same acquisition session as Rat 1; the two rats took turns on the task apparatus.",
    ),
    "PREFLI~1": dict(
        rat_id="Rat1",
        session_date=datetime(1998, 4, 14, 12, 53, tzinfo=timezone.utc),
        session_description="Preflight ground recording — Rat 1. Recorded 3 days before launch at Kennedy Space Center.",
    ),
    "PREFLI~2": dict(
        rat_id="Rat2",
        session_date=datetime(1998, 4, 13, 16, 37, tzinfo=timezone.utc),
        session_description="Preflight ground recording — Rat 2. Recorded 4 days before launch at Kennedy Space Center.",
    ),
    "PREFLI~3": dict(
        rat_id="Rat3",
        session_date=datetime(1998, 4, 14, 13, 49, tzinfo=timezone.utc),
        session_description="Preflight ground recording — Rat 3. Recorded 3 days before launch at Kennedy Space Center.",
    ),
}


# Subject folders recorded in the same acquisition session, which share one
# video tracker (same recording directory in the spike file headers).
SHARED_RECORDING: dict[str, str] = {
    "FD4RAT1": "FD4RAT2",
    "FD4RAT2": "FD4RAT1",
    "FD9RAT1": "FD9RAT2",
    "FD9RAT2": "FD9RAT1",
}


def get_session_metadata(subject_folder: str) -> SessionMetadata:
    info = SUBJECT_SESSION_MAP[subject_folder]
    return SessionMetadata(
        subject_folder=subject_folder,
        recording_directory="",
        **info,
    )


def parse_session_datetime_from_directory(directory: str) -> datetime | None:
    """Extract recording datetime from the original directory path.

    Example: '/data/SHUTTLE/e100-04.20.98-09:57/TT0' → 1998-04-20 09:57 UTC
    """
    m = re.search(r"e100-(\d{2})\.(\d{2})\.(\d{2})-(\d{2}):(\d{2})", directory)
    if not m:
        return None
    month, day, year_2d, hour, minute = (int(x) for x in m.groups())
    year = 1900 + year_2d if year_2d > 50 else 2000 + year_2d
    return datetime(year, month, day, hour, minute, tzinfo=timezone.utc)
