"""Recover the original long file names behind the archive's 8.3 short names.

The source files were copied from the lab's Unix server to a Windows NT disk
with WS_FTP in December 2001, and every copy since then carries the Windows
8.3 short names (ESCELL~1.RMA, ES2BC0~1.RMA, BL950D~1.CEL). The original
names (ES.cell1.rmap, BL1.cell5) carry the epoch and the cell number, and
each tetrode folder keeps the WS_FTP.LOG that lists them in the order the
files were written. Windows assigns a short name from the long name, the
order of creation, and a checksum of the long name, so replaying that
assignment over the log recovers the name of every file.

The checksum and the assignment follow RtlGenerate8dot3Name as implemented
in ReactOS (sdk/lib/rtl/dos8dot3.c).
"""

from __future__ import annotations

import re
from pathlib import Path

_ILLEGAL = set('"*+,/:;<=>?[\\]|')
_HEX = "0123456789ABCDEF"


def _checksum(name: str) -> int:
    chars = [ord(c) for c in name]
    n = len(chars)
    if n == 0:
        return 0
    if n == 1:
        return chars[0]
    value = ((chars[0] << 8) + chars[1]) & 0xFFFF
    if n == 2:
        return value
    saved = value
    i = 0
    length = 2
    while length < n:
        i += 2
        value = ((value << 7) + chars[i]) & 0xFFFF
        value = ((saved >> 1) + (value << 8)) & 0xFFFF
        if length + 1 < n:
            value = (value + chars[i + 1]) & 0xFFFF
        saved = value
        length += 2
    return value


def _clean(text: str) -> str:
    return "".join("_" if c in _ILLEGAL else c.upper() for c in text if c > " " and c != ".")


def _candidates(long_name: str):
    """Yield the short names Windows tries for a long name, in order."""
    dot = long_name.rfind(".")
    if dot < 0:
        dot = len(long_name)
    base = _clean(long_name[:dot])[:6]
    ext = "." + _clean(long_name[dot + 1 :])[:3] if dot < len(long_name) else ""

    checksum = _checksum(long_name)
    hex_digits = "".join(_HEX[(checksum >> (4 * k)) & 0xF] for k in range(4))

    if len(base) <= 2:
        for i in range(1, 10):
            yield f"{base}{hex_digits}~{i}{ext}"
        return
    for i in range(1, 5):
        yield f"{base}~{i}{ext}"
    for i in range(1, 10):
        yield f"{base[:2]}{hex_digits}~{i}{ext}"


def assign_short_names(long_names: list[str]) -> dict[str, str]:
    """Map short name to long name for files created in the given order."""
    short_to_long: dict[str, str] = {}
    for long_name in long_names:
        for candidate in _candidates(long_name):
            if candidate not in short_to_long:
                short_to_long[candidate] = long_name
                break
    return short_to_long


def read_long_names(directory: Path) -> dict[str, str]:
    """Map each short file name in a tetrode folder to its original long name."""
    log_text = (directory / "WS_FTP.LOG").read_text(errors="ignore")
    long_names = [line.split()[-1] for line in log_text.splitlines() if "<--" in line]
    return assign_short_names(long_names)


def parse_long_name(long_name: str) -> tuple[str, int]:
    """Split a long name such as BL2.cell5 or ES.cell1.rmap into epoch name and cell number."""
    match = re.fullmatch(r"([A-Z]+\d*)\.cell(\d+)(?:\.rmap)?", long_name)
    if match is None:
        raise ValueError(f"Unexpected file name: {long_name}")
    return match.group(1), int(match.group(2))
