"""Read run identity out of the FASTQ records themselves.

The parser below the marker is a verbatim copy of `wesingest/header.py` in the
lab's WES repository, which is its canonical home. It is duplicated rather than
imported so this workflow stays clone-and-run for deployments that have no
wesingest; `tests/test_vendored_parser.py` compares the two sources object by
object, so a divergence fails a test instead of going unnoticed.

Do not edit anything above the "pipeline-local" marker. Changes belong upstream
and then get re-copied here.

`read_first_header` is part of the vendored surface and is kept identical for
that comparison. The pipeline itself reads through `read_first_line`, because it
also needs the index, which lives in the header's comment field and is not part
of what upstream parses (see docs/TODO.md item 2 in the WES repo).
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

# --------------------------------------------------------------------------
# Vendored verbatim from /mnt/data/NGS/WES/scripts/wesingest/header.py
# --------------------------------------------------------------------------

# CASAVA 1.8+ header. Everything this repository has ever held is this shape; the
# pre-1.8 form (@inst:lane:tile:x:y#index/read) carries no flowcell and would defeat
# the point, so it is not parsed -- an unrecognised header leaves the fields blank
# rather than guessing.
HEADER_RE = re.compile(
    r"^@(?P<instrument>[^:\s]+):(?P<run_number>\d+):(?P<flowcell>[^:\s]+):"
    r"(?P<lane>\d+):\d+:[-\d]+:[-\d]+"
)


class ReadHeader:
    """The run-identifying fields of one read, plus the read's length."""

    __slots__ = ("instrument", "run_number", "flowcell", "lane", "read_length")

    def __init__(self, instrument: str, run_number: str, flowcell: str, lane: str,
                 read_length: int = 0):
        self.instrument = instrument
        self.run_number = run_number
        self.flowcell = flowcell
        self.lane = lane                    # "L001" form, to compare with filenames
        self.read_length = read_length

    def __repr__(self) -> str:
        return (f"ReadHeader({self.instrument}:{self.run_number}:"
                f"{self.flowcell}:{self.lane})")


def parse_header(line: str, read_length: int = 0) -> ReadHeader | None:
    m = HEADER_RE.match(line.strip())
    if not m:
        return None
    return ReadHeader(
        instrument=m.group("instrument"),
        run_number=m.group("run_number"),
        flowcell=m.group("flowcell"),
        lane=f"L{int(m.group('lane')):03d}",
        read_length=read_length,
    )


def read_first_header(path: Path) -> ReadHeader | None:
    """Parse the first record of a gzipped FASTQ. None if unreadable or unrecognised."""
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            line = fh.readline()
            if not line.startswith("@"):
                return None
            seq = fh.readline().strip()
        return parse_header(line, read_length=len(seq))
    except (OSError, EOFError):
        return None


# --------------------------------------------------------------------------
# Pipeline-local additions. Not part of the vendored surface.
# --------------------------------------------------------------------------

# The comment field of a CASAVA 1.8+ header: "<read>:<filtered>:<control>:<index>".
INDEX_RE = re.compile(r"^@\S+\s+[12]:[YN]:\d+:(?P<index>\S+)")

# A usable barcode is unambiguous bases, single or dual-indexed. Real first
# records routinely carry a miscalled base -- two of three files sampled from
# this cohort have an N in the index -- and one sequencing error must not be
# baked into every read group in the run.
VALID_INDEX_RE = re.compile(r"^[ACGT]+(\+[ACGT]+)?$")


def parse_index(line: str) -> str | None:
    """The barcode from a read header, or None if absent or ambiguous."""
    m = INDEX_RE.match(line.strip())
    if not m:
        return None
    index = m.group("index")
    return index if VALID_INDEX_RE.match(index) else None


def read_first_line(path: Path | str) -> str | None:
    """The first line of a gzipped FASTQ. None if unreadable or not a record.

    One gzip block per file, which is what makes reading every unit's header
    cheap enough to do on every DAG build.
    """
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            line = fh.readline()
    except (OSError, EOFError):
        return None
    return line if line.startswith("@") else None
