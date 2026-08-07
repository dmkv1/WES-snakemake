"""Read run identity out of the FASTQ records themselves.

A samplesheet may name nothing but `sample,fq1,fq2`, so flowcell, lane and
barcode have to be recoverable from the reads alone. That is what this module
does, and it is self-contained: no samplesheet, no sidecar, no ingest tool.

The parser fails open. A header it cannot read, a truncated payload or a file
that is not gzip at all degrades to a blank field rather than an exception --
unresolvable provenance becomes a positional unit token, is reported through
`rg_source`, and the pipeline still returns calls.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from pathlib import Path

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


# The comment field of a CASAVA 1.8+ header: "<read>:<filtered>:<control>:<index>".
INDEX_RE = re.compile(r"^@\S+\s+[12]:[YN]:\d+:(?P<index>\S+)")

# A usable barcode is unambiguous bases, single or dual-indexed. Index reads are
# the lowest-quality cycles on the flowcell and a miscalled base is common -- two
# of three files sampled from this cohort carry an N in the very first record --
# so a record with an ambiguous index is skipped rather than trusted.
VALID_INDEX_RE = re.compile(r"^[ACGT]+(\+[ACGT]+)?$")

# Unambiguous indexes wanted before the barcode vote is called, and the hard cap
# on records read to find them. One miscalled index base must not decide the read
# group for a whole unit, so the barcode is a vote rather than whatever the first
# record happened to say.
#
# The count that matters is votes, not records. A fixed window would be the wrong
# shape: the head of a FASTQ is the worst part of the run, because those records
# come from the first tile at the edge of the flowcell where index reads fail
# most. In this cohort one file's first unambiguous index is record 110 and
# another's is record 98, while both are over 90 percent clean further in. Any
# window short enough to be cheap is short enough to land entirely inside that
# bad head; reading until the evidence arrives costs nothing on a clean file and
# reads on only for the files that need it.
INDEX_MIN_VOTES = 200
INDEX_MAX_RECORDS = 5000

# A demultiplexed FASTQ holds one barcode, so the winning index takes essentially
# every unambiguous record. A lower share means the file mixes barcodes, and the
# consensus is then a majority rather than a fact about the file.
MIXED_INDEX_THRESHOLD = 0.9


def parse_index(line: str) -> str | None:
    """The barcode from a read header, or None if absent or ambiguous."""
    m = INDEX_RE.match(line.strip())
    if not m:
        return None
    index = m.group("index")
    return index if VALID_INDEX_RE.match(index) else None


def read_headers(path: Path | str, records: int,
                 min_votes: int | None = None) -> list[str]:
    """Header lines from the head of a gzipped FASTQ. Empty list if unreadable.

    Reads at most `records` of them. When `min_votes` is given, reading also
    stops early once that many records have carried an unambiguous index, which
    is what keeps the cost proportional to how dirty the file actually is.

    Stops at the first line that does not look like a header, which keeps a file
    whose four-line phase has drifted from being read as if it had not. A gzip
    error mid-stream keeps the records already collected, so a truncated file
    degrades to less evidence rather than to none.
    """
    lines: list[str] = []
    votes = 0
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for n, line in enumerate(fh):
                if n % 4:
                    continue
                if not line.startswith("@"):
                    break
                lines.append(line)
                if min_votes is not None and parse_index(line):
                    votes += 1
                    if votes >= min_votes:
                        break
                if len(lines) >= records:
                    break
    except (OSError, EOFError):
        pass
    return lines


def sample_headers(path: Path | str) -> list[str]:
    """Header lines carrying enough evidence to resolve a unit's read group."""
    return read_headers(path, INDEX_MAX_RECORDS, min_votes=INDEX_MIN_VOTES)


def sample_index(lines: list[str]) -> tuple[str | None, float]:
    """Consensus barcode over header lines, and its share of the clean ones.

    The share is measured against the records that yielded an unambiguous index,
    not against every record read, so it reports barcode agreement and not index
    read quality. It is 0.0 when nothing was usable.
    """
    indexes = [i for i in (parse_index(line) for line in lines) if i]
    if not indexes:
        return None, 0.0
    winner, hits = Counter(indexes).most_common(1)[0]
    return winner, hits / len(indexes)


def sample_lanes(lines: list[str]) -> set[str]:
    """Every distinct lane seen in the sampled headers.

    More than one is proof that the file spans lanes. One is not proof that it
    does not, because an externally merged file is a concatenation and its head
    is entirely the first lane -- which is why a lane is never taken from a read.
    """
    return {h.lane for h in (parse_header(line) for line in lines) if h}


def read_first_line(path: Path | str) -> str | None:
    """The first line of a gzipped FASTQ. None if unreadable or not a record."""
    lines = read_headers(path, 1)
    return lines[0] if lines else None
