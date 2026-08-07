"""Synthesise FASTQ files and samplesheets for the test suite.

Everything is built under a caller-supplied temp directory. Fixtures are
generated rather than committed because .gitignore excludes *.fastq.gz and
*.csv, so checked-in fixture data would silently never reach the repository.

FASTQ payloads are real four-line records carrying real Illumina headers, since
the read group is resolved by parsing them. The header variants reproduced here
are the ones observed in this cohort or expected from ported deployments:

    clean           CASAVA 1.8+, one lane, clean barcode
    dirty_barcode   CASAVA 1.8+, N in the index (2 of 3 sampled files have this)
    dirty_head      N in the index of the leading records only, clean after
    mixed_barcode   two barcodes in one file, i.e. not cleanly demultiplexed
    lane_drift      starts on one lane, ends on another -- an externally merged
                    file, where the first record's lane is a false claim
    non_casava      SRA-style, carries no flowcell or lane
    truncated       valid gzip header, payload cut mid-stream
"""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

INSTRUMENT = "NB501049"
RUN_NUMBER = 744
FLOWCELL = "HHFT3AFXC"
READ_LENGTH = 151
CLEAN_INDEX = "AGTCGCGA+CTGGTCTA"
DIRTY_INDEX = "AGTCGCGA+NTGGTCTA"
OTHER_INDEX = "TTGACCTG+CTGGTCTA"

SHEET_COLUMNS = [
    "ID", "sample", "sample_type", "gender", "capture_kit", "tumor_fraction",
    "fq1", "fq2",
]


def _record(rng: random.Random, *, lane: int, read: str, flowcell: str,
            instrument: str, run_number: int, index: str, length: int) -> str:
    """One four-line CASAVA 1.8+ record."""
    member = read[-1] if read[-1] in "12" else "1"
    tile = rng.choice([11101, 11102, 21101])
    x, y = rng.randint(1000, 30000), rng.randint(1000, 30000)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    header = (f"@{instrument}:{run_number}:{flowcell}:{lane}:{tile}:{x}:{y} "
              f"{member}:N:0:{index}")
    return f"{header}\n{seq}\n+\n{'I' * length}\n"


def fastq_text(n: int = 4, *, lane: int = 1, read: str = "R1",
               flowcell: str = FLOWCELL, instrument: str = INSTRUMENT,
               run_number: int = RUN_NUMBER, index: str = CLEAN_INDEX,
               length: int = READ_LENGTH, final_lane: int | None = None,
               head_index: str | None = None, head_records: int = 1,
               minority_index: str | None = None,
               minority_share: float = 0.0) -> str:
    """`n` records. `final_lane` puts the last record on a different lane.

    That last-record case is the externally-merged file: parsing only the first
    record reports `lane`, while the file actually spans `lane`..`final_lane`.

    `head_index` overrides the index of the leading `head_records` records, which
    builds the file whose first record carries a miscalled index while the rest
    of the unit is clean. `minority_index` gives that share of the records a
    second barcode, for the mixed-barcode case.
    """
    rng = random.Random(f"{instrument}{run_number}{flowcell}{lane}{read}")
    head = head_records if head_index is not None else 0
    minority = round(n * minority_share) if minority_index else 0
    out = []
    for i in range(n):
        rec_lane = final_lane if (final_lane and i == n - 1) else lane
        if i < head:
            rec_index = head_index
        elif i < head + minority:
            rec_index = minority_index
        else:
            rec_index = index
        out.append(_record(rng, lane=rec_lane, read=read, flowcell=flowcell,
                           instrument=instrument, run_number=run_number,
                           index=rec_index, length=length))
    return "".join(out)


def write_fastq(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def write_non_casava(path: Path, n: int = 4) -> Path:
    """Gzipped FASTQ whose headers carry no flowcell or lane (SRA-style)."""
    rng = random.Random("sra")
    out = []
    for i in range(n):
        seq = "".join(rng.choice("ACGT") for _ in range(36))
        out.append(f"@SRR001666.{i + 1} length=36\n{seq}\n+\n{'I' * 36}\n")
    return write_fastq(path, "".join(out))


def write_truncated(path: Path) -> Path:
    """Valid gzip magic, payload cut mid-stream -- read_first_header must not raise."""
    write_fastq(path, fastq_text(8))
    raw = path.read_bytes()
    path.write_bytes(raw[: len(raw) // 2])
    return path


def write_plain(path: Path) -> Path:
    """Not gzipped at all, despite the .gz suffix."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(fastq_text(2))
    return path


def pair(directory: Path, sample: str, *, lane: int | None = 1,
         lane_in_name: bool = True, snum: int = 1, **kwargs) -> tuple[Path, Path]:
    """Write an R1/R2 pair, returning both paths.

    `lane` is the lane written into the read headers. `lane_in_name` controls
    whether the filename carries the `_L00N_` token -- the two are deliberately
    separable so tests can build a file whose name and content disagree.
    """
    token = f"_L{lane:03d}" if (lane_in_name and lane is not None) else ""
    paths = []
    for read in ("R1", "R2"):
        p = directory / f"{sample}_S{snum}{token}_{read}_001.fastq.gz"
        write_fastq(p, fastq_text(read=read, lane=lane or 1, **kwargs))
        paths.append(p)
    return paths[0], paths[1]


def write_sheet(path: Path, rows: list[dict], columns: list[str] | None = None) -> Path:
    """A samplesheet CSV. Columns default to the current 8-column schema."""
    if columns is None:
        columns = list(SHEET_COLUMNS)
        for row in rows:
            for key in row:
                if key not in columns:
                    columns.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def sheet_row(sample: str, fq1: Path | str, fq2: Path | str, *, run: str = "P001",
              sample_type: str = "PT", gender: str = "f",
              capture_kit: str = "V6+UTR", tumor_fraction: str = "",
              **extra) -> dict:
    """One samplesheet row with the per-sample columns filled in."""
    row = {
        "ID": run, "sample": sample, "sample_type": sample_type,
        "gender": gender, "capture_kit": capture_kit,
        "tumor_fraction": tumor_fraction, "fq1": str(fq1), "fq2": str(fq2),
    }
    row.update(extra)
    return row
