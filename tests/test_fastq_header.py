"""Header parsing across the shapes a FASTQ can actually arrive in."""

from __future__ import annotations

from workflow.scripts.fastq_header import (
    parse_header,
    parse_index,
    read_first_header,
    read_first_line,
)


def test_parses_casava(fx, tmp_path):
    fq1, _ = fx.pair(tmp_path, "S1", lane=2)
    header = parse_header(read_first_line(fq1))
    assert header.flowcell == fx.FLOWCELL
    assert header.lane == "L002"
    assert header.instrument == fx.INSTRUMENT


def test_clean_index_is_returned(fx, tmp_path):
    fq1, _ = fx.pair(tmp_path, "S1", index=fx.CLEAN_INDEX)
    assert parse_index(read_first_line(fq1)) == fx.CLEAN_INDEX


def test_index_with_n_is_rejected(fx, tmp_path):
    """A miscalled base in one record must not reach every read group."""
    fq1, _ = fx.pair(tmp_path, "S1", index=fx.DIRTY_INDEX)
    assert parse_index(read_first_line(fq1)) is None


def test_non_casava_header_is_unparsed(fx, tmp_path):
    path = fx.write_non_casava(tmp_path / "sra_R1_001.fastq.gz")
    line = read_first_line(path)
    assert line is not None
    assert parse_header(line) is None
    assert parse_index(line) is None


def test_truncated_gzip_does_not_raise(fx, tmp_path):
    """A cut payload must degrade to a value, never an exception.

    The first record usually survives truncation because gzip decompresses
    incrementally, so this pins the contract rather than the exact result.
    """
    path = fx.write_truncated(tmp_path / "cut_R1_001.fastq.gz")
    line = read_first_line(path)
    assert line is None or line.startswith("@")
    header = read_first_header(path)
    assert header is None or header.flowcell == fx.FLOWCELL


def test_plain_text_masquerading_as_gzip(fx, tmp_path):
    path = fx.write_plain(tmp_path / "plain_R1_001.fastq.gz")
    assert read_first_line(path) is None
    assert read_first_header(path) is None


def test_missing_file(tmp_path):
    assert read_first_line(tmp_path / "absent.fastq.gz") is None
    assert read_first_header(tmp_path / "absent.fastq.gz") is None


def test_first_record_hides_a_lane_change(fx, tmp_path):
    """The merged-file case: the first record's lane is not the file's lane.

    This is why the resolution ladder never takes a lane from a read header.
    """
    path = tmp_path / "merged_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(6, lane=1, final_lane=4))
    assert parse_header(read_first_line(path)).lane == "L001"

    import gzip
    with gzip.open(path, "rt") as fh:
        lanes = {ln.split(":")[3] for ln in fh if ln.startswith("@")}
    assert lanes == {"1", "4"}
