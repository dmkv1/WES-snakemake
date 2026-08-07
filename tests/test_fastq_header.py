"""Header parsing across the shapes a FASTQ can actually arrive in."""

from __future__ import annotations

from workflow.scripts.fastq_header import (
    INDEX_MAX_RECORDS,
    INDEX_MIN_VOTES,
    parse_header,
    parse_index,
    read_first_header,
    read_first_line,
    sample_headers,
    sample_index,
    sample_lanes,
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


def test_sampling_recovers_a_barcode_the_first_record_lost(fx, tmp_path):
    """The fix: one miscalled index base must not blank the unit's barcode.

    Two of three files sampled from this cohort carry an N in the very first
    record, so reading only that record threw away a barcode the file states
    plainly from the second record on.
    """
    path = tmp_path / "dirty_head_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(20, index=fx.CLEAN_INDEX,
                                       head_index=fx.DIRTY_INDEX))
    assert parse_index(read_first_line(path)) is None

    index, share = sample_index(sample_headers(path))
    assert index == fx.CLEAN_INDEX
    assert share == 1.0


def test_every_index_ambiguous_yields_no_barcode(fx, tmp_path):
    """Sampling widens the evidence; it does not invent a barcode."""
    path = tmp_path / "all_dirty_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(20, index=fx.DIRTY_INDEX))
    assert sample_index(sample_headers(path)) == (None, 0.0)


def test_mixed_barcodes_report_the_majority_and_its_share(fx, tmp_path):
    path = tmp_path / "mixed_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(20, index=fx.CLEAN_INDEX,
                                       minority_index=fx.OTHER_INDEX,
                                       minority_share=0.25))
    index, share = sample_index(sample_headers(path))
    assert index == fx.CLEAN_INDEX
    assert share == 0.75


def test_sampling_stops_once_the_votes_are_in(fx, tmp_path):
    """A clean file costs the minimum: reading stops at the vote target."""
    path = tmp_path / "long_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(INDEX_MIN_VOTES + 500))
    assert len(sample_headers(path)) == INDEX_MIN_VOTES


def test_sampling_is_capped_when_the_votes_never_arrive(fx, tmp_path):
    """A file that never yields a clean index must still terminate."""
    path = tmp_path / "dirty_long_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(INDEX_MAX_RECORDS + 100,
                                       index=fx.DIRTY_INDEX))
    assert len(sample_headers(path)) == INDEX_MAX_RECORDS


def test_a_dirty_head_does_not_end_the_sample(fx, tmp_path):
    """The P020 case: no clean index until past record 100.

    The head of a run is its worst part, so a fixed window can land entirely
    inside the bad stretch. Reading is bounded by votes collected, not records
    read, which is what carries the sample past it.
    """
    path = tmp_path / "bad_head_R1_001.fastq.gz"
    fx.write_fastq(path, fx.fastq_text(INDEX_MIN_VOTES + 300,
                                       index=fx.CLEAN_INDEX,
                                       head_index=fx.DIRTY_INDEX,
                                       head_records=150))
    index, share = sample_index(sample_headers(path))
    assert index == fx.CLEAN_INDEX
    assert share == 1.0


def test_sample_headers_returns_only_headers(fx, tmp_path):
    fq1, _ = fx.pair(tmp_path, "S1", lane=3)
    lines = sample_headers(fq1)
    assert lines and all(line.startswith("@") for line in lines)
    assert all(parse_header(line).lane == "L003" for line in lines)


def test_sample_headers_degrade_on_unreadable_input(fx, tmp_path):
    assert sample_headers(tmp_path / "absent.fastq.gz") == []
    assert sample_headers(fx.write_plain(tmp_path / "plain_R1_001.fastq.gz")) == []
    # Truncation keeps whatever survived rather than discarding the file.
    lines = sample_headers(fx.write_truncated(tmp_path / "cut_R1_001.fastq.gz"))
    assert all(line.startswith("@") for line in lines)


def test_sampled_lanes_prove_a_file_spans_lanes(fx, tmp_path):
    """More than one lane in the window is proof; one lane is not proof of one.

    An externally merged file is a concatenation, so its head is entirely the
    first lane and the sample cannot see the change. That asymmetry is why the
    lane is still never taken from a read.
    """
    spanning = tmp_path / "interleaved_R1_001.fastq.gz"
    fx.write_fastq(spanning, fx.fastq_text(4, lane=1, final_lane=4))
    assert sample_lanes(sample_headers(spanning)) == {"L001", "L004"}

    concatenated = tmp_path / "merged_R1_001.fastq.gz"
    fx.write_fastq(concatenated,
                   fx.fastq_text(500, lane=1) + fx.fastq_text(500, lane=4))
    assert sample_lanes(sample_headers(concatenated)) == {"L001"}


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
