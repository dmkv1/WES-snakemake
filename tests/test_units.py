"""Read-group resolution: the ladder, the guards, and the collapse."""

from __future__ import annotations

import pandas as pd
import pytest

from workflow.scripts.units import (
    row_types,
    row_types_mqc_table,
    ROW_TYPE_TOOLTIP,
    sample_renames,
    units_mqc_table,
    SamplesheetError,
    build_rg_string,
    build_units,
    expand_units,
    unit_index,
    units_by_sample,
)


def _sheet(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def _resolve(rows: list[dict], **kwargs):
    return build_units(_sheet(rows), **kwargs)


# --- expansion -------------------------------------------------------------

def test_plain_path_is_one_unit(fx, tmp_path):
    """The existing 8-column sheet shape must keep working unchanged."""
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    units, samples, warnings = _resolve([fx.sheet_row("S1", fq1, fq2)])
    assert len(units) == 1
    assert len(samples) == 1
    assert samples.iloc[0]["n_units"] == 1


def test_glob_row_expands_to_all_lanes(fx, tmp_path):
    """wesingest emits glob rows, so a row yields N units."""
    for lane in (1, 2, 3, 4):
        fx.pair(tmp_path, "S1", lane=lane)
    units, samples, _ = _resolve([fx.sheet_row(
        "S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
        tmp_path / "S1_S1_L*_R2_001.fastq.gz")])
    assert list(units["unit"]) == ["L001", "L002", "L003", "L004"]
    assert samples.iloc[0]["n_units"] == 4


def test_r1_r2_count_mismatch_is_refused(fx, tmp_path):
    fx.pair(tmp_path, "S1", lane=1)
    (tmp_path / "S1_S1_L002_R1_001.fastq.gz").write_bytes(b"")
    with pytest.raises(SamplesheetError, match="count mismatch"):
        expand_units(tmp_path / "S1_S1_L*_R1_001.fastq.gz",
                     tmp_path / "S1_S1_L*_R2_001.fastq.gz", "P001/S1")


def test_empty_glob_is_refused(tmp_path):
    with pytest.raises(SamplesheetError, match="matched no files"):
        expand_units(tmp_path / "nope_R1.fastq.gz",
                     tmp_path / "nope_R2.fastq.gz", "P001/S1")


# --- the ladder ------------------------------------------------------------

def test_header_flowcell_plus_filename_lane(fx, tmp_path):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=2)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    row = units.iloc[0]
    assert row["rg_source"] == "header+lane"
    assert row["unit"] == "L002"
    assert row["flowcell"] == fx.FLOWCELL
    # ID identifies the read group and is keyed on the sample; PU names the
    # physical sequencing unit, which two multiplexed samples can share.
    assert row["rg_id"] == "S1.L002"
    assert row["rg_pu"] == f"{fx.FLOWCELL}.L002.{fx.CLEAN_INDEX}"


def test_lane_is_never_taken_from_the_header(fx, tmp_path):
    """A file with no lane in its name resolves to a positional unit token.

    Its reads say lane 3, but an externally merged file's first record says
    nothing about the rest of the file, so that claim is not propagated.
    """
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=3, lane_in_name=False)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    row = units.iloc[0]
    assert row["rg_source"] == "header_nolane"
    assert row["unit"] == "u1"
    assert row["lane"] == ""
    assert "L003" not in row["rg_id"]
    assert "L003" not in row["rg_pu"]


def test_filename_lane_without_a_parseable_header(fx, tmp_path):
    path1 = tmp_path / "S1_S1_L002_R1_001.fastq.gz"
    path2 = tmp_path / "S1_S1_L002_R2_001.fastq.gz"
    fx.write_non_casava(path1)
    fx.write_non_casava(path2)
    units, _, _ = _resolve([fx.sheet_row("S1", path1, path2)])
    row = units.iloc[0]
    assert row["rg_source"] == "filename"
    assert row["unit"] == "L002"
    assert row["rg_id"] == "S1.L002"
    assert row["rg_pu"] == "L002"


def test_no_provenance_at_all_is_positional(fx, tmp_path):
    path1 = tmp_path / "anon_R1.fastq.gz"
    path2 = tmp_path / "anon_R2.fastq.gz"
    fx.write_non_casava(path1)
    fx.write_non_casava(path2)
    units, _, _ = _resolve([fx.sheet_row("S1", path1, path2)])
    row = units.iloc[0]
    assert row["rg_source"] == "positional"
    assert row["unit"] == "u1"
    assert row["rg_id"] == "S1.u1"


def test_sheet_columns_override_the_data(fx, tmp_path):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1, lane_in_name=False)
    units, _, _ = _resolve([fx.sheet_row(
        "S1", fq1, fq2, flowcell=fx.FLOWCELL, lane="L007", library="lib9")])
    row = units.iloc[0]
    assert row["rg_source"] == "sheet"
    assert row["unit"] == "L007"
    assert row["library"] == "lib9"


@pytest.mark.parametrize("value,expected", [("2", "L002"), (2, "L002"), ("L002", "L002")])
def test_lane_column_accepts_bare_numbers(fx, tmp_path, value, expected):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=2, lane_in_name=False)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2, lane=value)])
    assert units.iloc[0]["unit"] == expected


def test_malformed_lane_column_is_refused(fx, tmp_path):
    fq1, fq2 = fx.pair(tmp_path, "S1")
    with pytest.raises(SamplesheetError, match="lane must be"):
        _resolve([fx.sheet_row("S1", fq1, fq2, lane="lane-two")])


# --- barcode ---------------------------------------------------------------

def test_ambiguous_barcode_is_dropped(fx, tmp_path):
    """Every record ambiguous means no barcode, not a guessed one."""
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1, index=fx.DIRTY_INDEX)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    row = units.iloc[0]
    assert row["barcode"] == ""
    assert row["rg_pu"] == f"{fx.FLOWCELL}.L001"


def test_barcode_survives_a_miscall_in_the_first_record(fx, tmp_path):
    """The whole point of sampling: it must reach PU, not just the parser.

    Lanes 3 and 4 of the cohort's real run lost their barcode this way, so their
    PU read `<flowcell>.L003` while lanes 1 and 2 carried the index.
    """
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1, n=20, index=fx.CLEAN_INDEX,
                       head_index=fx.DIRTY_INDEX)
    units, _, warnings = _resolve([fx.sheet_row("S1", fq1, fq2)])
    row = units.iloc[0]
    assert row["barcode"] == fx.CLEAN_INDEX
    assert row["rg_pu"] == f"{fx.FLOWCELL}.L001.{fx.CLEAN_INDEX}"
    assert warnings == []


def test_mixed_barcodes_warn_but_do_not_block(fx, tmp_path):
    """A file that was not cleanly demultiplexed still resolves to a majority.

    This stays a warning under strict: the majority barcode is the best answer
    available and the unit is usable, so blocking the run would be the wrong
    trade. What needs attention is the demultiplexing, not the pipeline.
    """
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1, n=20, index=fx.CLEAN_INDEX,
                       minority_index=fx.OTHER_INDEX, minority_share=0.4)
    units, _, warnings = _resolve([fx.sheet_row("S1", fq1, fq2)], strict=True)
    assert units.iloc[0]["barcode"] == fx.CLEAN_INDEX
    assert any("mixes barcodes" in w for w in warnings)


def test_sheet_barcode_overrides_the_consensus(fx, tmp_path):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1, index=fx.CLEAN_INDEX)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2, barcode="ACGTACGT")])
    assert units.iloc[0]["barcode"] == "ACGTACGT"


# --- cross-checks ----------------------------------------------------------

def test_lanes_within_the_sample_prove_a_span(fx, tmp_path):
    """Two lanes in the sampled window: the filename lane is a lie."""
    fq1 = tmp_path / "S1_S1_L001_R1_001.fastq.gz"
    fq2 = tmp_path / "S1_S1_L001_R2_001.fastq.gz"
    fx.write_fastq(fq1, fx.fastq_text(8, lane=1, final_lane=4, read="R1"))
    fx.write_fastq(fq2, fx.fastq_text(8, lane=1, final_lane=4, read="R2"))
    _, _, warnings = _resolve([fx.sheet_row("S1", fq1, fq2)])
    assert any("spans lanes" in w for w in warnings)
    with pytest.raises(SamplesheetError, match="spans lanes"):
        _resolve([fx.sheet_row("S1", fq1, fq2)], strict=True)

def test_header_lane_disagreeing_with_filename_warns(fx, tmp_path):
    """The signal that a file is lane-merged or was renamed."""
    fq1 = tmp_path / "S1_S1_L001_R1_001.fastq.gz"
    fq2 = tmp_path / "S1_S1_L001_R2_001.fastq.gz"
    fx.write_fastq(fq1, fx.fastq_text(lane=4, read="R1"))
    fx.write_fastq(fq2, fx.fastq_text(lane=4, read="R2"))
    _, _, warnings = _resolve([fx.sheet_row("S1", fq1, fq2)])
    assert any("may span lanes" in w for w in warnings)


def test_disagreement_is_fatal_under_strict(fx, tmp_path):
    fq1 = tmp_path / "S1_S1_L001_R1_001.fastq.gz"
    fq2 = tmp_path / "S1_S1_L001_R2_001.fastq.gz"
    fx.write_fastq(fq1, fx.fastq_text(lane=4, read="R1"))
    fx.write_fastq(fq2, fx.fastq_text(lane=4, read="R2"))
    with pytest.raises(SamplesheetError, match="may span lanes"):
        _resolve([fx.sheet_row("S1", fq1, fq2)], strict=True)


def test_sheet_flowcell_contradicting_the_reads_warns(fx, tmp_path):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    _, _, warnings = _resolve([fx.sheet_row("S1", fq1, fq2, flowcell="WRONGFC")])
    assert any("but reads carry" in w for w in warnings)


def test_explicit_lane_on_a_glob_row_is_refused(fx, tmp_path):
    for lane in (1, 2):
        fx.pair(tmp_path, "S1", lane=lane)
    with pytest.raises(SamplesheetError, match="expanded to 2 pairs"):
        _resolve([fx.sheet_row(
            "S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
            tmp_path / "S1_S1_L*_R2_001.fastq.gz", lane="L001")])


# --- guards ----------------------------------------------------------------

def test_repeated_unit_token_falls_back_to_positional(fx, tmp_path):
    """Chunked deliveries: two files claiming the same lane of one sample.

    The token names an intermediate file as well as a read group, so a repeat
    inside one sample would collide both.
    """
    fq1a, fq2a = fx.pair(tmp_path / "a", "S1", lane=1)
    fq1b, fq2b = fx.pair(tmp_path / "b", "S1", lane=1)
    units, _, warnings = _resolve([
        fx.sheet_row("S1", fq1a, fq2a),
        fx.sheet_row("S1", fq1b, fq2b),
    ])
    assert list(units["unit"]) == ["u1", "u2"]
    assert len(set(units["rg_id"])) == 2
    assert any("do not distinguish" in w for w in warnings)


def test_different_samples_may_share_a_flowcell(fx, tmp_path):
    """Multiplexed samples share a PU but must never share an RG ID."""
    fq1a, fq2a = fx.pair(tmp_path / "a", "S1", lane=1, lane_in_name=False)
    fq1b, fq2b = fx.pair(tmp_path / "b", "S2", lane=1, lane_in_name=False)
    units, _, warnings = _resolve([
        fx.sheet_row("S1", fq1a, fq2a),
        fx.sheet_row("S2", fq1b, fq2b),
    ])
    assert list(units["rg_id"]) == ["S1.u1", "S2.u1"]
    assert not warnings


def test_multi_unit_sample_without_library_warns(fx, tmp_path):
    for lane in (1, 2):
        fx.pair(tmp_path, "S1", lane=lane)
    _, _, warnings = _resolve([fx.sheet_row(
        "S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
        tmp_path / "S1_S1_L*_R2_001.fastq.gz")])
    assert any("assuming one library" in w for w in warnings)


def test_library_column_silences_the_warning(fx, tmp_path):
    fq1a, fq2a = fx.pair(tmp_path / "a", "S1", lane=1)
    fq1b, fq2b = fx.pair(tmp_path / "b", "S1", lane=2)
    _, samples, warnings = _resolve([
        fx.sheet_row("S1", fq1a, fq2a, library="lib1"),
        fx.sheet_row("S1", fq1b, fq2b, library="lib2"),
    ])
    assert not any("assuming one library" in w for w in warnings)
    assert samples.iloc[0]["library"] == "lib1,lib2"


def test_rg_string_is_shell_safe(fx, tmp_path):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    rg = units.iloc[0]["rg_string"]
    assert rg.startswith("@RG\\t")
    assert "\t" not in rg
    assert "SM:S1" in rg and "LB:S1" in rg and "PL:ILLUMINA" in rg


def test_rg_string_rejects_a_quote():
    with pytest.raises(SamplesheetError, match="cannot survive the shell"):
        build_rg_string({"rg_id": "S1.L001", "rg_pu": "FC.L001",
                         "sample": "S'1", "library": "lib"})


# --- collapse --------------------------------------------------------------

def test_disagreeing_gender_across_rows_is_refused(fx, tmp_path):
    fq1a, fq2a = fx.pair(tmp_path / "a", "S1", lane=1)
    fq1b, fq2b = fx.pair(tmp_path / "b", "S1", lane=2)
    with pytest.raises(SamplesheetError, match="disagree on gender"):
        _resolve([
            fx.sheet_row("S1", fq1a, fq2a, gender="f"),
            fx.sheet_row("S1", fq1b, fq2b, gender="m"),
        ])


def test_disagreeing_capture_kit_is_refused(fx, tmp_path):
    fq1a, fq2a = fx.pair(tmp_path / "a", "S1", lane=1)
    fq1b, fq2b = fx.pair(tmp_path / "b", "S1", lane=2)
    with pytest.raises(SamplesheetError, match="disagree on capture_kit"):
        _resolve([
            fx.sheet_row("S1", fq1a, fq2a, capture_kit="V6+UTR"),
            fx.sheet_row("S1", fq1b, fq2b, capture_kit="V8+UTR"),
        ])


def test_samples_meta_has_one_row_per_sample(fx, tmp_path):
    """The property combine_all.R's left_join depends on."""
    for lane in (1, 2, 3, 4):
        fx.pair(tmp_path, "S1", lane=lane)
    fq1, fq2 = fx.pair(tmp_path / "ctrl", "S1_CTRL", lane=1)
    _, samples, _ = _resolve([
        fx.sheet_row("S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
                     tmp_path / "S1_S1_L*_R2_001.fastq.gz"),
        fx.sheet_row("S1_CTRL", fq1, fq2, sample_type="CTRL"),
    ])
    assert len(samples) == 2
    assert not samples.duplicated(subset=["ID", "sample"]).any()
    assert "fq1" not in samples.columns


@pytest.mark.parametrize("value", ["0", "1.5", "-0.2", "abc"])
def test_invalid_tumor_fraction_is_refused(fx, tmp_path, value):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    with pytest.raises(SamplesheetError, match="tumor_fraction"):
        _resolve([fx.sheet_row("S1", fq1, fq2, tumor_fraction=value)])


@pytest.mark.parametrize("value", ["", "NA", "0.4", "1"])
def test_valid_tumor_fraction_is_accepted(fx, tmp_path, value):
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    _resolve([fx.sheet_row("S1", fq1, fq2, tumor_fraction=value)])


def test_missing_required_column(fx, tmp_path):
    with pytest.raises(SamplesheetError, match="missing required column"):
        build_units(pd.DataFrame([{"ID": "P001", "sample": "S1"}]))


# --- lookups ---------------------------------------------------------------

def test_single_unit_sample_drops_its_unit_token(fx, tmp_path):
    """A lone positional token says nothing and splits the sample in two rows."""
    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1, lane_in_name=False)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    assert sample_renames(units) == {"S1.u1": "S1"}


def test_multi_unit_sample_keeps_its_tokens(fx, tmp_path):
    """Collapsing these would overwrite one unit's row with another."""
    for lane in (1, 2):
        fx.pair(tmp_path, "S1", lane=lane)
    units, _, _ = _resolve([fx.sheet_row(
        "S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
        tmp_path / "S1_S1_L*_R2_001.fastq.gz")])
    assert sample_renames(units) == {}


def test_renames_are_per_sample(fx, tmp_path):
    """One sample collapsing must not collapse another that has real units."""
    fq1, fq2 = fx.pair(tmp_path / "solo", "S1", lane=1, lane_in_name=False)
    for lane in (1, 2):
        fx.pair(tmp_path / "multi", "S2", lane=lane)
    units, _, _ = _resolve([
        fx.sheet_row("S1", fq1, fq2),
        fx.sheet_row("S2", tmp_path / "multi/S2_S1_L*_R1_001.fastq.gz",
                     tmp_path / "multi/S2_S1_L*_R2_001.fastq.gz"),
    ])
    assert sample_renames(units) == {"S1.u1": "S1"}


def test_row_types_label_every_general_stats_row(fx, tmp_path):
    """One sample with four units, one with a single pair."""
    for lane in (1, 2, 3, 4):
        fx.pair(tmp_path / "multi", "S1", lane=lane)
    fq1, fq2 = fx.pair(tmp_path / "solo", "S2", lane=1, lane_in_name=False)
    units, _, _ = _resolve([
        fx.sheet_row("S1", tmp_path / "multi/S1_S1_L*_R1_001.fastq.gz",
                     tmp_path / "multi/S1_S1_L*_R2_001.fastq.gz"),
        fx.sheet_row("S2", fq1, fq2),
    ])
    rt = row_types(units)
    assert rt["S1"] == "sample"
    assert rt["S1.L001"] == "unit"
    assert rt["S1.L001_R1"] == "read"
    # The single-unit sample has no separate unit row: it merged onto the sample.
    assert rt["S2"] == "sample"
    assert "S2.u1" not in rt
    assert rt["S2_R1"] == "read"
    from collections import Counter
    assert Counter(rt.values()) == {"read": 10, "unit": 4, "sample": 2}


def test_row_type_header_is_valid_yaml(fx, tmp_path):
    import yaml

    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    text = row_types_mqc_table(units)
    header = "\n".join(l[2:] for l in text.splitlines() if l.startswith("# "))
    parsed = yaml.safe_load(header)
    assert parsed["plot_type"] == "generalstats"


def test_mqc_text_has_no_apostrophe():
    """Both sit in single-quoted YAML scalars; one apostrophe ends the scalar."""
    from workflow.scripts.units import MQC_DESCRIPTION, MQC_RG_SOURCE_TOOLTIP

    assert "'" not in MQC_DESCRIPTION
    assert "'" not in MQC_RG_SOURCE_TOOLTIP
    assert "'" not in ROW_TYPE_TOOLTIP


def test_mqc_header_is_valid_yaml(fx, tmp_path):
    """MultiQC parses the '# ' preamble as YAML.

    A stray apostrophe ends the single-quoted scalar and MultiQC fails -- on the
    last rule of the workflow, hours in. Cheap to pin here instead.
    """
    import yaml

    fq1, fq2 = fx.pair(tmp_path, "S1", lane=1)
    units, _, _ = _resolve([fx.sheet_row("S1", fq1, fq2)])
    text = units_mqc_table(units)

    header = "\n".join(
        line[2:] for line in text.splitlines() if line.startswith("# ")
    )
    parsed = yaml.safe_load(header)
    assert parsed["id"] == "read_groups"
    assert parsed["plot_type"] == "table"
    assert parsed["pconfig"]["id"] == "read_groups_table"
    # The rg_source legend is a column tooltip, so the section header stays one line.
    assert "positional" in parsed["headers"]["rg_source"]["description"]


def test_mqc_table_rows_match_units(fx, tmp_path):
    for lane in (1, 2):
        fx.pair(tmp_path, "S1", lane=lane)
    units, _, _ = _resolve([fx.sheet_row(
        "S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
        tmp_path / "S1_S1_L*_R2_001.fastq.gz")])
    body = [l for l in units_mqc_table(units).splitlines() if not l.startswith("#")]
    assert body[0].startswith("Sample\t")
    assert [r.split("\t")[0] for r in body[1:]] == ["S1.L001", "S1.L002"]


def test_lookup_helpers(fx, tmp_path):
    for lane in (1, 2):
        fx.pair(tmp_path, "S1", lane=lane)
    units, _, _ = _resolve([fx.sheet_row(
        "S1", tmp_path / "S1_S1_L*_R1_001.fastq.gz",
        tmp_path / "S1_S1_L*_R2_001.fastq.gz")])
    assert units_by_sample(units) == {("P001", "S1"): ["L001", "L002"]}
    assert unit_index(units)[("P001", "S1", "L002")]["unit"] == "L002"
