"""Resolve a samplesheet into alignment units and their read groups.

A *unit* is one expanded (fq1, fq2) pair -- the thing that gets aligned once,
with one read group. A samplesheet row yields one unit when it names a file pair
directly, and N units when its fq1/fq2 are globs covering a sample's lanes
(`..._S{n}_L*_R1_001.fastq.gz`). Both shapes are first-class.

Read groups are *derived* by default and *overridden* by samplesheet columns, so
the same code runs for a fully specified sheet and for the bare
`sample,fq1,fq2` minimum. Resolution, most confident first:

    flowcell = sheet.flowcell or header.flowcell or None
    lane     = sheet.lane     or filename _L(\\d{3})_ or None
    barcode  = sheet.barcode  or the consensus index over the first records

The lane is never taken from a read header. An externally lane-merged FASTQ
starts on one lane and ends on another -- a file in this cohort begins at
flowcell lane 1 and ends at lane 4 -- so the first record's lane is a claim
about the file that is simply false. The header's lane is still parsed: it flags
a merged or misnamed file when it disagrees with the filename, and more than one
lane across the sampled records proves the file spans lanes.

Failure is open: unresolvable provenance degrades to a positional unit token and
is reported through `rg_source`, rather than blocking a run. Operators who do
have curated metadata get the assertions instead, via strict=True.
"""

from __future__ import annotations

import glob
import os
import re
from typing import Any

import pandas as pd

from workflow.scripts.fastq_header import (
    MIXED_INDEX_THRESHOLD,
    parse_header,
    sample_headers,
    sample_index,
    sample_lanes,
)

# Lane token in a BaseSpace filename, e.g. ..._S3_L002_R1_001.fastq.gz.
FILENAME_LANE_RE = re.compile(r"_L(\d{3})_")

# An @RG line: tab-separated TAG:VALUE fields, tabs written as the two-character
# escape bwa expects. Printable ASCII only -- a brace would be eaten by
# Snakemake's shell formatting, a quote would break out of the shell argument.
RG_STRING_RE = re.compile(r"^@RG(\\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$")

# Columns that describe the sample, not the FASTQ pair. Every unit of a sample
# must agree on them; disagreement is an operator error, not something to
# silently resolve by picking a row.
SAMPLE_COLUMNS = ["sample_type", "gender", "capture_kit", "tumor_fraction"]

# Optional per-unit override columns. Absent in the current 8-column sheet.
UNIT_COLUMNS = ["flowcell", "lane", "barcode", "library"]

# Named provenance for the (flowcell source, lane source) pairs worth naming.
# Anything else is reported as "<flowcell source>+<lane source>" rather than
# being rounded to the nearest friendly label.
RG_SOURCE_LABELS = {
    ("sheet", "sheet"): "sheet",
    ("header", "filename"): "header+lane",
    ("header", "none"): "header_nolane",
    ("none", "filename"): "filename",
    ("none", "none"): "positional",
}

UNIT_TABLE_COLUMNS = [
    "ID", "sample", "unit", "fq1", "fq2", "flowcell", "lane", "barcode",
    "library", "rg_source", "rg_id", "rg_pu", "rg_string",
]


class SamplesheetError(ValueError):
    """A samplesheet that cannot be resolved into units."""


def _blank(value: Any) -> bool:
    return value is None or pd.isna(value) or str(value).strip() in ("", "NA")


def _text(value: Any) -> str | None:
    return None if _blank(value) else str(value).strip()


def _lane_token(lane: Any) -> str | None:
    """Normalise a lane to the L00N form, accepting 2, '2' or 'L002'."""
    if _blank(lane):
        return None
    raw = str(lane).strip()
    m = re.fullmatch(r"L?(\d{1,3})", raw)
    if not m:
        raise SamplesheetError(
            f"lane must be a number or an L00N token, got {raw!r}")
    return f"L{int(m.group(1)):03d}"


def expand_units(fq1_pattern: Any, fq2_pattern: Any, ident: str) -> list[dict]:
    """Expand one sheet row's fq1/fq2 into ordered (fq1, fq2, filename_lane).

    A plain path is a glob matching exactly one file, so single-file rows and
    glob rows go down the same path.
    """
    if _blank(fq1_pattern) or _blank(fq2_pattern):
        raise SamplesheetError(f"{ident}: fq1 and fq2 are both required")
    fq1s = sorted(glob.glob(str(fq1_pattern).strip()))
    fq2s = sorted(glob.glob(str(fq2_pattern).strip()))
    if not fq1s:
        raise SamplesheetError(
            f"{ident}: fq1 pattern matched no files: {fq1_pattern}")
    if len(fq1s) != len(fq2s):
        raise SamplesheetError(
            f"{ident}: R1/R2 file count mismatch ({len(fq1s)} vs {len(fq2s)}) "
            f"for patterns {fq1_pattern} / {fq2_pattern}")
    units = []
    for fq1, fq2 in zip(fq1s, fq2s):
        m = FILENAME_LANE_RE.search(os.path.basename(fq1))
        units.append({
            "fq1": fq1,
            "fq2": fq2,
            "filename_lane": f"L{int(m.group(1)):03d}" if m else None,
        })
    return units


def resolve_rg(unit: dict, row: dict, sample: str, n: int, *,
               resolve_headers: bool = True, strict: bool = False) -> dict:
    """Resolve one unit's read-group fields, provenance and warnings."""
    warnings: list[str] = []
    ident = f"{row.get('ID')}/{sample}"

    sheet_flowcell = _text(row.get("flowcell"))
    sheet_lane = _lane_token(row.get("lane"))
    sheet_barcode = _text(row.get("barcode"))
    sheet_library = _text(row.get("library"))

    header = index = None
    index_share = 0.0
    sampled_lanes: set[str] = set()
    if resolve_headers:
        lines = sample_headers(unit["fq1"])
        # The flowcell comes from the first record that parses, not strictly the
        # first record, so one malformed header does not blank the whole unit.
        header = next(
            (h for h in (parse_header(line) for line in lines) if h), None)
        index, index_share = sample_index(lines)
        sampled_lanes = sample_lanes(lines)

    header_flowcell = header.flowcell if header else None
    header_lane = header.lane if header else None
    filename_lane = unit["filename_lane"]

    def _disagree(message: str) -> None:
        if strict:
            raise SamplesheetError(f"{ident}: {message}")
        warnings.append(f"{ident}: {message}")

    if sheet_flowcell and header_flowcell and sheet_flowcell != header_flowcell:
        _disagree(f"samplesheet flowcell {sheet_flowcell} but reads carry "
                  f"{header_flowcell} ({os.path.basename(unit['fq1'])})")
    if sheet_lane and filename_lane and sheet_lane != filename_lane:
        _disagree(f"samplesheet lane {sheet_lane} but filename says "
                  f"{filename_lane} ({os.path.basename(unit['fq1'])})")
    if header_lane and filename_lane and header_lane != filename_lane:
        # Not an operator error: the file is lane-merged, or was renamed.
        _disagree(f"{os.path.basename(unit['fq1'])} is named {filename_lane} but "
                  f"its first read is {header_lane}; the file may span lanes")
    if len(sampled_lanes) > 1:
        _disagree(f"{os.path.basename(unit['fq1'])} carries lanes "
                  f"{'/'.join(sorted(sampled_lanes))} in its first records; it "
                  f"spans lanes and no single lane describes it")
    if index and index_share < MIXED_INDEX_THRESHOLD:
        # Not raised under strict: a majority barcode is still the best available
        # answer, and the file is usable. Demultiplexing is what needs looking at.
        warnings.append(
            f"{ident}: {os.path.basename(unit['fq1'])} mixes barcodes; "
            f"{index} is only {index_share:.0%} of the unambiguous records")

    flowcell = sheet_flowcell or header_flowcell
    lane = sheet_lane or filename_lane
    barcode = sheet_barcode or index

    flowcell_from = "sheet" if sheet_flowcell else ("header" if header_flowcell else "none")
    lane_from = "sheet" if sheet_lane else ("filename" if filename_lane else "none")
    rg_source = RG_SOURCE_LABELS.get(
        (flowcell_from, lane_from), f"{flowcell_from}+{lane_from}")

    return {
        "unit": lane if lane else f"u{n}",
        "position": n,
        "fq1": unit["fq1"],
        "fq2": unit["fq2"],
        "flowcell": flowcell or "",
        "lane": lane or "",
        "barcode": barcode or "",
        "library": sheet_library or sample,
        "rg_source": rg_source,
        "warnings": warnings,
    }


def _read_group_ids(row: dict) -> tuple[str, str]:
    """(ID, PU) for a resolved unit.

    They are not the same string and should not be. `PU` names the physical
    sequencing unit -- flowcell, lane, barcode -- which two samples multiplexed
    on one lane legitimately share. `ID` identifies the read group, must be
    unique, and is therefore keyed on the sample, whose names are unique across
    the cohort by construction.
    """
    rg_id = f"{row['sample']}.{row['unit']}"
    physical = ".".join(
        p for p in (row["flowcell"], row["lane"], row["barcode"]) if p)
    return rg_id, physical or rg_id


def build_rg_string(fields: dict) -> str:
    """The bwa `-R` argument, with the escaping asserted rather than hoped for.

    bwa's -R parser wants tabs written as the literal two characters `\\t` and
    rejects real tab characters. Metadata carrying a quote or a brace would
    break the shell argument or Snakemake's formatting, so the finished string
    is checked before it can reach a command line.
    """
    parts = [
        f"ID:{fields['rg_id']}",
        f"PU:{fields['rg_pu']}",
        f"SM:{fields['sample']}",
        f"LB:{fields['library']}",
        "PL:ILLUMINA",
    ]
    rg = "@RG\\t" + "\\t".join(parts)
    if "\t" in rg or "'" in rg or "{" in rg or "}" in rg:
        raise SamplesheetError(
            f"read group contains a character that cannot survive the shell: {rg!r}")
    if not RG_STRING_RE.match(rg):
        raise SamplesheetError(f"malformed read group string: {rg!r}")
    return rg


def _collapse_samples(units: pd.DataFrame, sheet_rows: dict) -> pd.DataFrame:
    """One validated row per (ID, sample), collapsed from its units.

    Every per-sample column must agree across a sample's rows. Resolving a
    disagreement by taking the first or last row is how a mislabelled sex or
    capture kit reaches the caller unnoticed.
    """
    records = []
    for (run, sample), group in units.groupby(["ID", "sample"], sort=False):
        record = {"ID": run, "sample": sample}
        for column in SAMPLE_COLUMNS:
            values = {
                _text(sheet_rows[idx].get(column))
                for idx in group["_row"].unique()
            }
            if len(values) > 1:
                shown = ", ".join(sorted(str(v) for v in values))
                raise SamplesheetError(
                    f"{run}/{sample}: rows disagree on {column} ({shown}). "
                    f"Every row of a sample must carry the same {column}.")
            record[column] = values.pop() if values else None
        libraries = sorted(set(group["library"]))
        record.update({
            "n_units": len(group),
            "units": ",".join(group["unit"]),
            "flowcells": ",".join(sorted({f for f in group["flowcell"] if f})),
            "library": ",".join(libraries),
            "rg_source": ",".join(sorted(set(group["rg_source"]))),
        })
        records.append(record)
    return pd.DataFrame.from_records(records)


def _validate_tumor_fraction(samples: pd.DataFrame) -> None:
    for _, row in samples.iterrows():
        value = row["tumor_fraction"]
        if _blank(value):
            continue
        try:
            tf = float(value)
        except (TypeError, ValueError):
            raise SamplesheetError(
                f"{row['sample']}: tumor_fraction must be a number in (0,1] or "
                f"NA/blank (got {value!r})") from None
        if not 0 < tf <= 1:
            raise SamplesheetError(
                f"{row['sample']}: tumor_fraction must be in (0,1] or NA/blank "
                f"(got {value!r}); 0 is not a valid tumor fraction")


def build_units(sheet: pd.DataFrame, *, resolve_headers: bool = True,
                strict: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Resolve a samplesheet into (units, samples_meta, warnings).

    `units` is one row per alignment unit, carrying its finished read group.
    `samples_meta` is one validated row per (ID, sample) -- the authority every
    per-sample consumer should read, in place of the raw samplesheet.
    """
    warnings: list[str] = []
    sheet_rows: dict[int, dict] = {}
    resolved: list[dict] = []

    for column in ("ID", "sample", "fq1", "fq2"):
        if column not in sheet.columns:
            raise SamplesheetError(f"samplesheet is missing required column {column!r}")

    per_sample_index: dict[tuple[str, str], int] = {}

    for idx, raw in sheet.iterrows():
        row = raw.to_dict()
        if _blank(row.get("ID")) or _blank(row.get("sample")):
            continue
        sheet_rows[idx] = row
        run, sample = str(row["ID"]).strip(), str(row["sample"]).strip()
        ident = f"{run}/{sample}"
        expanded = expand_units(row.get("fq1"), row.get("fq2"), ident)

        if len(expanded) > 1 and not _blank(row.get("lane")):
            raise SamplesheetError(
                f"{ident}: a row carrying an explicit lane must name one FASTQ "
                f"pair, but its pattern expanded to {len(expanded)} pairs. "
                f"Split it into one row per pair.")

        for unit in expanded:
            key = (run, sample)
            per_sample_index[key] = per_sample_index.get(key, 0) + 1
            fields = resolve_rg(unit, row, sample, per_sample_index[key],
                                resolve_headers=resolve_headers, strict=strict)
            warnings.extend(fields.pop("warnings"))
            fields.update({"ID": run, "sample": sample, "_row": idx})
            resolved.append(fields)

    if not resolved:
        raise SamplesheetError("samplesheet contains no usable rows")

    units = pd.DataFrame.from_records(resolved)
    warnings.extend(_ensure_unique_unit_tokens(units))
    warnings.extend(_warn_implicit_libraries(units, sheet_rows))

    # IDs are keyed on the unit token, so they are only assigned once the
    # tokens above are final.
    ids_and_pus = [_read_group_ids(row) for row in units.to_dict("records")]
    units["rg_id"] = [rg_id for rg_id, _ in ids_and_pus]
    units["rg_pu"] = [pu for _, pu in ids_and_pus]
    units["rg_string"] = [build_rg_string(row) for row in units.to_dict("records")]

    samples_meta = _collapse_samples(units, sheet_rows)
    _validate_tumor_fraction(samples_meta)

    units = units[UNIT_TABLE_COLUMNS].reset_index(drop=True)
    return units, samples_meta, warnings


def _ensure_unique_unit_tokens(units: pd.DataFrame) -> list[str]:
    """Force a sample's unit tokens apart, in place.

    The token is both the read group's discriminator and a filename component
    (`{sample}.{unit}_R1.fq.gz`), so a repeat inside one sample would collide
    two units' intermediates as well as merging their @RG records -- which GATK
    does without complaint, costing the distinction with nothing in the logs.
    Two lanes numbered L001 on different flowcells, and chunked `_001`/`_002`
    deliveries of one lane, both produce it.

    A sample whose lane tokens are not unique has lanes that do not identify its
    units, so the whole sample falls back to positional tokens rather than
    keeping a mix of real and invented ones.
    """
    warnings = []
    for (run, sample), group in units.groupby(["ID", "sample"], sort=False):
        tokens = list(group["unit"])
        if len(set(tokens)) == len(tokens):
            continue
        warnings.append(
            f"{run}/{sample}: unit tokens {sorted(set(tokens))} do not "
            f"distinguish its {len(tokens)} FASTQ pairs, so they are numbered "
            f"positionally instead. Set an explicit lane or library per row if "
            f"these are distinct sequencing units.")
        units.loc[group.index, "unit"] = [f"u{p}" for p in group["position"]]
    return warnings


def _warn_implicit_libraries(units: pd.DataFrame, sheet_rows: dict) -> list[str]:
    """Flag multi-unit samples that never said which library they came from.

    Defaulting LB to the sample name assumes one library. If two of those units
    are actually separate preps, MarkDuplicates treats their overlaps as PCR
    duplicates and discards real evidence -- and nothing in the data can say so.
    """
    warnings = []
    for (run, sample), group in units.groupby(["ID", "sample"], sort=False):
        if len(group) < 2:
            continue
        stated = any(
            not _blank(sheet_rows[idx].get("library"))
            for idx in group["_row"].unique()
        )
        if not stated:
            warnings.append(
                f"{run}/{sample}: {len(group)} units with no library column; "
                f"assuming one library (LB={sample}). If these are separate "
                f"preps, set library per row or duplicates will be overcalled.")
    return warnings


# MultiQC parses the '# '-prefixed preamble of a custom-content file as YAML.
# Apostrophes are avoided rather than escaped throughout: the text sits in a
# single-quoted YAML scalar inside a Python string, and one stray quote ends the
# scalar silently -- MultiQC then fails at the very end of a run.
MQC_DESCRIPTION = "Read group of each FASTQ pair, and the source of its flowcell and lane."

# The rg_source values name a pair of lookups: which source gave the flowcell,
# and which gave the lane. Kept as a column tooltip so the section header stays
# one line.
MQC_RG_SOURCE_TOOLTIP = (
    "sheet: the samplesheet. "
    "header+lane: flowcell from the reads, lane from the file name. "
    "header_nolane: flowcell from the reads, lane unknown, the file can span lanes. "
    "filename: lane from the file name only. "
    "positional: no source, units numbered in sequence."
)

MQC_COLUMNS = ["flowcell", "lane", "barcode", "library", "rg_source", "rg_id"]


def units_mqc_table(units: pd.DataFrame) -> str:
    """The units table as MultiQC custom content, header included."""
    header = [
        "# id: 'read_groups'",
        "# section_name: 'Read groups'",
        f"# description: '{MQC_DESCRIPTION}'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'read_groups_table'",
        "#     namespace: 'Read groups'",
        "# headers:",
        "#     rg_source:",
        f"#         description: '{MQC_RG_SOURCE_TOOLTIP}'",
    ]
    rows = ["Sample\t" + "\t".join(MQC_COLUMNS)]
    for row in units.to_dict("records"):
        rows.append(f"{row['sample']}.{row['unit']}\t"
                    + "\t".join(str(row[c]) for c in MQC_COLUMNS))
    return "\n".join(header + rows) + "\n"


def sample_renames(units: pd.DataFrame) -> dict[str, str]:
    """MultiQC row renames: drop the unit token from single-unit samples.

    The token stands in for a lane. A sample with one FASTQ pair has no lane to
    distinguish, so `{sample}.u1` says nothing and only splits that sample into
    two rows, one carrying its FASTQ metrics and one its BAM metrics.

    Done per sample rather than by stripping `.u<N>` everywhere: a sample that
    genuinely has several positional units needs them kept apart, or their rows
    would overwrite each other and lose data silently.
    """
    renames = {}
    for (_run, sample), group in units.groupby(["ID", "sample"], sort=False):
        if len(group) == 1:
            renames[f"{sample}.{group.iloc[0]['unit']}"] = sample
    return renames


ROW_TYPE_TOOLTIP = (
    "sample: metrics from the finished BAM. unit: one FASTQ pair. "
    "read: one FASTQ file. A sample with one FASTQ pair has its sample and unit "
    "rows merged."
)


def row_types(units: pd.DataFrame) -> dict[str, str]:
    """{MultiQC row name: sample|unit|read}.

    The General Statistics table holds three kinds of row, each with its own
    disjoint set of columns, because BAM metrics are per sample while FASTQ
    metrics are per unit and per read. This labels which is which so the table
    can be sorted into those groups.

    Names are the ones MultiQC ends up using, so sample_renames is applied
    first: a single-unit sample has no separate unit row to label.
    """
    renames = sample_renames(units)
    rows: dict[str, str] = {}
    for (_run, sample), group in units.groupby(["ID", "sample"], sort=False):
        rows[sample] = "sample"
        for unit in group["unit"]:
            base = renames.get(f"{sample}.{unit}", f"{sample}.{unit}")
            if base != sample:
                rows[base] = "unit"
            for read in ("R1", "R2"):
                rows[f"{base}_{read}"] = "read"
    return rows


def row_types_mqc_table(units: pd.DataFrame) -> str:
    """The row-type labels as MultiQC general-stats custom content."""
    header = [
        "# id: 'row_type'",
        "# plot_type: 'generalstats'",
        "# pconfig:",
        "#     - row_type:",
        "#         title: 'Row'",
        f"#         description: '{ROW_TYPE_TOOLTIP}'",
    ]
    rows = ["Sample\trow_type"]
    for name, kind in row_types(units).items():
        rows.append(f"{name}\t{kind}")
    return "\n".join(header + rows) + "\n"


def units_by_sample(units: pd.DataFrame) -> dict[tuple[str, str], list[str]]:
    """{(run, sample): [unit token, ...]} in samplesheet order."""
    out: dict[tuple[str, str], list[str]] = {}
    for row in units.itertuples():
        out.setdefault((row.ID, row.sample), []).append(row.unit)
    return out


def unit_index(units: pd.DataFrame) -> dict[tuple[str, str, str], dict]:
    """{(run, sample, unit): unit row} for input and params functions."""
    return {
        (row["ID"], row["sample"], row["unit"]): row
        for row in units.to_dict("records")
    }
