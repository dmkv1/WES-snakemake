import pandas as pd
import os
import re
import sys


configfile: "config.yaml"

bind_paths = ",".join([config["refs"]["path"], config["panel_of_normals"]["path"]])
os.environ["APPTAINER_BIND"] = bind_paths
os.environ["SINGULARITY_BIND"] = bind_paths

sys.path.insert(0, os.path.dirname(workflow.snakefile))
from workflow.scripts.units import (
    build_units,
    unit_index as _build_unit_index,
    units_by_sample as _build_units_by_sample,
)

sheet = pd.read_csv(config["samplesheet"])

# Run IDs feed the `{run}_{sample}` file-naming scheme (metrics, logs). An ID
# containing '_', '.' or '/' would make that split ambiguous, so reject it up
# front with a clear error instead of failing silently downstream.
bad_ids = sorted(
    {str(x) for x in sheet["ID"].dropna().unique() if re.search(r"[_./]", str(x))}
)
if bad_ids:
    raise ValueError(
        "Run ID(s) contain '_', '.' or '/', which collide with the "
        f"'{{run}}_{{sample}}' file naming: {', '.join(bad_ids)}"
    )

# Resolve the samplesheet into alignment units and a validated per-sample table.
# A row naming a file pair yields one unit; a row whose fq1/fq2 are globs (what
# wesingest emits) yields one per pair. Read groups are derived from the data
# and overridden by the optional flowcell/lane/library/barcode columns -- see
# workflow/scripts/units.py for the resolution ladder.
units, samples, rg_warnings = build_units(
    sheet, strict=config.get("params", {}).get("rg", {}).get("strict", True)
)
for _warning in rg_warnings:
    print(f"WARNING [read groups] {_warning}", file=sys.stderr)

unit_index = _build_unit_index(units)
units_by_sample = _build_units_by_sample(units)

valid_gender_values = {"f", "m"}
invalid_gender = samples[~samples["gender"].isin(valid_gender_values) & samples["gender"].notna()]
if not invalid_gender.empty:
    invalid_samples = invalid_gender[["ID", "sample", "gender"]].to_string(index=False)
    raise ValueError(
        f"Invalid gender values found. Must be 'f' or 'm':\n{invalid_samples}"
    )

# All four dicts below are built from the collapsed per-sample table, so a
# multi-unit sample contributes exactly one entry. Building them from the raw
# samplesheet would duplicate every multi-unit tumor through rule all and the
# result-combining rules, and would resolve conflicting per-sample columns by
# silently keeping whichever row happened to be last.

# {'ID': {'normal': 'CTRL', 'tumors': ['PT', 'PDX']}}
# For tumor-only runs, normal will be None
runs_dict = {}
for run in samples["ID"].unique():
    run_samples = samples[samples["ID"] == run]
    ctrl_samples = run_samples[run_samples["sample_type"] == "CTRL"]["sample"].tolist()
    tumor_samples = run_samples[run_samples["sample_type"] != "CTRL"]["sample"].tolist()

    runs_dict[run] = {
        "normal": ctrl_samples[0] if ctrl_samples else None,
        "tumors": tumor_samples,
    }

# {'ID': ['PDX']}
pdx_dict = {}
for run in samples["ID"].unique():
    run_samples = samples[samples["ID"] == run]
    pdx_samples = run_samples[run_samples["sample_type"] == "PDX"]["sample"].tolist()
    if pdx_samples:
        pdx_dict[run] = pdx_samples

probe_dict = {}
for _, row in samples.iterrows():
    if pd.notna(row["capture_kit"]):
        probe_dict.setdefault(row["ID"], {})[row["sample"]] = row["capture_kit"]

# tumor_fraction: orthogonal/measured tumor cell fraction (sorting, cytometry,
# PDX=1), used as ground truth for cnvkit purity when present. Blank/NA means
# unknown -> resolve_purity_source falls back to PureCN (see purecn.smk). Stored
# as float or None; build_units has already validated the range.
tumor_fraction_dict = {}
for _, row in samples.iterrows():
    val = row["tumor_fraction"]
    tf = None if pd.isna(val) or str(val).strip() in ("", "NA") else float(val)
    tumor_fraction_dict.setdefault(row["ID"], {})[row["sample"]] = tf

# Validate tumor-only runs have PON configured
for run, data in runs_dict.items():
    if data["normal"] is None:
        for sample in data["tumors"]:
            probe = probe_dict[run][sample]
            if probe not in config.get("panel_of_normals", {}).get("mutect2", {}):
                raise ValueError(f"No PON configured for probe '{probe}' (tumor-only run '{run}')")


wildcard_constraints:
    run="[^/._]+",  # Match anything except slashes, dots and underscores
    sample="[^/.]+",  # Match anything except slashes and dots
    # Unit token: a real lane (L001) where one is known, else positional (u1).
    # '.' separates it from the sample in {sample}.{unit} paths.
    unit="L[0-9]{3}|u[0-9]+",


# Import helper functions
from workflow.scripts.common import *

# Make dictionaries available to the common module
import workflow.scripts.common as common

common.units = units
common.unit_index = unit_index
common.units_by_sample = units_by_sample
common.runs_dict = runs_dict
common.pdx_dict = pdx_dict
common.probe_dict = probe_dict
common.tumor_fraction_dict = tumor_fraction_dict
common.config = config
common.samples = samples


# Rules
include: "workflow/rules/metadata.smk"
include: "workflow/rules/ref_index.smk"
include: "workflow/rules/host_read_filter.smk"
include: "workflow/rules/bam_mapping_gatk.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/somalier.smk"
include: "workflow/rules/snv_calling_mutect2.smk"
include: "workflow/rules/sv_calling_manta.smk"
include: "workflow/rules/sv_calling_delly.smk"
include: "workflow/rules/sv_merge.smk"
include: "workflow/rules/cnv_calling_cnvkit.smk"
include: "workflow/rules/purecn.smk"
include: "workflow/rules/integrate_results.smk"

def _tg_notify(msg):
    env = config.get("telegram_bot_env", "")
    if not env:
        return
    run_id = config.get("run_id", "")
    prefix = f"[{run_id}] " if run_id else ""
    shell(f"source {env} && "
          f"curl -s -d \"chat_id=$TELEGRAM_CHAT_ID&text={prefix}{msg}\" "
          f"\"https://api.telegram.org/bot$TELEGRAM_BOT_TOKEN/sendMessage\" > /dev/null")

onstart:
    _tg_notify("WES-snakemake started 🚀")

onsuccess:
    _tg_notify("WES-snakemake finished successfully ✅")

onerror:
    _tg_notify("WES-snakemake FAILED ❌ — check snakemake.log")


rule all:
    input:
        "results/qc/multiqc_report.html",
        # Table with results
        [
            f"results/{run}/{sample}/{sample}_results.xlsx"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        # Combined cross-sample tables for downstream analysis (Wesseract etc.)
        "results/combined/combined_snvs.tsv",
        "results/combined/combined_cnvs.tsv",
        "results/combined/combined_svs.tsv",
        "results/combined/combined_qc.tsv",
        # Cohort MAF for maftools (oncoplot/TMB/etc.)
        "results/combined/cohort.maf",
        # SNV/Indel VCFs (Mutect2 + VEP)
        [
            f"results/{run}/{sample}/{sample}.SNV.vcf"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        # CNVkit plots
        [
            f"results/{run}/{sample}/{sample}.scatter.png"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        [
            f"results/{run}/{sample}/{sample}.diagram.pdf"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        # PureCN purity/ploidy estimates — always built for paired tumors so
        # they're available for report comparison regardless of whether
        # params.cnv.purity_source actually uses them (see resolve_purity_source).
        [
            f"work/purecn/{run}/{sample}/{sample}.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
            if is_paired_run(run)
        ],
        # Somalier sample-swap / sex QC: cohort all-vs-all relate (HTML +
        # MultiQC), annotated relatedness table, and sample x sample heatmap.
        "results/qc/somalier/cohort/cohort.html",
        "results/combined/combined_relatedness.tsv",
        "results/combined/relatedness_heatmap.png",
