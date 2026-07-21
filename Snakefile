import pandas as pd
import os
import re


configfile: "config.yaml"

bind_paths = ",".join([config["refs"]["path"], config["panel_of_normals"]["path"]])
os.environ["APPTAINER_BIND"] = bind_paths
os.environ["SINGULARITY_BIND"] = bind_paths

samples = pd.read_csv(config["samplesheet"])

# Run IDs feed the `{run}_{sample}` file-naming scheme (metrics, logs). An ID
# containing '_', '.' or '/' would make that split ambiguous, so reject it up
# front with a clear error instead of failing silently downstream.
bad_ids = sorted(
    {str(x) for x in samples["ID"].dropna().unique() if re.search(r"[_./]", str(x))}
)
if bad_ids:
    raise ValueError(
        "Run ID(s) contain '_', '.' or '/', which collide with the "
        f"'{{run}}_{{sample}}' file naming: {', '.join(bad_ids)}"
    )

valid_gender_values = {"f", "m"}
invalid_gender = samples[~samples["gender"].isin(valid_gender_values) & samples["gender"].notna()]
if not invalid_gender.empty:
    invalid_samples = invalid_gender[["ID", "sample", "gender"]].to_string(index=False)
    raise ValueError(
        f"Invalid gender values found. Must be 'f' or 'm':\n{invalid_samples}"
    )

# Validate fastq files exist
missing_files = []
for _, row in samples.iterrows():
    for fq_col in ["fq1", "fq2"]:
        fq_path = row[fq_col]
        if pd.notna(fq_path) and not os.path.exists(fq_path):
            missing_files.append(f"  {row['ID']}/{row['sample']}: {fq_path}")
if missing_files:
    raise FileNotFoundError(
        f"FASTQ files not found:\n" + "\n".join(missing_files)
    )

# Create a dictionary of runs and their samples
# {'ID': {'normal': 'CTRL', 'tumors': ['PT', 'PDX']}}
# For tumor-only runs, normal will be None
runs_dict = {}
for run in samples["ID"].unique():
    if pd.isna(run):
        continue
    run_samples = samples[samples["ID"] == run]
    ctrl_samples = run_samples[run_samples["sample_type"] == "CTRL"]["sample"].tolist()
    tumor_samples = run_samples[run_samples["sample_type"] != "CTRL"]["sample"].tolist()

    runs_dict[run] = {
        "normal": ctrl_samples[0] if ctrl_samples else None,
        "tumors": tumor_samples,
    }

# Create a dictionary of PDX samples
# {'ID': ['PDX']}
pdx_dict = {}
for run in samples["ID"].unique():
    if pd.isna(run):
        continue
    run_samples = samples[samples["ID"] == run]
    pdx_samples = run_samples[run_samples["sample_type"] == "PDX"]["sample"].tolist()
    if pdx_samples:
        pdx_dict[run] = pdx_samples

# Create a dictionary to store fastq paths
fastq_dict = {}
for _, row in samples.iterrows():
    if pd.notna(row["ID"]):
        if row["ID"] not in fastq_dict:
            fastq_dict[row["ID"]] = {}
        fastq_dict[row["ID"]][row["sample"]] = {
            "fq1": row["fq1"],
            "fq2": row["fq2"],
        }

# Create probe configuration dictionary
probe_dict = {}
for _, row in samples.iterrows():
    if pd.notna(row["ID"]) and pd.notna(row["probes"]):
        if row["ID"] not in probe_dict:
            probe_dict[row["ID"]] = {}
        probe_dict[row["ID"]][row["sample"]] = row["probes"]

purity_dict = {}
for _, row in samples.iterrows():
    if pd.notna(row["ID"]):
        if row["ID"] not in purity_dict:
            purity_dict[row["ID"]] = {}
        purity_dict[row["ID"]][row["sample"]] = row["purity"]

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


# Import helper functions
from workflow.scripts.common import *

# Make dictionaries available to the common module
import sys
import workflow.scripts.common as common

common.fastq_dict = fastq_dict
common.runs_dict = runs_dict
common.pdx_dict = pdx_dict
common.probe_dict = probe_dict
common.purity_dict = purity_dict
common.config = config
common.samples = samples


# Rules
include: "workflow/rules/ref_index.smk"
include: "workflow/rules/host_read_filter.smk"
include: "workflow/rules/bam_mapping_gatk.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/somalier.smk"
include: "workflow/rules/snv_calling_mutect2.smk"
include: "workflow/rules/sv_calling_manta.smk"
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
        # Somalier sample-swap / sex QC report
        [f"results/qc/somalier/{run}/{run}.html" for run in runs_dict],
