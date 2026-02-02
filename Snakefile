import pandas as pd
import os


configfile: "config.yaml"


samples = pd.read_csv(config["samplesheet"])

valid_sex_values = {"XX", "XY"}
invalid_sex = samples[~samples["Chr_sex"].isin(valid_sex_values) & samples["Chr_sex"].notna()]
if not invalid_sex.empty:
    invalid_samples = invalid_sex[["ID", "sample", "Chr_sex"]].to_string(index=False)
    raise ValueError(
        f"Invalid Chr_sex values found. Must be 'XX' or 'XY':\n{invalid_samples}"
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
            if probe not in config.get("tumor_only", {}).get("pon", {}):
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


# Rules
include: "workflow/rules/ref_index.smk"
include: "workflow/rules/host_read_filter.smk"
include: "workflow/rules/bam_mapping_gatk.smk"
include: "workflow/rules/snv_calling_mutect2.smk"
include: "workflow/rules/sv_calling_manta.smk"
include: "workflow/rules/cnv_calling_cnvkit.smk"
include: "workflow/rules/integrate_results.smk"

rule all:
    input:
        "results/qc/multiqc_report.html",
        # Table with results
        [
            f"results/{run}/{sample}/{sample}_results.xlsx"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        # SNV/Indel VCFs (Mutect2)
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
