import pandas as pd


configfile: "config.yaml"


samples = pd.read_excel("mock_fastq/fastq.groups.xlsx")

# Create a dictionary of runs and their samples
runs_dict = {}
for run in samples["run"].unique():
    if pd.isna(run):
        continue
    run_samples = samples[samples["run"] == run]
    runs_dict[run] = {
        "normal": run_samples[run_samples["Type"] == "CTRL"]["samplename"].iloc[0],
        "tumors": run_samples[run_samples["Type"] != "CTRL"]["samplename"].tolist(),
    }

# Create a dictionary of PDX samples
pdx_dict = {}
for run in samples["run"].unique():
    if pd.isna(run):
        continue
    run_samples = samples[samples["run"] == run]
    pdx_samples = run_samples[run_samples["Type"] == "PDX"]["samplename"].tolist()
    if pdx_samples:  # Only add the run if it has PDX samples
        pdx_dict[run] = pdx_samples

# Create a dictionary to store fastq paths
fastq_dict = {}
for _, row in samples.iterrows():
    if pd.notna(row["run"]):
        if row["run"] not in fastq_dict:
            fastq_dict[row["run"]] = {}
        fastq_dict[row["run"]][row["samplename"]] = {
            "fq1": row["fq1"],
            "fq2": row["fq2"],
        }


wildcard_constraints:
    sample="[^.]+",  # Match anything except dots
    run="[^.]+",  # Match anything except dots


# Import helper functions
from workflow.scripts.common import *

# Make dictionaries available to the common module
import sys
import workflow.scripts.common as common

common.fastq_dict = fastq_dict
common.runs_dict = runs_dict
common.pdx_dict = pdx_dict
common.config = config


# Rules
include: "workflow/rules/ref_index.smk"
include: "workflow/rules/bam_mapping_gatk.smk"
include: "workflow/rules/xenofilter.smk"
include: "workflow/rules/snv_calling_mutect2.smk"
include: "workflow/rules/snv_calling_varscan2.smk"
include: "workflow/rules/snv_calling_strelka2.smk"


rule all:
    input:
        [
            f"bam/{run}/{sample}.bam"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
        [f"vcf/{run}/{run}.mutect2.vcf" for run in runs_dict],
        [f"vcf/{run}/{run}.varscan.vcf" for run in runs_dict],
        [
            f"vcf/{run}/{sample}.strelka.vcf"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
