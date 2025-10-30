import gzip
from typing import Dict


# Getter functions
def get_fastq1(wildcards):
    run = wildcards.run
    return fastq_dict[run][wildcards.sample]["fq1"]


def get_fastq2(wildcards):
    run = wildcards.run
    return fastq_dict[run][wildcards.sample]["fq2"]


def get_tumor_bams(wildcards):
    return expand("results/bam/{sample}.bam", sample=runs_dict[wildcards.run]["tumors"])


def get_library_prep_for_sample(sample_name):
    for run, samples in probe_dict.items():
        if sample_name in samples:
            probe_version = samples[sample_name]
            return config["probe_configs"][probe_version]["library_prep"]
    raise ValueError(f"Sample {sample_name} not found in probe_dict")


def get_probe_version(wildcards):
    return probe_dict[wildcards.run][wildcards.sample]

def get_purity(wildcards):
    return purity_dict[wildcards.run][wildcards.sample]


def parse_fastq_header(fastq_path: str, sample_name: str) -> Dict[str, str]:
    # Read first line of the FASTQ file
    with (
        gzip.open(fastq_path, "rt")
        if fastq_path.endswith(".gz")
        else open(fastq_path, "r")
    ) as f:
        header = f.readline().strip()

    # Parse the instrument string (part before the space)
    try:
        parts = header.lstrip("@").split()[0].split(":")
        rgid = f"{parts[0]}:{parts[1]}"
        platform_unit = parts[2]
    except (IndexError, AttributeError):
        raise ValueError(f"Could not parse header in {fastq_path}: {header}")

    return {
        "RGID": rgid,
        "RGPU": platform_unit,
        "RGSM": sample_name,
        "RGPL": "ILLUMINA",
        "RGLB": f"{sample_name}_{get_library_prep_for_sample(sample_name)}",
    }


def get_read_group_params(wildcards) -> Dict[str, str]:
    fq1_path = fastq_dict[wildcards.run][wildcards.sample]["fq1"]
    return parse_fastq_header(fq1_path, wildcards.sample)
