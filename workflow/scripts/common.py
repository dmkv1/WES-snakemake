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


def is_paired_run(run):
    return runs_dict[run]["normal"] is not None


def is_purecn_eligible(wildcards):
    """Paired tumor/normal runs only (PDX included). Tumor-only runs have no
    matched-normal het-SNP track for PureCN's BAF-based purity/ploidy fit."""
    return is_paired_run(wildcards.run)


def _get_gender(run, sample):
    sample_row = samples[(samples["ID"] == run) & (samples["sample"] == sample)]
    return sample_row["gender"].iloc[0]


def get_purecn_normaldb(wildcards):
    probe_version = probe_dict[wildcards.run][wildcards.sample]
    sex_key = "normaldb_m" if _get_gender(wildcards.run, wildcards.sample) == "m" else "normaldb_f"
    return config["panel_of_normals"]["purecn"][sex_key][probe_version]


def get_purecn_mapping_bias(wildcards):
    probe_version = probe_dict[wildcards.run][wildcards.sample]
    sex_key = "mapping_bias_m" if _get_gender(wildcards.run, wildcards.sample) == "m" else "mapping_bias_f"
    return config["panel_of_normals"]["purecn"][sex_key][probe_version]


def get_purity_ploidy_args(wildcards, input):
    """Single source of truth for cnvkit call's --purity/--ploidy, read from
    resolve_purity_source's sidecar so cnvkit_call and combine_results can
    never disagree on which purity value was actually used."""
    import csv

    with open(input.purity_csv) as fh:
        row = next(csv.DictReader(fh))
    if row["source"] == "purecn":
        return f"--purity {row['purity']} --ploidy {row['ploidy']}"
    # No --ploidy -> cnvkit's own default (2), matching pre-PureCN behaviour.
    return f"--purity {row['purity']}"


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


def is_tumor_only(wildcards):
    """Check if run has no matched normal"""
    return runs_dict[wildcards.run]["normal"] is None


def get_pon_path(wildcards):
    """Get PON VCF path for tumor-only samples"""
    probe = probe_dict[wildcards.run][wildcards.sample]
    return config["panel_of_normals"]["mutect2"][probe]
