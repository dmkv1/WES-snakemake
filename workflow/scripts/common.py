import gzip
from typing import Dict


# Getter functions
def get_fastq1(wildcards):
    return fastq_dict[wildcards.run][wildcards.sample]["fq"][wildcards.lane]["fq1"]


def get_fastq2(wildcards):
    return fastq_dict[wildcards.run][wildcards.sample]["fq"][wildcards.lane]["fq2"]


def get_lanes(run, sample):
    """Ordered lane tokens for a sample (['L001', 'L002', ...])."""
    return fastq_dict[run][sample]["lanes"]


def get_lane_bams(wildcards):
    """Per-lane coordinate-sorted BAMs gathered by MarkDuplicates."""
    lanes = fastq_dict[wildcards.run][wildcards.sample]["lanes"]
    return [
        f"results/{wildcards.run}/{wildcards.sample}/bam/lanes/"
        f"{wildcards.sample}.{lane}.sorted.bam"
        for lane in lanes
    ]


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

def get_known_purity(wildcards):
    """Orthogonal/measured tumor fraction from the samplesheet (tumor_fraction
    column), or None when unknown. Ground truth for cnvkit purity; when None,
    resolve_purity_source falls back to PureCN."""
    return tumor_fraction_dict[wildcards.run][wildcards.sample]


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
    # resolve_purity_source always writes an integer ploidy (PureCN's rounded
    # estimate, else diploid), so pass it regardless of the purity source.
    return f"--purity {row['purity']} --ploidy {row['ploidy']}"


def get_lane_read_group(wildcards) -> str:
    """Build a per-lane bwa '-R' @RG line from the lane's FASTQ header.

    Illumina header: @instrument:run:flowcell:lane:tile:x:y ...
    ID/PU are flowcell.lane (unique per lane -> no collision across lanes of a
    multi-lane sample, TODO #15). LB is per-sample (lanes share one library, so
    MarkDuplicates still detects cross-lane PCR duplicates). Tabs are the literal
    two-character escape '\\t': bwa's -R parser requires escaped tabs and rejects
    real <tab> characters, expanding '\\t' to tabs itself.
    """
    fq1 = fastq_dict[wildcards.run][wildcards.sample]["fq"][wildcards.lane]["fq1"]
    with (
        gzip.open(fq1, "rt") if str(fq1).endswith(".gz") else open(fq1, "r")
    ) as f:
        header = f.readline().strip()

    try:
        parts = header.lstrip("@").split()[0].split(":")
        flowcell, lane = parts[2], parts[3]
    except IndexError:
        raise ValueError(f"Could not parse flowcell/lane from header in {fq1}: {header}")

    sample = wildcards.sample
    library = f"{sample}_{get_library_prep_for_sample(sample)}"
    fields = [
        f"ID:{flowcell}.{lane}",
        f"PU:{flowcell}.{lane}",
        f"SM:{sample}",
        f"LB:{library}",
        "PL:ILLUMINA",
    ]
    return "@RG\\t" + "\\t".join(fields)


def is_tumor_only(wildcards):
    """Check if run has no matched normal"""
    return runs_dict[wildcards.run]["normal"] is None


def get_pon_path(wildcards):
    """Get PON VCF path for tumor-only samples"""
    probe = probe_dict[wildcards.run][wildcards.sample]
    return config["panel_of_normals"]["mutect2"][probe]
