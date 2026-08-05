# Getter functions
def get_fastq1(wildcards):
    return unit_index[(wildcards.run, wildcards.sample, wildcards.unit)]["fq1"]


def get_fastq2(wildcards):
    return unit_index[(wildcards.run, wildcards.sample, wildcards.unit)]["fq2"]


def get_units(run, sample):
    """Ordered unit tokens for a sample (['L001', 'L002', ...] or ['u1', ...])."""
    return units_by_sample[(run, sample)]


def get_unit_bams(wildcards):
    """Per-unit query-grouped BAMs gathered by MarkDuplicates."""
    return [
        f"results/{wildcards.run}/{wildcards.sample}/bam/units/"
        f"{wildcards.sample}.{unit}.qgrp.bam"
        for unit in get_units(wildcards.run, wildcards.sample)
    ]


def get_read_group(wildcards):
    """The unit's finished bwa '-R' argument.

    Resolved and validated once at load time by units.build_units, so nothing
    here parses a FASTQ or builds a string that could reach a shell malformed.
    """
    return unit_index[(wildcards.run, wildcards.sample, wildcards.unit)]["rg_string"]


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
    # `samples` is the collapsed per-sample table, which build_units validated
    # as one row per (ID, sample), so this match is unique by construction.
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


def is_tumor_only(wildcards):
    """Check if run has no matched normal"""
    return runs_dict[wildcards.run]["normal"] is None


def get_pon_path(wildcards):
    """Get PON VCF path for tumor-only samples"""
    probe = probe_dict[wildcards.run][wildcards.sample]
    return config["panel_of_normals"]["mutect2"][probe]
