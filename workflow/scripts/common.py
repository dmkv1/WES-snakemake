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


def get_unit_display_name(wildcards):
    """The name this unit's QC rows carry in the MultiQC report.

    Mirrors units.sample_renames, which does the same collapse for rows MultiQC
    names after a filename. This one is for rows named from log content, which
    that renaming never reaches.
    """
    if len(get_units(wildcards.run, wildcards.sample)) == 1:
        return wildcards.sample
    return f"{wildcards.sample}.{wildcards.unit}"


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


def get_known_ploidy(wildcards):
    """Orthogonal/measured ploidy from the samplesheet (known_ploidy column),
    or None when unknown. Ground truth for cnvkit ploidy (karyotype, flow DNA
    index); when None, resolve_purity_source falls back to PureCN, then to
    diploid."""
    return known_ploidy_dict[wildcards.run][wildcards.sample]


# Below this, PureCN's own fit ("LOW PURITY" comment) counts as unreliable
# enough to flag downstream, but is still used — see resolve_purity below.
PURITY_CONFIDENCE_THRESHOLD = 0.30


def resolve_purity(known, purecn_purity, purecn_failed, purecn_available):
    """Purity + source + confidence for one sample.

    known: orthogonal ground-truth purity (float in (0,1]) or None.
    purecn_purity: PureCN's raw Purity field ("", "NA", None, or numeric string).
    purecn_failed: PureCN's own Failed flag.
    purecn_available: PureCN enabled, eligible, and its output was produced.

    Priority: known > purecn (any numeric fit, regardless of flag content —
    a flagged purity is still the best estimate we have; forcing purity 1
    instead only erases real CNA signal) > assumed_pure (purity 1, used only
    when no measurement exists at all).

    purity_confidence flags samples for cautious interpretation without
    suppressing them: 'unknown' for assumed_pure (no real measurement),
    'low_purity' for any resolved purity below PURITY_CONFIDENCE_THRESHOLD,
    else 'high'.
    """
    if known is not None:
        purity, source = str(known), "known"
    elif purecn_available and not purecn_failed and purecn_purity not in ("", "NA", None):
        purity, source = str(purecn_purity), "purecn"
    else:
        purity, source = "1", "assumed_pure"

    if source == "assumed_pure":
        confidence = "unknown"
    elif float(purity) < PURITY_CONFIDENCE_THRESHOLD:
        confidence = "low_purity"
    else:
        confidence = "high"

    return purity, source, confidence


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
