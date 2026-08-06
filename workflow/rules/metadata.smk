# Provenance for how each FASTQ pair's read group was decided.
#
# Read groups are resolved once at load time (workflow/scripts/units.py) from
# the samplesheet and the reads themselves. These rules serialise that decision
# so it is inspectable without running an aligner, and so a non-expert reading
# the MultiQC report can see which rung of the ladder each unit landed on --
# `header_nolane` in particular means the file may be externally lane-merged and
# its read group is per file rather than per lane.
#
# Written by rules rather than at load time so nothing is emitted during a
# --dry-run, and so editing the samplesheet retriggers the consumers.

rule write_unit_table:
    """Every alignment unit, including the exact @RG line bwa will be given."""
    input:
        samplesheet=config["samplesheet"],
    output:
        tsv="results/metadata/units.tsv",
    run:
        units.to_csv(output.tsv, sep="\t", index=False)


rule write_sample_table:
    """The collapsed per-sample authority.

    One row per (ID, sample), validated at load time. Every consumer that wants
    per-sample metadata reads this instead of the samplesheet, which has one row
    per FASTQ pair and would fan results out by unit count on a join.
    """
    input:
        samplesheet=config["samplesheet"],
    output:
        tsv="results/metadata/samples.tsv",
    run:
        samples.to_csv(output.tsv, sep="\t", index=False)


rule row_types_multiqc_table:
    """Label each General Statistics row as sample, unit or read.

    That table carries three kinds of row with disjoint columns, because BAM
    metrics are per sample and FASTQ metrics are per unit and per read. The
    label makes the table sortable into those groups.
    """
    input:
        samplesheet=config["samplesheet"],
    output:
        tsv="results/qc/row_type_mqc.tsv",
    run:
        with open(output.tsv, "w") as fh:
            fh.write(row_types_mqc_table(units))


rule multiqc_renames:
    """Collapse single-unit samples back onto one MultiQC row.

    Generated per cohort rather than written into multiqc_config.yaml, because
    which samples have one unit depends on the samplesheet.
    """
    input:
        samplesheet=config["samplesheet"],
    output:
        yaml="results/qc/multiqc_renames.yaml",
    run:
        import yaml as _yaml

        with open(output.yaml, "w") as fh:
            _yaml.safe_dump({"sample_names_replace": sample_renames(units)}, fh)


rule units_multiqc_table:
    """MultiQC custom-content view of the unit table."""
    input:
        tsv="results/metadata/units.tsv",
    output:
        tsv="results/qc/units_rg_mqc.tsv",
    run:
        # Built by units.units_mqc_table so the YAML preamble is covered by the
        # test suite: a malformed one only surfaces when MultiQC runs, which is
        # the last rule of the workflow.
        with open(output.tsv, "w") as fh:
            fh.write(units_mqc_table(units))
