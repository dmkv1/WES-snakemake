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


rule units_multiqc_table:
    """MultiQC custom-content view of the unit table."""
    input:
        tsv="results/metadata/units.tsv",
    output:
        tsv="results/qc/units_rg_mqc.tsv",
    run:
        columns = [
            "ID", "sample", "unit", "flowcell", "lane", "barcode", "library",
            "rg_source", "rg_id",
        ]
        with open(output.tsv, "w") as fh:
            fh.write("# id: 'read_groups'\n")
            fh.write("# section_name: 'Read groups'\n")
            fh.write(
                "# description: 'How each FASTQ pair's read group was resolved. "
                "rg_source is the provenance: sheet (samplesheet asserted it), "
                "header+lane (flowcell from the reads, lane from the filename), "
                "header_nolane (flowcell known, lane unknown -- the file may span "
                "lanes), filename (lane from the filename only), positional "
                "(neither; units numbered in order).'\n"
            )
            fh.write("# plot_type: 'table'\n")
            fh.write("# pconfig:\n")
            fh.write("#     id: 'read_groups_table'\n")
            fh.write("#     namespace: 'Read groups'\n")
            fh.write("Sample\t" + "\t".join(columns[2:]) + "\n")
            for row in units.to_dict("records"):
                name = f"{row['sample']}.{row['unit']}"
                fh.write(name + "\t" + "\t".join(
                    str(row[c]) for c in columns[2:]) + "\n")
