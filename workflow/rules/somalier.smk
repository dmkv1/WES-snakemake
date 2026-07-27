# Genotype-based sample-swap and sex QC (TODO #20). Report-only: results feed
# MultiQC's native somalier module for human review, no DAG gating on mismatch.


rule somalier_extract:
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        refg=config["refs"]["genome_human"],
        sites=config["refs"]["somalier_sites"],
    output:
        "work/somalier/{run}/{sample}.somalier",
    params:
        out_dir=lambda w: f"work/somalier/{w.run}",
    log:
        "work/logs/somalier_extract_{run}_{sample}.log",
    conda:
        "../envs/somalier.yaml"
    shell:
        """
        somalier extract \
            -d {params.out_dir} \
            --sites {input.sites} \
            -f {input.refg} \
            {input.bam} \
            > {log} 2>&1
        """


def all_somalier_extracts():
    """Every sample's .somalier across all runs (flat list)."""
    return [
        f"work/somalier/{run}/{s}.somalier"
        for run in runs_dict
        for s in get_all_samples_for_run(run)
    ]


rule somalier_ped_cohort:
    """Cohort PED: concat of every run's samples (family_id=run) so the
    all-vs-all `relate` can label within-run vs cross-run pairs and still flag
    sex against the samplesheet gender."""
    output:
        ped="work/somalier/cohort.ped",
    run:
        sex_code = {"m": "1", "f": "2"}
        with open(output.ped, "w") as fh:
            for run in runs_dict:
                for sample in get_all_samples_for_run(run):
                    gender = samples.loc[samples["sample"] == sample, "gender"].iloc[0]
                    fh.write(f"{run}\t{sample}\t0\t0\t{sex_code.get(gender, '0')}\t0\n")


rule somalier_relate_cohort:
    """All-vs-all relate over every sample in the cohort in one pass. Cross-run
    pairs (invisible to the per-run relate) surface patient-swap errors; this is
    also the pairs table MultiQC ingests for its relatedness scatter."""
    input:
        extracted=all_somalier_extracts(),
        ped="work/somalier/cohort.ped",
    output:
        html="results/qc/somalier/cohort/cohort.html",
        samples="results/qc/somalier/cohort/cohort.samples.tsv",
        pairs="results/qc/somalier/cohort/cohort.pairs.tsv",
    params:
        prefix="results/qc/somalier/cohort/cohort",
    log:
        "work/logs/somalier_relate_cohort.log",
    conda:
        "../envs/somalier.yaml"
    shell:
        """
        somalier relate \
            -o {params.prefix} \
            --ped {input.ped} \
            {input.extracted} \
            > {log} 2>&1
        """


rule somalier_summarise:
    """Annotate the cohort pairs table with samplesheet-derived expectations
    (same_patient) and PASS/WARN/FAIL verdicts, and render a full sample x
    sample relatedness heatmap. Report-only; no DAG gating."""
    input:
        pairs="results/qc/somalier/cohort/cohort.pairs.tsv",
        samplesheet=config["samplesheet"],
    output:
        tsv="results/combined/combined_relatedness.tsv",
        heatmap_png="results/combined/relatedness_heatmap.png",
    params:
        ibs0_frac_max=config["params"]["somalier"]["ibs0_frac_max"],
        relatedness_warn=config["params"]["somalier"]["relatedness_warn"],
        unrelated_max=config["params"]["somalier"]["unrelated_max"],
    log:
        "work/logs/somalier_summarise.log",
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/somalier_summarise.R"
