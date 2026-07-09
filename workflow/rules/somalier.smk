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


rule somalier_ped:
    """Per-run PED (family=run, sex from samplesheet Chr_sex) so `somalier
    relate` flags sex mismatches against the same Chr_sex value that drives
    CNVkit reference selection."""
    output:
        ped="work/somalier/{run}/{run}.ped",
    params:
        run_samples=lambda w: get_all_samples_for_run(w.run),
    run:
        sex_code = {"XY": "1", "XX": "2"}
        with open(output.ped, "w") as fh:
            for sample in params.run_samples:
                chr_sex = samples.loc[samples["sample"] == sample, "Chr_sex"].iloc[0]
                fh.write(f"{wildcards.run}\t{sample}\t0\t0\t{sex_code.get(chr_sex, '0')}\t0\n")


rule somalier_relate:
    input:
        extracted=lambda w: [
            f"work/somalier/{w.run}/{s}.somalier" for s in get_all_samples_for_run(w.run)
        ],
        ped="work/somalier/{run}/{run}.ped",
    output:
        html="results/qc/somalier/{run}/{run}.html",
        samples="results/qc/somalier/{run}/{run}.samples.tsv",
        pairs="results/qc/somalier/{run}/{run}.pairs.tsv",
    params:
        prefix=lambda w: f"results/qc/somalier/{w.run}/{w.run}",
    log:
        "work/logs/somalier_relate_{run}.log",
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
