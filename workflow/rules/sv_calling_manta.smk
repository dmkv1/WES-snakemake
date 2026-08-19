rule cap_manta_bam_quality:
    # Manta-only input: caps base qualities before calling, does not touch the
    # BAM used by any other rule. See TODO.md #22.
    #
    # temp(): run_manta is the only consumer, and a run's normal is held until
    # the last tumour of that run has been called. Re-deriving one costs a
    # single-threaded pysam pass over the full BAM (~30-40 min for 10 GB), so
    # pass --notemp when a manta rerun over already-called samples is expected.
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
    output:
        bam=temp("work/manta/{run}/{sample}/{sample}.bqcap.bam"),
        bai=temp("work/manta/{run}/{sample}/{sample}.bqcap.bam.bai"),
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/cap_base_quality.py"


rule run_manta:
    input:
        tumor="work/manta/{run}/{sample}/{sample}.bqcap.bam",
        # configManta.py resolves the indices itself, but they must be declared
        # so temp() keeps them alive until this job runs.
        tumor_bai="work/manta/{run}/{sample}/{sample}.bqcap.bam.bai",
        normal=lambda w: (
            f"work/manta/{w.run}/{runs_dict[w.run]['normal']}/{runs_dict[w.run]['normal']}.bqcap.bam"
            if not is_tumor_only(w) else []
        ),
        normal_bai=lambda w: (
            f"work/manta/{w.run}/{runs_dict[w.run]['normal']}/{runs_dict[w.run]['normal']}.bqcap.bam.bai"
            if not is_tumor_only(w) else []
        ),
        refg=config["refs"]["genome_human"],
        regions_bed=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz",
        regions_tbi=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz.tbi",
    output:
        workspace_dir=directory("work/manta/{run}/{sample}/workspace"),
        workflow="work/manta/{run}/{sample}/runWorkflow.py",
        workflow_pickle="work/manta/{run}/{sample}/runWorkflow.py.config.pickle",
        vcf_sv="work/manta/{run}/{sample}/results/variants/sv_output.vcf.gz",
        tbi_sv="work/manta/{run}/{sample}/results/variants/sv_output.vcf.gz.tbi",
    params:
        res_dir="work/manta/{run}/{sample}",
        normal_arg=lambda w: (
            f"--normalBam work/manta/{w.run}/{runs_dict[w.run]['normal']}/{runs_dict[w.run]['normal']}.bqcap.bam"
            if not is_tumor_only(w) else ""
        ),
        # Tumor-only outputs tumorSV.vcf.gz, paired outputs somaticSV.vcf.gz
        sv_output_name=lambda w: "tumorSV" if is_tumor_only(w) else "somaticSV",
    threads: config["resources"]["threads"]
    resources:
        mem_gb=config["resources"]["manta_max_gb"],
        mem_mb=config["resources"]["manta_mem_mb"],
    conda:
        "../envs/manta.yaml"
    log:
        "work/logs/Manta_{run}_{sample}.log",
    shell:
        """
        configManta.py --exome \
        --referenceFasta {input.refg} \
        {params.normal_arg} \
        --tumorBam {input.tumor} \
        --runDir {params.res_dir} \
        --callRegions {input.regions_bed} \
        > {log} 2>&1

        {output.workflow} \
        --mode local \
        --jobs {threads} --memGb {resources.mem_gb} \
        >> {log} 2>&1

        # Fold paired-BND inversion records into proper INV records (TODO.md #22)
        convertInversion.py $(command -v samtools) {input.refg} \
            {params.res_dir}/results/variants/{params.sv_output_name}.vcf.gz \
            > {params.res_dir}/results/variants/{params.sv_output_name}.converted.vcf \
            2>> {log}

        bgzip -c {params.res_dir}/results/variants/{params.sv_output_name}.converted.vcf > {output.vcf_sv}
        tabix -p vcf {output.vcf_sv}
        """


rule filter_manta_variants:
    input:
        vcf="work/manta/{run}/{sample}/results/variants/sv_output.vcf.gz",
        tbi="work/manta/{run}/{sample}/results/variants/sv_output.vcf.gz.tbi",
    output:
        vcf="work/manta/{run}/{sample}/{sample}.SV.filtered.vcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} | bcftools sort -Ov -o {output.vcf}
        """