rule configure_manta:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        refg=config["paths"]["refs"]["genome_human"],
        regions_bed="refs/regions/regions.bed.gz",
        regions_tbi="refs/regions/regions.bed.gz.tbi",
    output:
        workflow="vcf/{run}/{sample}/manta/runWorkflow.py",
    params:
        res_dir="vcf/{run}/{sample}/manta",
    conda:
        "../envs/strelka.yaml"
    shell:
        """
        configManta.py --exome \
        --referenceFasta {input.refg} \
        --normalBam {input.normal} \
        --tumorBam {input.tumor} \
        --runDir {params.res_dir} \
        --callRegions {input.regions_bed}
        """


rule run_manta:
    input:
        workflow="vcf/{run}/{sample}/manta/runWorkflow.py",
    output:
        vcf="vcf/{run}/{sample}/manta/results/variants/candidateSmallIndels.vcf.gz",
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["strelka_max_gb"],
    conda:
        "../envs/strelka.yaml"
    shell:
        "{input.workflow} --mode local --jobs {resources.threads} --memGb {resources.mem_gb}"


rule configure_strelka:
    input:
        normal_bam=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        normal_bai=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bai",
        tumor_bam=lambda w: f"bam/{w.run}/{w.sample}.bam",
        tumor_bai=lambda w: f"bam/{w.run}/{w.sample}.bai",
        refg=config["paths"]["refs"]["genome_human"],
        manta_vcf="vcf/{run}/{sample}/manta/results/variants/candidateSmallIndels.vcf.gz",
        regions_bed="refs/regions/regions.bed.gz",
        regions_tbi="refs/regions/regions.bed.gz.tbi",
    output:
        workflow="vcf/{run}/{sample}/strelka/runWorkflow.py",
    params:
        res_dir="vcf/{run}/{sample}/strelka",
    conda:
        "../envs/strelka.yaml"
    shell:
        """
        configureStrelkaSomaticWorkflow.py --exome \
        --runDir {params.res_dir} \
        --referenceFasta {input.refg} \
        --normalBam {input.normal_bam} \
        --tumorBam {input.tumor_bam} \
        --callRegions {input.regions_bed} \
        --indelCandidates {input.manta_vcf}
        """


rule run_strelka:
    input:
        workflow="vcf/{run}/{sample}/strelka/runWorkflow.py",
    output:
        vcf="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz",
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["strelka_max_gb"],
    conda:
        "../envs/strelka.yaml"
    shell:
        "{input.workflow} --mode local --jobs 8 --memGb {resources.mem_gb}"
