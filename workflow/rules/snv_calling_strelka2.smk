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
        "../envs/manta.yaml"
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
        vcf_SmallIndels="vcf/{run}/{sample}/manta/results/variants/candidateSmallIndels.vcf.gz",
        tbi_SmallIndels="vcf/{run}/{sample}/manta/results/variants/candidateSmallIndels.vcf.gz.tbi",
        vcf_candidateSV="vcf/{run}/{sample}/manta/results/variants/candidateSV.vcf.gz",
        tbi_candidateSV="vcf/{run}/{sample}/manta/results/variants/candidateSV.vcf.gz.tbi",
        vcf_diploidSV="vcf/{run}/{sample}/manta/results/variants/diploidSV.vcf.gz",
        tbi_diploidSV="vcf/{run}/{sample}/manta/results/variants/diploidSV.vcf.gz.tbi",
        vcf_somaticSV="vcf/{run}/{sample}/manta/results/variants/somaticSV.vcf.gz",
        tbi_somaticSV="vcf/{run}/{sample}/manta/results/variants/somaticSV.vcf.gz.tbi",
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["strelka_max_gb"],
    conda:
        "../envs/manta.yaml"
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
        run_path="vcf/{run}/{sample}/strelka",
        ref_path=config["paths"]["refs"]["path"],
        strelka_image=config["tools"]["strelka_image"],
        strelka_ver=config["tools"]["strelka_version"],
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.strelka_image}:{params.strelka_ver} \
        configureStrelkaSomaticWorkflow.py --exome \
        --runDir {params.run_path} \
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
        vcf_snv="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz",
        tbi_snv="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz.tbi",
        vcf_indel="vcf/{run}/{sample}/strelka/results/variants/somatic.indels.vcf.gz",
        tbi_indel="vcf/{run}/{sample}/strelka/results/variants/somatic.indels.vcf.gz.tbi",
    params:
        run_path="vcf/{run}/{sample}/strelka",
        ref_path=config["paths"]["refs"]["path"],
        strelka_image=config["tools"]["strelka_image"],
        strelka_ver=config["tools"]["strelka_version"],
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["strelka_max_gb"],
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.strelka_image}:{params.strelka_ver} \
        {input.workflow} --mode local --jobs 8 --memGb {resources.mem_gb}
        """


rule concat_strelka_vcfs:
    input:
        vcf_snv="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz",
        tbi_snv="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz.tbi",
        vcf_indel="vcf/{run}/{sample}/strelka/results/variants/somatic.indels.vcf.gz",
        tbi_indel="vcf/{run}/{sample}/strelka/results/variants/somatic.indels.vcf.gz.tbi",
    output:
        vcf="vcf/{run}/{sample}/{sample}.strelka.vcf",
    conda:
        "../envs/manta.yaml"
    shell:
        """
        bcftools concat -a --allow-overlaps {input.vcf_snv} {input.vcf_indel} | \
        bcftools sort -Ov -o {output.vcf}
        """
