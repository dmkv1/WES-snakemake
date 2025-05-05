rule run_manta:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        refg=config["paths"]["refs"]["genome_human"],
        regions_bed="refs/regions/regions.bed.gz",
        regions_tbi="refs/regions/regions.bed.gz.tbi",
    output:
        workflow="vcf/{run}/{sample}/manta/runWorkflow.py",
        workflow_pickle="vcf/{run}/{sample}/manta/runWorkflow.py.config.pickle",
        results_dir=directory("vcf/{run}/{sample}/manta/results"),
        workspace_dir=directory("vcf/{run}/{sample}/manta/workspace"),
        vcf_SmallIndels="vcf/{run}/{sample}/manta/results/variants/candidateSmallIndels.vcf.gz",
        tbi_SmallIndels="vcf/{run}/{sample}/manta/results/variants/candidateSmallIndels.vcf.gz.tbi",
        vcf_candidateSV="vcf/{run}/{sample}/manta/results/variants/candidateSV.vcf.gz",
        tbi_candidateSV="vcf/{run}/{sample}/manta/results/variants/candidateSV.vcf.gz.tbi",
        vcf_diploidSV="vcf/{run}/{sample}/manta/results/variants/diploidSV.vcf.gz",
        tbi_diploidSV="vcf/{run}/{sample}/manta/results/variants/diploidSV.vcf.gz.tbi",
        vcf_somaticSV="vcf/{run}/{sample}/manta/results/variants/somaticSV.vcf.gz",
        tbi_somaticSV="vcf/{run}/{sample}/manta/results/variants/somaticSV.vcf.gz.tbi",
    params:
        res_dir="vcf/{run}/{sample}/manta",
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["strelka_max_gb"],
    conda:
        "../envs/manta.yaml"
    log:
        "logs/{run}/{sample}/Manta.log",
    shell:
        """
        configManta.py --exome \
        --referenceFasta {input.refg} \
        --normalBam {input.normal} \
        --tumorBam {input.tumor} \
        --runDir {params.res_dir} \
        --callRegions {input.regions_bed} \
        > {log} 2>&1

        {output.workflow} \
        --mode local \
        --jobs {resources.threads} --memGb {resources.mem_gb} \
        >> {log} 2>&1
        """


rule run_strelka:
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
        workflow_pickle="vcf/{run}/{sample}/strelka/runWorkflow.py.config.pickle",
        results_dir=directory("vcf/{run}/{sample}/strelka/results"),
        workspace_dir=directory("vcf/{run}/{sample}/strelka/workspace"),
        vcf_snv="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz",
        tbi_snv="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz.tbi",
        vcf_indel="vcf/{run}/{sample}/strelka/results/variants/somatic.indels.vcf.gz",
        tbi_indel="vcf/{run}/{sample}/strelka/results/variants/somatic.indels.vcf.gz.tbi",
    params:
        run_path="vcf/{run}/{sample}/strelka",
        ref_path=config["paths"]["refs"]["path"],
        strelka_image=config["tools"]["strelka"]["image"],
        strelka_ver=config["tools"]["strelka"]["version"],
    resources:
        threads=config["resources"]["threads"],
        mem_gb=config["resources"]["strelka_max_gb"],
    log:
        "logs/{run}/{sample}/Strelka.log",
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
        --indelCandidates {input.manta_vcf} \
        > {log} 2>&1

        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.strelka_image}:{params.strelka_ver} \
        {output.workflow} --mode local --jobs 8 --memGb {resources.mem_gb} \
        >> {log} 2>&1
        """


rule sort_strelka_calls:
    input:
        vcf_gz="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz",
        tbi="vcf/{run}/{sample}/strelka/results/variants/somatic.snvs.vcf.gz.tbi",
        bed="refs/regions/regions.bed.gz",
    output:
        vcf="vcf/{run}/{sample}/strelka/vcf_sorted/{sample}.strelka.sorted.vcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -R {input.bed} {input.vcf_gz} | bcftools sort -Ov -o {output.vcf}
        """
