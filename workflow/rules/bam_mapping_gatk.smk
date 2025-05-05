rule bwa_map:
    input:
        fq1=lambda wildcards: (
            f"fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.xengsort-graft.1.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else get_fastq1(wildcards)
        ),
        fq2=lambda wildcards: (
            f"fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.xengsort-graft.2.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else get_fastq2(wildcards)
        ),
    output:
        temp("bam/{run}/{sample}.raw.bam"),
    params:
        refg=get_ref_path(config["refs"]["genome_human"], use_container=False),
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    log:
        "logs/{run}/{sample}/bwamem.log",
    shell:
        "(bwa mem -M -t {threads} {params.refg} {input.fq1} {input.fq2} | samtools view -Sb - > {output}) 2> {log}"


rule add_read_groups:
    input:
        "bam/{run}/{sample}.raw.bam",
    output:
        temp("bam/{run}/{sample}.rg.bam"),
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        tmp_dir="tmp",
        rg=get_read_group_params,
    resources:
        memory_min_gb=config["resources"]["memory_min_gb"],
        memory_max_gb=config["resources"]["memory_max_gb"],
    log:
        "logs/{run}/{sample}/AddOrReplaceReadGroups.log",
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        AddOrReplaceReadGroups \
        -I {input} \
        -O {output} \
        -ID {params.rg[RGID]} \
        -SM {params.rg[RGSM]} \
        -PU {params.rg[RGPU]} \
        -LB {params.rg[RGLB]} \
        -PL {params.rg[RGPL]} \
        -TMP_DIR {params.tmp_dir} \
        > {log} 2>&1
        """


rule fix_mate_info:
    input:
        "bam/{run}/{sample}.rg.bam",
    output:
        temp("bam/{run}/{sample}.fixmate.bam"),
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        tmp_dir="tmp",
    resources:
        memory_min_gb=config["resources"]["memory_min_gb"],
        memory_max_gb=config["resources"]["memory_max_gb"],
    log:
        "logs/{run}/{sample}/FixMateInformation.log",
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        FixMateInformation \
        -I {input} \
        -O {output} \
        -SO coordinate \
        -VALIDATION_STRINGENCY SILENT \
        -TMP_DIR {params.tmp_dir} \
        > {log} 2>&1
        """


rule mark_duplicates:
    input:
        "bam/{run}/{sample}.fixmate.bam",
    output:
        bam=temp("bam/{run}/{sample}.md.bam"),
        metrics="metrics/{run}/{sample}/{sample}.dupl_metrics.txt",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        tmp_dir="tmp",
    resources:
        memory_max_gb=config["resources"]["memory_max_gb"],
        memory_min_gb=config["resources"]["memory_min_gb"],
    log:
        "logs/{run}/{sample}/MarkDuplicates.log",
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        MarkDuplicates \
        -I {input} \
        -O {output.bam} \
        -M {output.metrics} \
        --CREATE_INDEX false \
        --TMP_DIR tmp \
        > {log} 2>&1
        """


rule create_base_recalibration:
    input:
        bam="bam/{run}/{sample}.md.bam",
    output:
        recal_data="metrics/{run}/{sample}/{sample}.recal_data.table",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        tmp_dir="tmp",
        refg=get_ref_path(config["refs"]["genome_human"], use_container=True),
        known_sites=lambda _: " ".join(
            f"--known-sites {get_ref_path(site, use_container= True)}"
            for site in config["refs"]["known_sites"]
        ),
    resources:
        memory_max_gb=config["resources"]["memory_max_gb"],
        memory_min_gb=config["resources"]["memory_min_gb"],
    log:
        "logs/{run}/{sample}/BaseRecalibrator.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        BaseRecalibrator \
        -I {input.bam} \
        -O {output.recal_data} \
        -R {params.refg} \
        {params.known_sites} \
        --tmp-dir {params.tmp_dir} \
        > {log} 2>&1
        """


rule apply_base_recalibration:
    input:
        bam="bam/{run}/{sample}.md.bam",
        bsqr_recal="metrics/{run}/{sample}/{sample}.recal_data.table",
    output:
        bam="bam/{run}/{sample}.bam",
        bai="bam/{run}/{sample}.bai",
    params:
        refg=get_ref_path(config["refs"]["genome_human"], use_container=True),
        tmp_dir="tmp",
    resources:
        memory_max_gb=config["resources"]["memory_max_gb"],
        memory_min_gb=config["resources"]["memory_min_gb"],
    log:
        "logs/{run}/{sample}/ApplyBQSR.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        ApplyBQSR \
        -R {params.refg} \
        -I {input.bam} \
        --bqsr-recal-file {input.bsqr_recal} \
        -O {output.bam} \
        --tmp-dir {params.tmp_dir} \
        > {log} 2>&1
        """
