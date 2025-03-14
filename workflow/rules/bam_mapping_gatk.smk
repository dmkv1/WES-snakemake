rule bwa_map:
    input:
        refg=config["paths"]["refs"]["genome_human"],
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
    threads: config["resources"]["threads"]
    log:
        "logs/{run}/{sample}/bwamem.log",
    shell:
        "(bwa mem -M -t {threads} {input.refg} {input.fq1} {input.fq2} | samtools view -Sb - > {output}) 2> {log}"


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
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    log:
        "logs/{run}/{sample}/AddOrReplaceReadGroups.log",
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
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
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    log:
        "logs/{run}/{sample}/FixMateInformation.log",
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
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
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "logs/{run}/{sample}/MarkDuplicates.log",
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
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
        refg=config["paths"]["refs"]["genome_human"],
    output:
        recal_data="metrics/{run}/{sample}/{sample}.recal_data.table",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        tmp_dir="tmp",
        ref_path=config["paths"]["refs"]["path"],
        known_sites=lambda _: " ".join(
            f"--known-sites {site}" for site in config["paths"]["refs"]["known_sites"]
        ),
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "logs/{run}/{sample}/BaseRecalibrator.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        BaseRecalibrator \
        -I {input.bam} \
        -O {output.recal_data} \
        -R {input.refg} \
        {params.known_sites} \
        --tmp-dir {params.tmp_dir} \
        > {log} 2>&1
        """


rule apply_base_recalibration:
    input:
        refg=config["paths"]["refs"]["genome_human"],
        bam="bam/{run}/{sample}.md.bam",
        bsqr_recal="metrics/{run}/{sample}/{sample}.recal_data.table",
    output:
        bam="bam/{run}/{sample}.bam",
        bai="bam/{run}/{sample}.bai",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        tmp_dir="tmp",
        ref_path=config["paths"]["refs"]["path"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "logs/{run}/{sample}/ApplyBQSR.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        ApplyBQSR \
        -R {input.refg} \
        -I {input.bam} \
        --bqsr-recal-file {input.bsqr_recal} \
        -O {output.bam} \
        --tmp-dir {params.tmp_dir} \
        > {log} 2>&1
        """
