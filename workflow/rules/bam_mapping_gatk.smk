rule fastp_trim:
    input:
        fq1=lambda wildcards: get_fastq1(wildcards),
        fq2=lambda wildcards: get_fastq2(wildcards),
    output:
        fq1=temp("work/fastq/{run}/{sample}/{sample}.{lane}_R1.fq.gz"),
        fq2=temp("work/fastq/{run}/{sample}/{sample}.{lane}_R2.fq.gz"),
        html="results/qc/fastp/{run}/{sample}.{lane}_fastp.html",
        json="results/qc/fastp/{run}/{sample}.{lane}_fastp.json",
    threads: 4
    params:
        detect_adapter=(
            "--detect_adapter_for_pe"
            if config["params"]["fastp"]["detect_adapter_for_pe"]
            else ""
        ),
    conda:
        "../envs/fastp.yaml"
    log:
        "work/logs/fastp_{run}_{sample}.{lane}.log",
    shell:
        "fastp -i {input.fq1} -I {input.fq2} "
        "-o {output.fq1} -O {output.fq2} "
        "-h {output.html} -j {output.json} "
        "{params.detect_adapter} "
        "-w {threads} > {log} 2>&1"


# Per-lane alignment. Each lane is mapped separately with its own read group
# (bwa's inline -R, built from the lane's FASTQ header) and coordinate-sorted,
# so multi-lane samples carry distinct RGs (TODO #15). Single-lane samples are
# just the one-lane case of the same rule. PDX samples map the xengsort graft
# reads. Lanes are gathered by mark_duplicates below.
rule bwa_map:
    input:
        refg=config["refs"]["genome_human"],
        fq1=lambda wildcards: (
            f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.lane}.xengsort-graft.1.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.lane}_R1.fq.gz"
        ),
        fq2=lambda wildcards: (
            f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.lane}.xengsort-graft.2.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.lane}_R2.fq.gz"
        ),
    output:
        temp("results/{run}/{sample}/bam/lanes/{sample}.{lane}.sorted.bam"),
    params:
        rg=get_lane_read_group,
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    log:
        "work/logs/bwamem_{run}_{sample}.{lane}.log",
    shell:
        "(bwa mem -Y -t {threads} -R '{params.rg}' {input.refg} {input.fq1} {input.fq2} "
        "| samtools sort -@ {threads} -T tmp/{wildcards.sample}.{wildcards.lane}.sort -o {output} -) 2> {log}"


# Gather all lanes of a sample and mark duplicates in one pass. MarkDuplicates
# merges the coordinate-sorted per-lane inputs and, because the lanes share one
# library (LB), detects PCR duplicates across lanes. Replaces the former
# AddOrReplaceReadGroups (RG now set inline by bwa) and FixMateInformation
# (coordinate sort now done in bwa_map).
rule mark_duplicates:
    input:
        get_lane_bams,
    output:
        bam=temp("results/{run}/{sample}/bam/{sample}.md.bam"),
        bai=temp("results/{run}/{sample}/bam/{sample}.md.bai"),
        metrics="results/metrics/dupl_metrics_{run}_{sample}.txt",
    params:
        tmp_dir="tmp",
        inputs=lambda wildcards, input: " ".join(f"-I {b}" for b in input),
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
        mem_mb=config["resources"]["mem_mb"],
    log:
        "work/logs/MarkDuplicates_{run}_{sample}.log",
    container:
        config["containers"]["gatk"]
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        MarkDuplicates \
        {params.inputs} \
        -O {output.bam} \
        -M {output.metrics} \
        --CREATE_INDEX true \
        --TMP_DIR tmp \
        > {log} 2>&1
        """


rule create_base_recalibration:
    input:
        bam="results/{run}/{sample}/bam/{sample}.md.bam",
        bai="results/{run}/{sample}/bam/{sample}.md.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda wildcards: config["probe_configs"][
            probe_dict[wildcards.run][wildcards.sample]
        ]["covered_bedfile"],
    output:
        recal_data="results/metrics/{run}_{sample}.recal_data.table",
    params:
        tmp_dir="tmp",
        ref_path=config["refs"]["path"],
        known_sites=lambda _: " ".join(
            f"--known-sites {site}" for site in config["refs"]["known_sites"]
        ),
        interval_padding=config["params"]["bqsr"]["interval_padding"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
        mem_mb=config["resources"]["mem_mb"],
    log:
        "work/logs/BaseRecalibrator_{run}_{sample}.log",
    container:
        config["containers"]["gatk"]
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        BaseRecalibrator \
        -I {input.bam} \
        -O {output.recal_data} \
        -R {input.refg} \
        {params.known_sites} \
        --intervals {input.regions} \
        --interval-padding {params.interval_padding} \
        --tmp-dir {params.tmp_dir} \
        > {log} 2>&1
        """


rule apply_base_recalibration:
    input:
        refg=config["refs"]["genome_human"],
        bam="results/{run}/{sample}/bam/{sample}.md.bam",
        bsqr_recal="results/metrics/{run}_{sample}.recal_data.table",
    output:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
    params:
        tmp_dir="tmp",
        ref_path=config["refs"]["path"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
        mem_mb=config["resources"]["mem_mb"],
    log:
        "work/logs/ApplyBQSR_{run}_{sample}.log",
    container:
        config["containers"]["gatk"]
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        ApplyBQSR \
        -R {input.refg} \
        -I {input.bam} \
        --bqsr-recal-file {input.bsqr_recal} \
        -O {output.bam} \
        --tmp-dir {params.tmp_dir} \
        > {log} 2>&1
        """


rule mosdepth:
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        regions_bed=lambda wildcards: config["probe_configs"][
            probe_dict[wildcards.run][wildcards.sample]
        ]["covered_bedfile"],
    output:
        summary="results/metrics/{run}/{sample}.mosdepth.summary.txt",
        region_dist="results/metrics/{run}/{sample}.mosdepth.region.dist.txt",
        thresholds="results/metrics/{run}/{sample}.thresholds.bed.gz",
    params:
        prefix="results/metrics/{run}/{sample}",
    conda:
        "../envs/qc.yaml"
    log:
        "work/logs/mosdepth_{run}_{sample}.log",
    shell:
        "mosdepth --by {input.regions_bed} --thresholds 10,20,30,50 {params.prefix} {input.bam} > {log} 2>&1"
