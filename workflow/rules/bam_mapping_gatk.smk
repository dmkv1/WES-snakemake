rule fastp_trim:
    input:
        fq1=lambda wildcards: get_fastq1(wildcards),
        fq2=lambda wildcards: get_fastq2(wildcards),
    output:
        fq1=temp("work/fastq/{run}/{sample}/{sample}.{unit}_R1.fq.gz"),
        fq2=temp("work/fastq/{run}/{sample}/{sample}.{unit}_R2.fq.gz"),
        html="results/qc/fastp/{run}/{sample}.{unit}_fastp.html",
        json="results/qc/fastp/{run}/{sample}.{unit}_fastp.json",
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
        "work/logs/fastp_{run}_{sample}.{unit}.log",
    shell:
        "fastp -i {input.fq1} -I {input.fq2} "
        "-o {output.fq1} -O {output.fq2} "
        "-h {output.html} -j {output.json} "
        "{params.detect_adapter} "
        "-w {threads} > {log} 2>&1"


# Per-unit alignment. Each unit is mapped with its own read group, resolved and
# validated at load time (workflow/scripts/units.py), so multi-lane samples
# carry distinct RGs and BQSR can model each unit's error profile separately
# (TODO #15). A single-FASTQ sample is just the one-unit case. PDX samples map
# the xengsort graft reads. Units are gathered by mark_duplicates below.
#
# The output is deliberately NOT coordinate-sorted: bwa emits query-grouped
# reads, samtools fixmate requires that grouping, and MarkDuplicates wants it
# too (see below). The single coordinate sort happens once, after duplicates
# are marked.
#
# -K fixes bwa's input chunk size so alignment is independent of thread count,
# which is what makes two runs of the same FASTQ comparable.
rule bwa_map:
    input:
        refg=config["refs"]["genome_human"],
        fq1=lambda wildcards: (
            f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.unit}.xengsort-graft.1.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.unit}_R1.fq.gz"
        ),
        fq2=lambda wildcards: (
            f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.unit}.xengsort-graft.2.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.{wildcards.unit}_R2.fq.gz"
        ),
    output:
        temp("results/{run}/{sample}/bam/units/{sample}.{unit}.qgrp.bam"),
    params:
        rg=get_read_group,
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    log:
        "work/logs/bwamem_{run}_{sample}.{unit}.log",
    shell:
        "(bwa mem -Y -K 100000000 -t {threads} -R '{params.rg}' "
        "{input.refg} {input.fq1} {input.fq2} "
        "| samtools fixmate -m -O bam,level=1 - {output}) 2> {log}"


# Gather all units of a sample and mark duplicates in one pass. Units sharing a
# library (LB) have their PCR duplicates detected across units; units from
# different libraries are treated separately, which is what preserves the
# independent evidence of two preps of one sample.
#
# ASSUME_SORT_ORDER queryname consumes bwa's query-grouped output directly, so
# secondary and supplementary reads are marked correctly. This is the ordering
# GATK Best Practices uses (MarkDuplicates before the coordinate sort), and it
# is why sort_bam below exists as its own rule.
rule mark_duplicates:
    input:
        get_unit_bams,
    output:
        bam=temp("results/{run}/{sample}/bam/{sample}.md.qgrp.bam"),
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
        --ASSUME_SORT_ORDER queryname \
        --TMP_DIR tmp \
        > {log} 2>&1
        """


# samtools sort -m is a per-thread buffer, so the real footprint is
# sort_mem x threads. Without this the scheduler sees the 4000 MB default and
# can start a full-size GATK job beside the sort, which overcommits the host.
def _sort_mem_mb():
    raw = str(config["resources"]["sort_mem"]).strip()
    unit = raw[-1].upper()
    if unit in ("K", "M", "G"):
        value = float(raw[:-1])
        mb = {"K": value / 1024, "M": value, "G": value * 1024}[unit]
    else:
        mb = float(raw) / (1024 * 1024)  # a bare number is bytes, as in samtools
    # 15% above the buffers themselves, for the merge pass and the BGZF output.
    return int(mb * config["resources"]["threads"] * 1.15)


# The one coordinate sort. --write-index with the ##idx## form names the index
# {sample}.md.bai rather than {sample}.md.bam.bai, which is the filename the
# BQSR rules below already expect.
rule sort_bam:
    input:
        "results/{run}/{sample}/bam/{sample}.md.qgrp.bam",
    output:
        bam=temp("results/{run}/{sample}/bam/{sample}.md.bam"),
        bai=temp("results/{run}/{sample}/bam/{sample}.md.bai"),
    params:
        sort_mem=config["resources"]["sort_mem"],
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    resources:
        mem_mb=_sort_mem_mb(),
    log:
        "work/logs/sort_bam_{run}_{sample}.log",
    shell:
        "samtools sort -@ {threads} -m {params.sort_mem} --write-index "
        "-T tmp/{wildcards.run}_{wildcards.sample}.sort "
        "-o {output.bam}##idx##{output.bai} {input} 2> {log}"


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
        "mosdepth --no-per-base --by {input.regions_bed} --thresholds 10,20,30,50 {params.prefix} {input.bam} > {log} 2>&1"
