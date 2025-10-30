rule fastp_trim:
    input:
        fq1=lambda wildcards: get_fastq1(wildcards),
        fq2=lambda wildcards: get_fastq2(wildcards),
    output:
        fq1=temp("work/fastq/{run}/{sample}/{sample}.trimmed.1.fq.gz"),
        fq2=temp("work/fastq/{run}/{sample}/{sample}.trimmed.2.fq.gz"),
        html="results/qc/fastp/{run}/{sample}_fastp.html",
        json="results/qc/fastp/{run}/{sample}_fastp.json",
    threads: 4
    conda:
        "../envs/fastp.yaml"
    log:
        "work/logs/fastp_{run}_{sample}.log",
    shell:
        "fastp -i {input.fq1} -I {input.fq2} "
        "-o {output.fq1} -O {output.fq2} "
        "-h {output.html} -j {output.json} "
        "-w {threads} > {log} 2>&1"


rule bwa_map:
    input:
        refg=config["refs"]["genome_human"],
        fq1=lambda wildcards: (
            f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.xengsort-graft.1.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.trimmed.1.fq.gz"
        ),
        fq2=lambda wildcards: (
            f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.xengsort-graft.2.fq.gz"
            if wildcards.run in pdx_dict
            and wildcards.sample in pdx_dict[wildcards.run]
            else f"work/fastq/{wildcards.run}/{wildcards.sample}/{wildcards.sample}.trimmed.2.fq.gz"
        ),
    output:
        temp("results/{run}/{sample}/bam/{sample}.raw.bam"),
    conda:
        "../envs/bwamem.yaml"
    threads: config["resources"]["threads"]
    log:
        "work/logs/bwamem_{run}_{sample}.log",
    shell:
        "(bwa mem -M -t {threads} {input.refg} {input.fq1} {input.fq2} | samtools view -Sb - > {output}) 2> {log}"


rule add_read_groups:
    input:
        "results/{run}/{sample}/bam/{sample}.raw.bam",
    output:
        temp("results/{run}/{sample}/bam/{sample}.rg.bam"),
    params:
        tmp_dir="tmp",
        rg=get_read_group_params,
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    log:
        "work/logs/AddOrReplaceReadGroups_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
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
        "results/{run}/{sample}/bam/{sample}.rg.bam",
    output:
        temp("results/{run}/{sample}/bam/{sample}.fixmate.bam"),
    params:
        tmp_dir="tmp",
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    log:
        "work/logs/FixMateInformation_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
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
        "results/{run}/{sample}/bam/{sample}.fixmate.bam",
    output:
        bam=temp("results/{run}/{sample}/bam/{sample}.md.bam"),
        metrics="results/metrics/dupl_metrics_{run}_{sample}.txt",
    params:
        tmp_dir="tmp",
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "work/logs/MarkDuplicates_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
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
        bam="results/{run}/{sample}/bam/{sample}.md.bam",
        refg=config["refs"]["genome_human"],
    output:
        recal_data="results/metrics/{run}_{sample}.recal_data.table",
    params:
        tmp_dir="tmp",
        ref_path=config["refs"]["path"],
        known_sites=lambda _: " ".join(
            f"--known-sites {site}" for site in config["refs"]["known_sites"]
        ),
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "work/logs/BaseRecalibrator_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
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
    log:
        "work/logs/ApplyBQSR_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
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
        ]["regions_bedfile"],
    output:
        summary="results/metrics/{run}_{sample}.mosdepth.summary.txt",
        thresholds="results/metrics/{run}_{sample}.thresholds.bed.gz",
    params:
        prefix="results/metrics/{run}_{sample}",
    conda:
        "../envs/qc.yaml"
    log:
        "work/logs/mosdepth_{run}_{sample}.log",
    shell:
        "mosdepth --by {input.regions_bed} --thresholds 10,20,30,50 {params.prefix} {input.bam} > {log} 2>&1"


rule fastqc:
    input:
        "results/{run}/{sample}/bam/{sample}.bam",
    output:
        html="results/qc/fastqc/{run}/{sample}_fastqc.html",
        zip="results/qc/fastqc/{run}/{sample}_fastqc.zip",
    conda:
        "../envs/qc.yaml"
    threads: 2
    shell:
        "fastqc {input} -o results/qc/fastqc/{run} -t {threads}"


rule multiqc:
    input:
        [
            f"results/qc/fastqc/{run}/{sample}_fastqc.zip"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
        [
            f"results/qc/fastp/{run}/{sample}_fastp.zip"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
        [
            f"results/metrics/{run}_{sample}.mosdepth.summary.txt"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
        [
            f"results/metrics/{run}_{sample}.mosdepth.region.dist.txt"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
        [
            f"results/metrics/{run}_{sample}.recal_data.table"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
        [
            f"results/metrics/dupl_metrics_{run}_{sample}.txt"
            for run in runs_dict
            for sample in ([runs_dict[run]["normal"]] + runs_dict[run]["tumors"])
        ],
    output:
        "results/qc/multiqc_report.html",
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc results/ work/logs/ -o results/qc/ --force"
