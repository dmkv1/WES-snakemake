# QC metric sources (mirrors WES-PON-smk/workflow/rules/qc.smk):
#   - FastQC reads the trimmed FASTQs (not the BQSR'd BAM, whose recalibrated
#     base qualities make FastQC's chemistry plots meaningless).
#   - CollectHsMetrics reports capture efficiency / on-target rate / fold
#     enrichment from the final recalibrated BAM, against the kit's bait+target
#     intervals (built once per probe kit by bed_to_interval_list).
#   - MultiQC depends on the raw metric files it actually parses (fastp .json,
#     fastqc .zip, Picard metrics, mosdepth summary, hs_metrics), not rendered HTML.


def get_all_samples_for_run(run):
    """All samples (normal + tumors) for a run, excluding None."""
    normal = runs_dict[run]["normal"]
    tumors = runs_dict[run]["tumors"]
    return ([normal] if normal else []) + tumors


# Picard HsMetrics distinguishes baits (where probes hybridize -> the "Covered"
# BED) from targets (regions of interest -> the "Regions" BED).
_BED_FOR_KIND = {"bait": "covered_bedfile", "target": "target_regions_bedfile"}


rule bed_to_interval_list:
    input:
        bed=lambda wildcards: config["probe_configs"][wildcards.probes][
            _BED_FOR_KIND[wildcards.kind]
        ],
        refg=config["refs"]["genome_human"],
    output:
        "work/intervals/{probes}.{kind}.interval_list",
    log:
        "work/logs/BedToIntervalList_{probes}_{kind}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    wildcard_constraints:
        kind="bait|target",
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            BedToIntervalList \
            -I {input.bed} -O {output} \
            -SD {input.refg} \
            -TMP_DIR {params.tmp_dir} \
            > {log} 2>&1
        """


rule fastqc:
    input:
        fq1="work/fastq/{run}/{sample}/{sample}_R1.fq.gz",
        fq2="work/fastq/{run}/{sample}/{sample}_R2.fq.gz",
    output:
        zip1="results/qc/fastqc/{run}/{sample}_R1_fastqc.zip",
        zip2="results/qc/fastqc/{run}/{sample}_R2_fastqc.zip",
    log:
        "work/logs/fastqc_{run}_{sample}.log",
    conda:
        "../envs/qc.yaml"
    threads: 2
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.zip1),
    shell:
        "fastqc {input.fq1} {input.fq2} -o {params.out_dir} -t {threads} > {log} 2>&1"


rule collect_hs_metrics:
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        refg=config["refs"]["genome_human"],
        baits=lambda wildcards: f"work/intervals/{probe_dict[wildcards.run][wildcards.sample]}.bait.interval_list",
        targets=lambda wildcards: f"work/intervals/{probe_dict[wildcards.run][wildcards.sample]}.target.interval_list",
    output:
        metrics="results/metrics/{run}_{sample}.hs_metrics.txt",
    log:
        "work/logs/CollectHsMetrics_{run}_{sample}.log",
    container:
        config["containers"]["gatk"]
    resources:
        java_min_gb=config["resources"]["java_min_gb"],
        java_max_gb=config["resources"]["java_max_gb"],
    params:
        tmp_dir="tmp",
    shell:
        """
        gatk --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
            CollectHsMetrics \
            -I {input.bam} -O {output.metrics} \
            -R {input.refg} \
            --BAIT_INTERVALS {input.baits} \
            --TARGET_INTERVALS {input.targets} \
            -TMP_DIR {params.tmp_dir} \
            > {log} 2>&1
        """


rule multiqc:
    input:
        fastp_json=[
            f"results/qc/fastp/{run}/{sample}_fastp.json"
            for run in runs_dict
            for sample in get_all_samples_for_run(run)
        ],
        fastqc_zip=[
            f"results/qc/fastqc/{run}/{sample}_R{read}_fastqc.zip"
            for run in runs_dict
            for sample in get_all_samples_for_run(run)
            for read in (1, 2)
        ],
        dupl_metrics=[
            f"results/metrics/dupl_metrics_{run}_{sample}.txt"
            for run in runs_dict
            for sample in get_all_samples_for_run(run)
        ],
        mosdepth=[
            f"results/metrics/{run}_{sample}.mosdepth.summary.txt"
            for run in runs_dict
            for sample in get_all_samples_for_run(run)
        ],
        hs_metrics=[
            f"results/metrics/{run}_{sample}.hs_metrics.txt"
            for run in runs_dict
            for sample in get_all_samples_for_run(run)
        ],
        config="multiqc_config.yaml",
    output:
        "results/qc/multiqc_report.html",
    log:
        "work/logs/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc -c {input.config} results/ work/logs/ -o results/qc/ --force > {log} 2>&1"
