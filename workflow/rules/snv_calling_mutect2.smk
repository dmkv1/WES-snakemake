def get_normal_sample(wildcards):
    """Get normal sample name for a given run"""
    return runs_dict[wildcards.run]["normal"]


rule run_mutect2:
    input:
        normal=lambda w: f"results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"results/{w.run}/{w.sample}/bam/{w.sample}.bam",
        refg=config["refs"]["genome_human"],
        regions=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed",
    output:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf",
        idx="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.idx",
        stats="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.stats",
    params:
        normal_name=get_normal_sample,
        germline=config["refs"]["germline_resource"],
        ref_path=config["refs"]["path"],
    threads: config["resources"]["threads"]
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "work/logs/Mutect2_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        Mutect2 \
        --native-pair-hmm-threads {threads} \
        -R {input.refg} \
        --germline-resource {params.germline} \
        --intervals {input.regions} \
        -I {input.normal} \
        -I {input.tumor} \
        -normal {params.normal_name} \
        -O {output.vcf} \
        > {log} 2>&1
        """


rule filter_mutect2_calls:
    input:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf",
        idx="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.idx",
        stats="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.stats",
        refg=config["refs"]["genome_human"],
    output:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.filtered.vcf",
        idx="work/mutect2/{run}/{sample}/{sample}.mutect2.filtered.vcf.idx",
        filtering_stats="results/metrics/mutect2_filteringStats_{run}_{sample}.tsv",
    params:
        ref_path=config["refs"]["path"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "work/logs/FilterMutectCalls_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        FilterMutectCalls \
        -R {input.refg} \
        -V {input.vcf} \
        -O {output.vcf} \
        --stats {input.stats} \
        --filtering-stats {output.filtering_stats} \
        > {log} 2>&1
        """


rule filter_and_sort_mutect2_calls:
    input:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.filtered.vcf",
        regions=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz",
    output:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.final.vcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} | bcftools sort -Ov -o {output.vcf}
        """


rule funcotator:
    input:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.final.vcf",
        refg=config["refs"]["genome_human"],
    output:
        vcf="results/{run}/{sample}/{sample}.SNV.vcf",
        idx="results/{run}/{sample}/{sample}.SNV.vcf.idx",
    params:
        genome_ver=config["refs"]["funcotator_data_sources"]["genome_version"],
        data_sources=config["refs"]["funcotator_data_sources"]["path"],
        ref_path=config["refs"]["path"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "work/logs/Funcotator_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        Funcotator \
        --reference {input.refg} \
        --ref-version {params.genome_ver} \
        --data-sources-path {params.data_sources} \
        --output-file-format VCF \
        --variant {input.vcf} \
        --output {output.vcf} \
        > {log} 2>&1
        """
