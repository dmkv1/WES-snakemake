def get_normal_sample(wildcards):
    """Get normal sample name for a given run"""
    return runs_dict[wildcards.run]["normal"]


rule run_mutect2:
    input:
        tumor=lambda w: f"results/{w.run}/{w.sample}/bam/{w.sample}.bam",
        normal=lambda w: (
            f"results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bam"
            if not is_tumor_only(w) else []
        ),
        refg=config["refs"]["genome_human"],
        regions=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed",
        pon=lambda w: get_pon_path(w) if is_tumor_only(w) else [],
    output:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf",
        idx="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.idx",
        stats="work/mutect2/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.stats",
    params:
        normal_args=lambda w: (
            f"-I results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bam "
            f"-normal {runs_dict[w.run]['normal']}"
            if not is_tumor_only(w) else ""
        ),
        pon_arg=lambda w: f"--panel-of-normals {get_pon_path(w)}" if is_tumor_only(w) else "",
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
        {params.pon_arg} \
        --intervals {input.regions} \
        {params.normal_args} \
        -I {input.tumor} \
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


rule filter_population_af:
    """Filter variants with population AF > threshold (tumor-only only)"""
    input:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.final.vcf",
        gnomad=config["refs"]["germline_resource"],
    output:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.af_filtered.vcf",
    params:
        af_threshold=config.get("tumor_only", {}).get("af_threshold", 0.001),
        is_tumor_only=lambda w: is_tumor_only(w),
    conda:
        "../envs/bcftools.yaml"
    log:
        "work/logs/filter_population_af_{run}_{sample}.log",
    shell:
        """
        if [ "{params.is_tumor_only}" = "True" ]; then
            # Bgzip and index input for annotation
            bgzip -c {input.vcf} > {output.vcf}.tmp.gz 2> {log}
            tabix -p vcf {output.vcf}.tmp.gz 2>> {log}
            # Annotate with gnomAD AF and filter
            bcftools annotate -a {input.gnomad} -c INFO/gnomAD_AF:=INFO/AF \
                {output.vcf}.tmp.gz 2>> {log} | \
            bcftools filter -e 'INFO/gnomAD_AF > {params.af_threshold}' \
                -Ov -o {output.vcf} 2>> {log}
            rm -f {output.vcf}.tmp.gz {output.vcf}.tmp.gz.tbi
        else
            cp {input.vcf} {output.vcf}
        fi
        """


rule funcotator:
    input:
        vcf="work/mutect2/{run}/{sample}/{sample}.mutect2.af_filtered.vcf",
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
