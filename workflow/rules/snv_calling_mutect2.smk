def get_normal_sample(wildcards):
    """Get normal sample name for a given run"""
    return runs_dict[wildcards.run]["normal"]


rule run_mutect2:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        refg=config["paths"]["refs"]["genome_human"],
        regions="refs/regions/regions.bed",
    output:
        vcf=temp("vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf"),
        idx=temp("vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.idx"),
        stats="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.stats",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        normal_name=get_normal_sample,
        germline=config["paths"]["refs"]["germline_resource"],
        ref_path=config["paths"]["refs"]["path"],
    threads: config["resources"]["threads"]
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "logs/{run}/{sample}/Mutect2.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
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
        vcf="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf",
        idx="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.idx",
        stats="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.stats",
        refg=config["paths"]["refs"]["genome_human"],
    output:
        vcf="vcf/{run}/{sample}/mutect/{sample}.mutect2.filtered.vcf",
        idx="vcf/{run}/{sample}/mutect/{sample}.mutect2.filtered.vcf.idx",
        filtering_stats="metrics/{run}/{sample}/{sample}.mutect2.filtered.filteringStats.tsv",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        ref_path=config["paths"]["refs"]["path"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "logs/{run}/{sample}/FilterMutectCalls.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        FilterMutectCalls \
        -R {input.refg} \
        -V {input.vcf} \
        -O {output.vcf} \
        --stats {input.stats} \
        --filtering-stats {output.filtering_stats} \
        > {log} 2>&1
        """


rule sort_mutect2_calls:
    input:
        vcf="vcf/{run}/{sample}/mutect/{sample}.mutect2.filtered.vcf",
        bed="refs/regions/regions.bed.gz",
    output:
        vcf_gz=temp("vcf/{run}/{sample}/mutect/{sample}.mutect2.filtered.vcf.gz"),
        vcf="vcf/{run}/{sample}/mutect/{sample}.mutect2.sorted.vcf",
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        bcftools view -R {input.bed} {output.vcf_gz} | bcftools sort -Ov -o {output.vcf}
        """
