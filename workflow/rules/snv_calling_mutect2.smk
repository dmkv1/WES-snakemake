# def get_tumor_inputs(wildcards):
#     """Generate -I flags for tumor samples"""
#     tumors = runs_dict[wildcards.run]["tumors"]
#     return [f"-I bam/{wildcards.run}/{sample}.bam" for sample in tumors]


def get_normal_sample(wildcards):
    """Get normal sample name for a given run"""
    return runs_dict[wildcards.run]["normal"]


rule run_mutect2:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        # Multisample mode - incompartible with the rest of the callers
        #tumors=lambda w: expand(
        #    "bam/{run}/{sample}.bam", run=w.run, sample=runs_dict[w.run]["tumors"]
        #),
        refg=config["paths"]["refs"]["genome_human"],
        regions="refs/regions/regions.bed",
    output:
        vcf=temp("vcf/{run}/{sample}/{sample}.mutect2.unfiltered.vcf"),
        idx=temp("vcf/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.idx"),
        stats="vcf/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.stats",
    params:
        gatk_image=config["tools"]["gatk_image"],
        gatk_ver=config["tools"]["gatk_version"],
        normal_name=get_normal_sample,
        # Multisample mode - incompartible with the rest of the callers
        # tumor_inputs=get_tumor_inputs,
        germline=config["paths"]["refs"]["germline_resource"],
        ref_path=config["paths"]["refs"]["path"],
    threads: config["resources"]["threads"]
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
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
        -O {output.vcf}
        """


rule filter_mutect2_calls:
    input:
        vcf="vcf/{run}/{sample}/{sample}.mutect2.unfiltered.vcf",
        stats="vcf/{run}/{sample}/{sample}.mutect2.unfiltered.vcf.stats",
        refg=config["paths"]["refs"]["genome_human"],
    output:
        vcf="vcf/{run}/{sample}/{sample}.mutect2.vcf",
        idx="vcf/{run}/{sample}/{sample}.mutect2.vcf.idx",
        filtering_stats="metrics/{run}/{sample}/{sample}.mutect2.filteringStats.tsv",
    params:
        gatk_image=config["tools"]["gatk_image"],
        gatk_ver=config["tools"]["gatk_version"],
        ref_path=config["paths"]["refs"]["path"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
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
        --filtering-stats {output.filtering_stats}
        """
