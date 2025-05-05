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
        refg=config["refs"]["genome_human"],
        regions="refs/regions/regions.bed",
    output:
        vcf=temp("vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf"),
        idx=temp("vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.idx"),
        stats="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.stats",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        normal_name=get_normal_sample,
        # Multisample mode - incompartible with the rest of the callers
        # tumor_inputs=get_tumor_inputs,
        germline=config["refs"]["germline_resource"],
        refdir=config["refs"]["path"],
        pon=config["refs"]["panel_of_normals"],
    threads: config["resources"]["threads"]
    resources:
        memory_max_gb=config["resources"]["memory_max_gb"],
        memory_min_gb=config["resources"]["memory_min_gb"],
    log:
        "logs/{run}/{sample}/Mutect2.log",
    singularity:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        Mutect2 \
        --native-pair-hmm-threads {threads} \
        -R {params.refdir}/{input.refg} \
        --germline-resource {params.germline} \
        --intervals {input.regions} \
        -I {input.tumor} \
        {'' if params.normal_name == params.default_normal else '-I ' + input.normal} \
        -normal {'null' if params.normal_name == params.default_normal else params.normal_name} \
        {'-pon ' + params.refdir + '/' + params.pon if params.normal_name == params.default_normal else ''} \
        -O {output.vcf} \
        > {log} 2>&1
        """


rule filter_mutect2_calls:
    input:
        vcf="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf",
        idx="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.idx",
        stats="vcf/{run}/{sample}/mutect/{sample}.mutect2.unfiltered.vcf.stats",
        refg=config["refs"]["genome_human"],
    output:
        vcf="vcf/{run}/{sample}/mutect/{sample}.mutect2.filtered.vcf",
        idx="vcf/{run}/{sample}/mutect/{sample}.mutect2.filtered.vcf.idx",
        filtering_stats="metrics/{run}/{sample}/{sample}.mutect2.filtered.filteringStats.tsv",
    params:
        refdir=config["refs"]["path"],
    resources:
        memory_max_gb=config["resources"]["memory_max_gb"],
        memory_min_gb=config["resources"]["memory_min_gb"],
    log:
        "logs/{run}/{sample}/FilterMutectCalls.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.memory_min_gb}G -Xmx{resources.memory_max_gb}G" \
        FilterMutectCalls \
        -R {params.refdir}/{input.refg} \
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
        "../envs/bcftools.yaml"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        bcftools view -R {input.bed} {output.vcf_gz} | bcftools sort -Ov -o {output.vcf}
        """
