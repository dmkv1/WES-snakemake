rule funcotator:
    input:
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.merged.dp.vcf",
        refg=config["paths"]["refs"]["genome_human"],
        data_sources=config["paths"]["refs"]["funcotator_data_sources"],
    output:
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.annotated.vcf",
        idx="vcf/{run}/{sample}/somaticseq/{sample}.consensus.annotated.vcf.idx",
    params:
        gatk_image=config["tools"]["gatk"]["image"],
        gatk_ver=config["tools"]["gatk"]["version"],
        ref_path=config["paths"]["refs"]["path"],
        genome_ver=config["params"]["genome_version"],
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "logs/{run}/{sample}/Funcotator.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.gatk_image}:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        Funcotator \
        --reference {input.refg} \
        --ref-version {params.genome_ver} \
        --data-sources-path {input.data_sources} \
        --output-file-format VCF \
        --variant {input.vcf} \
        --output {output.vcf} \
        > {log} 2>&1
        """


rule somaticseq_filter:
    input:
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.annotated.vcf",
    output:
        vcf_final="results/{run}/{sample}/{sample}.snv_indels.vcf",
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} | bcftools sort -Ov -o {output.vcf_final}
        """
