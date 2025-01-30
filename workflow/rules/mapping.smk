rule bwa_map:
    input:
        ref=config["paths"]["refs"]["genome_human"],
        fq1=get_fastq1,
        fq2=get_fastq2
    output:
        "bam/{run}/{sample}.raw.bam",
    threads: 
        config["resources"]["bwa_threads"]
    shell:
        "bwa mem -M -t {threads} {input.ref} {input.fq1} {input.fq2} | samtools view -Sb - > {output}"

rule add_read_groups:
    input:
        "bam/{run}/{sample}.raw.bam"
    output:
        "bam/{run}/{sample}.rg.bam"
    params:
        rg=get_read_group_params,
        tmp_dir="tmp",
        gatk_ver=config["tools"]["gatk_version"]
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    container:
        "docker://broadinstitute/gatk:{params.gatk_ver}"
    shell:
        """
        docker run --rm \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        broadinstitute/gatk:{params.gatk_ver} gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        AddOrReplaceReadGroups \
        -I {input} \
        -O {output} \
        -ID {params.rg[RGID]} \
        -SM {params.rg[RGSM]} \
        -PU {params.rg[RGPU]} \
        -LB {params.rg[RGLB]} \
        -PL {params.rg[RGPL]} \
        -TMP_DIR {params.tmp_dir}
        """

