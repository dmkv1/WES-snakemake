rule bwa_map_host:
    input:
        refg=config["paths"]["refs"]["genome_host"],
        fq1=get_fastq1,
        fq2=get_fastq2,
    output:
        temp("bam/{run}/{sample}.host.unsorted.bam"),
    threads: config["resources"]["threads"]
    shell:
        "bwa mem -M -t {threads} {input.refg} {input.fq1} {input.fq2} | samtools view -Sb - > {output}"

rule sort_host_bam:
    input:
        "bam/{run}/{sample}.host.unsorted.bam"
    output:
        temp("bam/{run}/{sample}.host.bam")
    threads: config["resources"]["threads"]
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule run_xenofilter:
    input:
        graft_bam="bam/{run}/{sample}.md.bam",
        host_bam="bam/{run}/{sample}.host.bam"
    output:
        input_graft_bai=temp("bam/{run}/{sample}.md.bam.bai"),
        bam=temp("bam/{run}/{sample}.xenofiltered.bam"),
        bai=temp("bam/{run}/{sample}.xenofiltered.bam.bai"),
        log="metrics/{run}/{sample}.xenofilter.log"
    script:
        "../scripts/run_xenofilter.R"
