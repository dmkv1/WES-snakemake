rule bwa_map_mouse:
    input:
        refg=config["paths"]["refs"]["genome_mouse"],
        fq1=get_fastq1,
        fq2=get_fastq2,
    output:
        temp("bam/{run}/{sample}.mouse.unsorted.bam"),
    threads: config["resources"]["bwa_threads"]
    shell:
        "bwa mem -M -t {threads} {input.refg} {input.fq1} {input.fq2} | samtools view -Sb - > {output}"

rule sort_mouse_bam:
    input:
        "bam/{run}/{sample}.mouse.unsorted.bam"
    output:
        temp("bam/{run}/{sample}.mouse.bam")
    threads: config["resources"]["bwa_threads"]
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule run_xenofilter:
    input:
        human_bam="bam/{run}/{sample}.md.bam",
        mouse_bam="bam/{run}/{sample}.mouse.bam"
    output:
        temp("bam/{run}/{sample}.xenofiltered.bam")
    script:
        "../scripts/run_xenofilter.R"
