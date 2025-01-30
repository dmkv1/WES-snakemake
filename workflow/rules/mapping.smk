rule bwa_map:
    input:
        ref=config["paths"]["refs"]["genome_human"],
        fq1=get_fastq1,
        fq2=get_fastq2
    output:
        "bam/{run}/{sample}.bam",
    threads: 
        config["resources"]["bwa_threads"]
    shell:
        "bwa mem -M -t {threads} {input.ref} {input.fq1} {input.fq2} | samtools view -Sb - > {output}"
