# Getter functions
def get_fastq1(wildcards):
    run = wildcards.run
    return fastq_dict[run][wildcards.sample]["fq1"]


def get_fastq2(wildcards):
    run = wildcards.run
    return fastq_dict[run][wildcards.sample]["fq2"]


def get_tumor_bams(wildcards):
    return expand("results/bam/{sample}.bam", sample=runs_dict[wildcards.run]["tumors"])
