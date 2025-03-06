rule build_xengsort_index:
    input:
        refg_human=config["paths"]["refs"]["genome_human"],
        refg_host=config["paths"]["refs"]["genome_host"],
    output:
        info="refs/xengsort/xengsort-k25.info",
        hash="refs/xengsort/xengsort-k25.hash",
    log:
        "logs/xengsort/xengsort-k25.log",
    params:
        index="refs/xengsort/xengsort-k25",
        size=4500000000,
        fill=0.88,
        k=25,
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        xengsort index --index {params.index} \
            -G {input.refg_human} \
            -H {input.refg_host} \
            -k {params.k} -n {params.size} \
            --fill {params.fill} -W {threads} &> {log}
        """


rule run_xengsort:
    input:
        index_info="refs/xengsort/xengsort-k25.info",
        index_hash="refs/xengsort/xengsort-k25.hash",
        fq1=get_fastq1,
        fq2=get_fastq2,
    output:
        graft1=temp("fastq/{run}/{sample}/{sample}.xengsort-graft.1.fq.gz"),
        graft2=temp("fastq/{run}/{sample}/{sample}.xengsort-graft.2.fq.gz"),
    params:
        index="refs/xengsort/xengsort-k25",
        outprefix="fastq/{run}/{sample}/{sample}.xengsort",
        chunksize=16.0,
        prefetch=1,
    log:
        "logs/xengsort/{run}/{sample}.log",
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        xengsort classify --out {params.outprefix} \
            --index {params.index} \
            --threads {threads} \
            --fastq {input.fq1} --pairs {input.fq2} \
            --chunksize {params.chunksize} \
            --prefetchlevel {params.prefetch} &> {log}
        """
