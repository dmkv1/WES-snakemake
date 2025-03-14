rule build_xengsort_index:
    input:
        refg_human=config["paths"]["refs"]["genome_human"],
        refg_host=config["paths"]["refs"]["genome_host"],
    output:
        info="refs/xengsort/xengsort-index.info",
        hash="refs/xengsort/xengsort-index.hash",
    params:
        index="refs/xengsort/xengsort-index",
        size=config["params"]["xengsort"]["index_size"],
        fill=config["params"]["xengsort"]["index_fill"],
        k=config["params"]["xengsort"]["index_k"],
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    log:
        "logs/xengsort/xengsort-index.log",
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
        index_info="refs/xengsort/xengsort-index.info",
        index_hash="refs/xengsort/xengsort-index.hash",
        fq1=get_fastq1,
        fq2=get_fastq2,
    output:
        graft1=temp("fastq/{run}/{sample}/{sample}.xengsort-graft.1.fq.gz"),
        graft2=temp("fastq/{run}/{sample}/{sample}.xengsort-graft.2.fq.gz"),
    params:
        index="refs/xengsort/xengsort-index",
        outprefix="fastq/{run}/{sample}/{sample}.xengsort",
        chunksize=config["params"]["xengsort"]["chunksize"],
        prefetch=config["params"]["xengsort"]["prefetch"],
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    log:
        "logs/{run}/{sample}/xengsort.log",
    shell:
        """
        xengsort classify --out {params.outprefix} \
            --index {params.index} \
            --threads {threads} \
            --fastq {input.fq1} --pairs {input.fq2} \
            --chunksize {params.chunksize} \
            --prefetchlevel {params.prefetch} &> {log}
        """
