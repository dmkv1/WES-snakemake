rule build_xengsort_index:
    input:
        refg_human=config["refs"]["genome_human"],
        refg_host=config["refs"]["genome_host"],
    output:
        info="work/refs/xengsort/xengsort-index.info",
        hash="work/refs/xengsort/xengsort-index.hash",
    params:
        index="work/refs/xengsort/xengsort-index",
        size=config["params"]["xengsort"]["index_size"],
        fill=config["params"]["xengsort"]["index_fill"],
        k=config["params"]["xengsort"]["index_k"],
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    log:
        "work/logs/xengsort-index.log",
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
        index_info="work/refs/xengsort/xengsort-index.info",
        index_hash="work/refs/xengsort/xengsort-index.hash",
        fq1="work/fastq/{run}/{sample}/{sample}.trimmed.1.fq.gz",
        fq2="work/fastq/{run}/{sample}/{sample}.trimmed.2.fq.gz",
    output:
        graft1=temp("work/fastq/{run}/{sample}/{sample}.xengsort-graft.1.fq.gz"),
        graft2=temp("work/fastq/{run}/{sample}/{sample}.xengsort-graft.2.fq.gz"),
    params:
        index="work/refs/xengsort/xengsort-index",
        outprefix="work/fastq/{run}/{sample}/{sample}.xengsort",
        chunksize=config["params"]["xengsort"]["chunksize"],
        prefetch=config["params"]["xengsort"]["prefetch"],
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    log:
        "work/logs/xengsort_{run}_{sample}.log",
    shell:
        """
        xengsort classify --out {params.outprefix} \
            --index {params.index} \
            --threads {threads} \
            --fastq {input.fq1} --pairs {input.fq2} \
            --chunksize {params.chunksize} \
            --prefetchlevel {params.prefetch} &> {log}
        """
