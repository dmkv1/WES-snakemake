rule build_xengsort_index:
    output:
        index_info="reference/xengsort/xengsort-index.info",
        index_hash="reference/xengsort/xengsort-index.hash",
    params:
        index="reference/xengsort/xengsort-index",
        refg_human=get_ref_path(config["refs"]["genome_human"], use_container=False),
        refg_host=get_ref_path(config["refs"]["genome_host"], use_container=False),
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
            -G {params.refg_human} \
            -H {params.refg_host} \
            -k {params.k} -n {params.size} \
            --fill {params.fill} -W {threads} &> {log}
        """


rule run_xengsort:
    input:
        index_info="reference/xengsort/xengsort-index.info",
        index_hash="reference/xengsort/xengsort-index.hash",
        fq1=get_fastq1,
        fq2=get_fastq2,
    output:
        graft1=temp("fastq/{run}/{sample}/{sample}.xengsort-graft.1.fq.gz"),
        graft2=temp("fastq/{run}/{sample}/{sample}.xengsort-graft.2.fq.gz"),
        ambiguous1="fastq/{run}/{sample}/{sample}.xengsort-ambiguous.1.fq.gz",
        ambiguous2="fastq/{run}/{sample}/{sample}.xengsort-ambiguous.2.fq.gz",
        both1="fastq/{run}/{sample}/{sample}.xengsort-both.1.fq.gz",
        both2="fastq/{run}/{sample}/{sample}.xengsort-both.2.fq.gz",
        host1="fastq/{run}/{sample}/{sample}.xengsort-host.1.fq.gz",
        host2="fastq/{run}/{sample}/{sample}.xengsort-host.2.fq.gz",
        neither1="fastq/{run}/{sample}/{sample}.xengsort-neither.1.fq.gz",
        neither2="fastq/{run}/{sample}/{sample}.xengsort-neither.2.fq.gz",
    params:
        index="reference/xengsort/xengsort-index",
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
