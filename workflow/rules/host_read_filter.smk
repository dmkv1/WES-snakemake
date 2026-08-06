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
        fq1="work/fastq/{run}/{sample}/{sample}.{unit}_R1.fq.gz",
        fq2="work/fastq/{run}/{sample}/{sample}.{unit}_R2.fq.gz",
    output:
        graft1=temp("work/fastq/{run}/{sample}/{sample}.{unit}.xengsort-graft.1.fq.gz"),
        graft2=temp("work/fastq/{run}/{sample}/{sample}.{unit}.xengsort-graft.2.fq.gz"),
    params:
        index="work/refs/xengsort/xengsort-index",
        outprefix="work/fastq/{run}/{sample}/{sample}.{unit}.xengsort",
        chunksize=config["params"]["xengsort"]["chunksize"],
        prefetch=config["params"]["xengsort"]["prefetch"],
        display_name=get_unit_display_name,
    conda:
        "../envs/xengsort.yaml"
    threads: config["resources"]["threads"]
    log:
        "work/logs/xengsort_{run}_{sample}.{unit}.log",
    shell:
        """
        xengsort classify --out {params.outprefix} \
            --index {params.index} \
            --threads {threads} \
            --fastq {input.fq1} --pairs {input.fq2} \
            --chunksize {params.chunksize} \
            --prefetchlevel {params.prefetch} &> {log}
        # MultiQC's xengsort module names the report row after the --out path printed in
        # this log, not the log filename, so the name is rewritten here. It carries the
        # unit only when the sample has more than one, matching how every other QC row
        # is named. Because that name comes from the log content and not a filename,
        # MultiQC sample renaming never sees it.
        # Braces are avoided in this comment on purpose: Snakemake formats the whole
        # shell block, comments included, and a bare placeholder is a NameError.
        sed -i "s#{params.outprefix}#{params.display_name}#g" {log}
        """
