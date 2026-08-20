rule run_delly:
    # Orthogonal SV caller. DELLY needs sorted, indexed, dup-marked BAMs (the
    # final analysis BAM already is), so unlike Manta no base-quality capping is
    # applied. Paired runs pass tumor + matched normal (tumor first); tumor-only
    # runs call on the tumor alone.
    input:
        tumor="results/{run}/{sample}/bam/{sample}.bam",
        tumor_bai="results/{run}/{sample}/bam/{sample}.bai",
        normal=lambda w: (
            f"results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bam"
            if not is_tumor_only(w) else []
        ),
        normal_bai=lambda w: (
            f"results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bai"
            if not is_tumor_only(w) else []
        ),
        refg=config["refs"]["genome_human"],
    output:
        bcf="work/delly/{run}/{sample}/{sample}.delly.bcf",
        csi="work/delly/{run}/{sample}/{sample}.delly.bcf.csi",
    params:
        normal_arg=lambda w: (
            f"results/{w.run}/{runs_dict[w.run]['normal']}/bam/{runs_dict[w.run]['normal']}.bam"
            if not is_tumor_only(w) else ""
        ),
        exclude_arg=lambda w: (
            f"-x {config['params']['delly']['exclude']}"
            if config["params"]["delly"]["exclude"] else ""
        ),
        map_qual=config["params"]["delly"]["map_qual"],
        mad_cutoff=config["params"]["delly"]["mad_cutoff"],
        min_clique_size=config["params"]["delly"]["min_clique_size"],
    # DELLY's OpenMP parallelism is over the input BAMs, so it can use one thread
    # per sample and no more: two for a paired run, one tumor-only.
    threads: lambda wildcards: 1 if is_tumor_only(wildcards) else 2
    resources:
        mem_mb=config["resources"]["delly_mem_mb"],
    benchmark:
        "work/benchmarks/run_delly/{run}_{sample}.tsv",
    container:
        config["containers"]["delly"]
    log:
        "work/logs/Delly_{run}_{sample}.log",
    shell:
        """
        export OMP_NUM_THREADS={threads}
        delly sr \
        -g {input.refg} \
        {params.exclude_arg} \
        -q {params.map_qual} \
        -s {params.mad_cutoff} \
        -z {params.min_clique_size} \
        -o {output.bcf} \
        {input.tumor} {params.normal_arg} \
        > {log} 2>&1
        """


rule filter_delly:
    # Paired runs: somatic filtering against the matched normal (samples.tsv maps
    # the tumor sample to 'tumor' and the normal to 'control', both by their VCF
    # sample IDs, which equal the sample names since RGSM == sample). Tumor-only
    # runs: no somatic filter is possible, so the raw germline calls pass through.
    input:
        bcf="work/delly/{run}/{sample}/{sample}.delly.bcf",
        csi="work/delly/{run}/{sample}/{sample}.delly.bcf.csi",
    output:
        bcf="work/delly/{run}/{sample}/{sample}.delly.filtered.bcf",
    params:
        tumor_only=lambda w: is_tumor_only(w),
        tumor=lambda w: w.sample,
        normal=lambda w: runs_dict[w.run]["normal"] or "",
        samples_tsv="work/delly/{run}/{sample}/{sample}.samples.tsv",
    container:
        config["containers"]["delly"]
    log:
        "work/logs/DellyFilter_{run}_{sample}.log",
    shell:
        """
        if [ "{params.tumor_only}" = "True" ]; then
            cp {input.bcf} {output.bcf}
            echo "Tumor-only: no somatic filter applied" > {log}
        else
            printf '%s\\ttumor\\n%s\\tcontrol\\n' {params.tumor} {params.normal} > {params.samples_tsv}
            delly filter -f somatic -s {params.samples_tsv} -o {output.bcf} {input.bcf} > {log} 2>&1
        fi
        """


rule delly_to_vcf:
    # Convert DELLY BCF to VCF, keep PASS, and restrict to the capture-kit target
    # regions (DELLY has no exome mode; Manta is already exome-restricted via
    # --callRegions, so this keeps the consensus in target space).
    input:
        bcf="work/delly/{run}/{sample}/{sample}.delly.filtered.bcf",
        regions_bed=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz",
        regions_tbi=lambda w: f"work/refs/regions/{get_probe_version(w)}/regions.bed.gz.tbi",
    output:
        vcf="work/delly/{run}/{sample}/{sample}.SV.filtered.vcf",
    params:
        tmp_vcf="work/delly/{run}/{sample}/{sample}.delly.tmp.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view {input.bcf} -Oz -o {params.tmp_vcf}
        tabix -p vcf {params.tmp_vcf}
        bcftools view -f PASS -R {input.regions_bed} {params.tmp_vcf} \
        | bcftools sort -Ov -o {output.vcf}
        rm -f {params.tmp_vcf} {params.tmp_vcf}.tbi
        """
