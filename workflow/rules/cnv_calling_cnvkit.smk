def get_cnvkit_reference(wildcards):
    probe_version = probe_dict[wildcards.run][wildcards.sample]
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    gender = sample_row["gender"].iloc[0]

    ref_key = "cnvkit_ref_m" if gender == "m" else "cnvkit_ref_f"
    return config["panel_of_normals"]["cnvkit"][ref_key][probe_version]


def get_normal_sample_name(wildcards):
    """Get normal sample name, returns None for tumor-only"""
    return runs_dict[wildcards.run]["normal"]


def get_normal_sample_name_or_empty(wildcards):
    """Get normal sample name, returns empty string for tumor-only (for shell params)"""
    normal = runs_dict[wildcards.run]["normal"]
    return normal if normal else ""


def get_sample_sex(wildcards):
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    gender = sample_row["gender"].iloc[0]
    return "male" if gender == "m" else "female"


def is_male_reference(wildcards):
    return "--male-reference" if get_sample_sex(wildcards) == "male" else ""


rule cnvkit_coverage_and_fix:
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        targets=lambda w: config["probe_configs"][get_probe_version(w)][
            "cnvkit_targets"
        ],
        antitargets=lambda w: config["probe_configs"][get_probe_version(w)][
            "cnvkit_antitargets"
        ],
        reference=get_cnvkit_reference,
    output:
        cnr="work/cnvkit/{run}/{sample}/{sample}.cnr",
    params:
        prefix="work/cnvkit/{run}/{sample}/{sample}",
    threads: config["resources"]["threads"]
    resources:
        # pysam BAM traversal, one worker process per -p thread.
        mem_mb=8192,
    benchmark:
        "work/benchmarks/cnvkit_coverage_and_fix/{run}_{sample}.tsv",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_coverage_fix_{run}_{sample}.log",
    shell:
        """
        cnvkit.py coverage {input.bam} {input.targets} \
            -o {params.prefix}.targetcoverage.cnn \
            -p {threads} > {log} 2>&1
        
        cnvkit.py coverage {input.bam} {input.antitargets} \
            -o {params.prefix}.antitargetcoverage.cnn \
            -p {threads} >> {log} 2>&1
        
        cnvkit.py fix {params.prefix}.targetcoverage.cnn \
            {params.prefix}.antitargetcoverage.cnn \
            {input.reference} \
            -o {output.cnr} >> {log} 2>&1
        
        rm {params.prefix}.targetcoverage.cnn {params.prefix}.antitargetcoverage.cnn
        """


def _is_paired_tumor(w):
    """True for a tumour sample in paired mode (i.e. not tumour-only and not
    the run's matched normal)."""
    return not is_tumor_only(w) and w.sample != runs_dict[w.run]["normal"]


rule gatk_haplotypecaller:
    # Germline SNP discovery. Only ever requested for the run's matched normal:
    # it genotypes the exome so extract_normal_het_sites can select het BAF
    # anchors. The paired tumour no longer runs HaplotypeCaller — its ref/alt
    # depth at these known sites comes from CollectAllelicCounts (pileup, no
    # assembly) via tumor_baf_allelic_counts / tumor_baf_vcf (#24). Tumour-only
    # runs never request this output (no BAF track).
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda w: f"work/refs/regions/{probe_dict[w.run][w.sample]}/regions.bed",
    output:
        vcf="work/haplotypecaller/{run}/{sample}/{sample}.germline.vcf",
        idx="work/haplotypecaller/{run}/{sample}/{sample}.germline.vcf.idx",
    params:
        ref_path=config["refs"]["path"],
    threads: config["resources"]["threads"]
    resources:
        java_max_gb=config["resources"]["gatk"]["medium"]["java_max_gb"],
        java_min_gb=config["resources"]["gatk"]["medium"]["java_min_gb"],
        mem_mb=config["resources"]["gatk"]["medium"]["mem_mb"],
    benchmark:
        "work/benchmarks/gatk_haplotypecaller/{run}_{sample}.tsv",
    log:
        "work/logs/HaplotypeCaller_{run}_{sample}.log",
    container:
        config["containers"]["gatk"]
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        HaplotypeCaller \
        -R {input.refg} \
        -I {input.bam} \
        -L {input.regions} \
        -O {output.vcf} \
        --native-pair-hmm-threads {threads} \
        > {log} 2>&1
        """


rule tumor_baf_allelic_counts:
    # #24: tumour-side BAF depth at the matched normal's het sites via pileup
    # counting instead of HaplotypeCaller's per-active-region re-assembly + HMM.
    # We only need ref/alt depth at *known* positions, which is exactly what
    # CollectAllelicCounts (GATK's own allelic-CN tool) provides — much cheaper
    # than force-calling with HaplotypeCaller --alleles. Paired tumour only.
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        refg=config["refs"]["genome_human"],
        het_sites="work/haplotypecaller/{run}/normal_het_sites.vcf.gz",
        het_sites_tbi="work/haplotypecaller/{run}/normal_het_sites.vcf.gz.tbi",
    output:
        tsv="work/haplotypecaller/{run}/{sample}/{sample}.allelicCounts.tsv",
    resources:
        java_max_gb=config["resources"]["gatk"]["light"]["java_max_gb"],
        java_min_gb=config["resources"]["gatk"]["light"]["java_min_gb"],
        mem_mb=config["resources"]["gatk"]["light"]["mem_mb"],
    benchmark:
        "work/benchmarks/tumor_baf_allelic_counts/{run}_{sample}.tsv",
    log:
        "work/logs/CollectAllelicCounts_{run}_{sample}.log",
    container:
        config["containers"]["gatk"]
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        CollectAllelicCounts \
        -R {input.refg} \
        -I {input.bam} \
        -L {input.het_sites} \
        --minimum-base-quality 20 \
        -O {output.tsv} \
        > {log} 2>&1
        """


rule tumor_baf_vcf:
    # #24: convert CollectAllelicCounts output into the single-sample tumour VCF
    # that cnvkit_merge_germline_and_filter_hetsnp expects in place of the old
    # HaplotypeCaller force-call output. REF/ALT and contig header are taken
    # verbatim from the normal het-site VCF so bcftools merge aligns records
    # exactly; AD/DP come from the pileup counts; GT is a genuine call from the
    # observed BAF (mirrors the shape of the former force-call: every record
    # carries a real genotype + AD).
    input:
        het_sites="work/haplotypecaller/{run}/normal_het_sites.vcf.gz",
        tsv="work/haplotypecaller/{run}/{sample}/{sample}.allelicCounts.tsv",
    output:
        vcf="work/haplotypecaller/{run}/{sample}/{sample}.tumor_baf.vcf",
    params:
        sample=lambda w: w.sample,
    run:
        import gzip

        # Normal het sites: authoritative REF/ALT + contig header lines, in order.
        contigs = []
        sites = []  # (chrom, pos, ref, alt)
        with gzip.open(input.het_sites, "rt") as fh:
            for line in fh:
                if line.startswith("##contig="):
                    contigs.append(line.rstrip("\n"))
                elif line.startswith("#"):
                    continue
                else:
                    f = line.split("\t")
                    sites.append((f[0], f[1], f[3], f[4]))

        # CollectAllelicCounts: skip SAM (@) header and the column header row.
        counts = {}  # (chrom, pos) -> (ref_count, alt_count, ref_nuc, alt_nuc)
        with open(input.tsv) as fh:
            for line in fh:
                if line.startswith("@") or line.startswith("CONTIG\t"):
                    continue
                c = line.rstrip("\n").split("\t")
                counts[(c[0], c[1])] = (int(c[2]), int(c[3]), c[4], c[5])

        def genotype(ref_ad, alt_ad):
            dp = ref_ad + alt_ad
            if dp == 0:
                return "./."
            baf = alt_ad / dp
            if baf < 0.1:
                return "0/0"
            if baf > 0.9:
                return "1/1"
            return "0/1"

        with open(output.vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            for c in contigs:
                out.write(c + "\n")
            out.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            )
            out.write(
                '##FORMAT=<ID=AD,Number=R,Type=Integer,'
                'Description="Allelic depths for the ref and alt alleles">\n'
            )
            out.write(
                '##FORMAT=<ID=DP,Number=1,Type=Integer,'
                'Description="Approximate read depth">\n'
            )
            out.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + params.sample
                + "\n"
            )
            for chrom, pos, ref, alt in sites:
                rc, ac, ref_nuc, alt_nuc = counts.get((chrom, pos), (0, 0, ref, alt))
                ref_ad = rc
                # CollectAllelicCounts reports the single highest-count alt base;
                # only credit it when it is the normal's alt (at near-hom sites a
                # noise base may be reported, but its count is ~0 either way).
                alt_ad = ac if alt_nuc == alt else 0
                dp = ref_ad + alt_ad
                gt = genotype(ref_ad, alt_ad)
                out.write(
                    "\t".join(
                        [
                            chrom,
                            pos,
                            ".",
                            ref,
                            alt,
                            ".",
                            ".",
                            ".",
                            "GT:AD:DP",
                            f"{gt}:{ref_ad},{alt_ad}:{dp}",
                        ]
                    )
                    + "\n"
                )


rule extract_normal_het_sites:
    """Heterozygous biallelic SNP sites from the matched normal. Seeds the
    tumour CollectAllelicCounts -L (#24) so BAF sites survive in CNV/LoH
    regions. Mirrors the normal-side criteria of
    cnvkit_merge_germline_and_filter_hetsnp."""
    input:
        normal=lambda w: (
            f"work/haplotypecaller/{w.run}/{runs_dict[w.run]['normal']}/"
            f"{runs_dict[w.run]['normal']}.germline.vcf"
        ),
    output:
        vcf="work/haplotypecaller/{run}/normal_het_sites.vcf.gz",
        tbi="work/haplotypecaller/{run}/normal_het_sites.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    log:
        "work/logs/extract_normal_het_sites_{run}.log",
    shell:
        """
        bcftools filter -i "
                TYPE='snp' &&
                N_ALT==1 &&
                STRLEN(REF)==1 &&
                STRLEN(ALT)==1 &&
                FORMAT/DP[0]>10 &&
                FORMAT/AD[0:1]>3 &&
                GT[0]=='0/1'" {input.normal} -Oz -o {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule cnvkit_merge_germline_and_filter_hetsnp:
    input:
        normal=lambda w: (
            f"work/haplotypecaller/{w.run}/{runs_dict[w.run]['normal']}/{runs_dict[w.run]['normal']}.germline.vcf"
            if not is_tumor_only(w) else []
        ),
        tumor=lambda w: (
            f"work/haplotypecaller/{w.run}/{w.sample}/{w.sample}.tumor_baf.vcf"
            if not is_tumor_only(w) else []
        ),
    output:
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    params:
        normal=lambda w: runs_dict[w.run]["normal"] if not is_tumor_only(w) else "",
        tumor=lambda w: w.sample,
        is_tumor_only=lambda w: is_tumor_only(w),
    conda:
        "../envs/bcftools.yaml"
    log:
        "work/logs/merge_filter_germline_{run}_{sample}.log",
    shell:
        """
        if [ "{params.is_tumor_only}" = "True" ]; then
            # Tumor-only: create empty VCF (no BAF analysis)
            echo "##fileformat=VCFv4.2" > {output.vcf}
            echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> {output.vcf}
            echo "Tumor-only mode: skipping BAF analysis" > {log}
        else
            # Paired mode: bgzip into sample-specific temp files to avoid race conditions
            # when multiple tumors share the same control
            TMPDIR=$(dirname {output.vcf})
            NORMAL_GZ="$TMPDIR/{params.normal}.germline.vcf.gz"
            TUMOR_GZ="$TMPDIR/{params.tumor}.tumor_baf.vcf.gz"

            bgzip -c {input.normal} > "$NORMAL_GZ" 2> {log}
            tabix -f -p vcf "$NORMAL_GZ" 2>> {log}

            bgzip -c {input.tumor} > "$TUMOR_GZ" 2>> {log}
            tabix -f -p vcf "$TUMOR_GZ" 2>> {log}

            bcftools merge -m none "$NORMAL_GZ" "$TUMOR_GZ" -Ou 2>> {log} | \
            bcftools view -s {params.normal},{params.tumor} -Ou | \
            bcftools filter -i "
                    TYPE='snp' &&
                    N_ALT==1 &&
                    STRLEN(REF)==1 &&
                    STRLEN(ALT)==1 &&
                    FORMAT/DP[0]>10 &&
                    FORMAT/AD[0:1]>3 &&
                    GT[0]=='0/1' &&
                    FORMAT/DP[1]>10" -Ov -o {output.vcf} 2>> {log}
        fi
        """


rule filter_chromosomes:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.cnr",
    output:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
    params:
        chromosomes=[str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"],
    run:
        def check_chromosome(chrom, allowed_chroms):
            """Check if chromosome (with or without chr prefix) is allowed"""
            return chrom.replace("chr", "") in allowed_chroms


        def filter_file(file_in, file_out, allowed_chroms):
            """Filter file keeping only allowed chromosomes"""
            with open(file_in) as in_handle, open(file_out, "w") as out_handle:
                header = True
                for line in in_handle:
                    if header:
                        out_handle.write(line)
                        header = False
                    else:
                        chrom = line.split("\t")[0]
                        if check_chromosome(chrom, allowed_chroms):
                            out_handle.write(line)


        filter_file(input.cnr, output.cnr, params.chromosomes)


rule cnvkit_segment:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.cns",
    params:
        normal=get_normal_sample_name_or_empty,
        tumor=lambda w: w.sample,
        vcf_arg=lambda w: f"-v work/cnvkit/{w.run}/{w.sample}/{w.sample}.hetsnp.vcf" if not is_tumor_only(w) else "",
        normal_arg=lambda w: f"-n {runs_dict[w.run]['normal']}" if not is_tumor_only(w) else "",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_segment_{run}_{sample}.log",
    shell:
        """
        cnvkit.py segment {input.cnr} \
            -o {output.cns} \
            -m cbs \
            {params.vcf_arg} \
            {params.normal_arg} \
            -i {params.tumor} \
            > {log} 2>&1
        """


rule cnvkit_segmetrics:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns="work/cnvkit/{run}/{sample}/{sample}.cns",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.segmetrics.cns",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_segmetrics_{run}_{sample}.log",
    shell:
        """
        cnvkit.py segmetrics {input.cnr} \
            -s {input.cns} \
            --ci --alpha 0.5 --smooth-bootstrap 10 \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_filter_ci:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.segmetrics.cns",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.ci_filtered.cns",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_filter_ci_{run}_{sample}.log",
    shell:
        """
        cnvkit.py call {input.cns} \
            --method none \
            --filter ci \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_add_ttest:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns=lambda w: (
            f"work/cnvkit/{w.run}/{w.sample}/{w.sample}.ci_filtered.cns"
            if config["params"]["cnvkit"]["filter_ci"]
            else f"work/cnvkit/{w.run}/{w.sample}/{w.sample}.segmetrics.cns"
        ),
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.ttest.cns",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_ttest_{run}_{sample}.log",
    shell:
        """
        cnvkit.py segmetrics {input.cnr} \
            -s {input.cns} \
            --t-test \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_call:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.ttest.cns",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
        purity_csv="work/purity/{run}/{sample}/{sample}.purity.csv",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.call.cns",
    params:
        normal=get_normal_sample_name_or_empty,
        tumor=lambda w: w.sample,
        sample_sex=get_sample_sex,
        male_reference=is_male_reference,
        purity_ploidy_args=get_purity_ploidy_args,
        vcf_arg=lambda w: f"-v work/cnvkit/{w.run}/{w.sample}/{w.sample}.hetsnp.vcf" if not is_tumor_only(w) else "",
        normal_arg=lambda w: f"-n {runs_dict[w.run]['normal']}" if not is_tumor_only(w) else "",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_call_{run}_{sample}.log",
    shell:
        """
        cnvkit.py call {input.cns} \
            {params.vcf_arg} \
            {params.normal_arg} \
            -i {params.tumor} \
            --sample-sex {params.sample_sex} \
            {params.male_reference} \
            {params.purity_ploidy_args} \
            -m clonal \
            --center median \
            --drop-low-coverage \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_plots:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns="work/cnvkit/{run}/{sample}/{sample}.call.cns",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    output:
        scatter="results/{run}/{sample}/{sample}.scatter.png",
        diagram="results/{run}/{sample}/{sample}.diagram.pdf",
    params:
        normal=get_normal_sample_name_or_empty,
        tumor=lambda w: w.sample,
        vcf_arg=lambda w: f"-v work/cnvkit/{w.run}/{w.sample}/{w.sample}.hetsnp.vcf" if not is_tumor_only(w) else "",
        normal_arg=lambda w: f"-n {runs_dict[w.run]['normal']}" if not is_tumor_only(w) else "",
    container:
        config["containers"]["cnvkit"]
    log:
        "work/logs/cnvkit_plots_{run}_{sample}.log",
    shell:
        """
        cnvkit.py scatter {input.cnr} \
            -s {input.cns} \
            {params.vcf_arg} \
            {params.normal_arg} \
            -i {params.tumor} \
            --title {params.tumor} \
            --y-max 4 --y-min -4 \
            --fig-size 10 5 \
            -o {output.scatter} \
            > {log} 2>&1

        cnvkit.py diagram {input.cnr} \
            -s {input.cns} \
            --no-gene-labels \
            --title {params.tumor} \
            -o {output.diagram} \
            >> {log} 2>&1
        """
