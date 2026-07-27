rule merge_sv_callers:
    # Consensus of Manta + DELLY via SURVIVOR. The merged VCF carries SUPP and
    # SUPP_VEC in INFO; SUPP_VEC bit order follows the file-list order below
    # (bit 1 = Manta, bit 2 = DELLY), so each SV records which caller(s) found it.
    input:
        manta="work/manta/{run}/{sample}/{sample}.SV.filtered.vcf",
        delly="work/delly/{run}/{sample}/{sample}.SV.filtered.vcf",
    output:
        vcf="work/sv_merge/{run}/{sample}/{sample}.SV.merged.vcf",
    params:
        file_list="work/sv_merge/{run}/{sample}/{sample}.sv_vcf_list.txt",
        max_dist=config["params"]["delly"]["survivor_max_dist"],
        min_size=config["params"]["delly"]["survivor_min_size"],
    conda:
        "../envs/survivor.yaml"
    log:
        "work/logs/SurvivorMerge_{run}_{sample}.log",
    shell:
        """
        # Bit order in SUPP_VEC follows this list: Manta then DELLY.
        # Args: max_dist, min_support=1 (union), type=1 (require same SVTYPE),
        # strand=0 (ignore strand: Manta and DELLY encode BND strands
        # differently, so strand matching under-merges across callers),
        # estimate_dist=0, min_size.
        printf '%s\\n%s\\n' {input.manta} {input.delly} > {params.file_list}
        SURVIVOR merge {params.file_list} {params.max_dist} 1 1 0 0 {params.min_size} {output.vcf} \
        > {log} 2>&1
        """


rule filter_merged_sv_size:
    # Drop spurious multi-Mb intra-chromosomal SVs (see params.delly.max_sv_size).
    # Translocations (BND/TRA) have no meaningful SVLEN and are always kept.
    input:
        vcf="work/sv_merge/{run}/{sample}/{sample}.SV.merged.vcf",
    output:
        vcf="work/sv_merge/{run}/{sample}/{sample}.SV.merged.sizefilt.vcf",
    params:
        max_size=config["params"]["delly"]["max_sv_size"],
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        if [ "{params.max_size}" -eq 0 ]; then
            cp {input.vcf} {output.vcf}
        else
            bcftools view \
            -e '(SVTYPE!="BND" && SVTYPE!="TRA") && (SVLEN>{params.max_size} || SVLEN<-{params.max_size})' \
            {input.vcf} -Ov -o {output.vcf}
        fi
        """


rule annotate_merged_sv:
    # AnnotSV over the consensus VCF. SURVIVOR's SUPP/SUPP_VEC ride through in the
    # INFO column. Empty-VCF guard mirrors the former Manta annotate_sv.
    input:
        vcf="work/sv_merge/{run}/{sample}/{sample}.SV.merged.sizefilt.vcf",
    output:
        tsv="work/sv_merge/{run}/{sample}/{sample}.SV.annotated.tsv",
    params:
        annotations=config["refs"]["annotsv_annotations"]["path"],
        genome_build=config["refs"]["annotsv_annotations"]["genome_build"],
        output_prefix="work/sv_merge/{run}/{sample}/{sample}.SV.annotated",
    threads: 1
    conda:
        "../envs/annotsv.yaml"
    log:
        "work/logs/AnnotSV_{run}_{sample}.log",
    shell:
        """
        variant_count=$(bcftools view -H {input.vcf} | wc -l)

        if [ "$variant_count" -eq 0 ]; then
            touch {output.tsv}
            echo "No variants to annotate" > {log}
        else
            AnnotSV \
            -SVinputFile {input.vcf} \
            -outputFile {params.output_prefix} \
            -genomeBuild {params.genome_build} \
            -annotationsDir {params.annotations} \
            > {log} 2>&1
        fi
        """
