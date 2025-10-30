rule combine_results:
    input:
        snv_vcf="results/{run}/{sample}/{sample}.SNV.vcf",
        sv_tsv="work/manta/{run}/{sample}/{sample}.SV.annotated.tsv",
        cnv_cns="work/cnvkit/{run}/{sample}/{sample}.call.cns"
    output:
        xlsx="results/{run}/{sample}/{sample}_results.xlsx",
    params:
        purity=get_purity,
        sample_sex=get_sample_sex,
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_results.R"