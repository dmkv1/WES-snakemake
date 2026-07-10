rule combine_results:
    input:
        snv_vcf="results/{run}/{sample}/{sample}.SNV.vcf",
        sv_tsv="work/manta/{run}/{sample}/{sample}.SV.annotated.tsv",
        cnv_cns="work/cnvkit/{run}/{sample}/{sample}.call.cns"
    output:
        xlsx="results/{run}/{sample}/{sample}_results.xlsx",
        snv_csv="results/{run}/{sample}/{sample}_SNVs.csv",
        cnv_csv="results/{run}/{sample}/{sample}_CNVs.csv",
        sv_csv="results/{run}/{sample}/{sample}_SVs.csv",
    params:
        purity=get_purity,
        sample_sex=get_sample_sex,
        normal=lambda w: runs_dict[w.run]["normal"] if runs_dict[w.run]["normal"] else "",
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_results.R"


rule merge_results:
    """Row-bind every sample's per-type CSVs into three cross-cohort TSVs."""
    input:
        snv_csvs=[
            f"results/{run}/{sample}/{sample}_SNVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        cnv_csvs=[
            f"results/{run}/{sample}/{sample}_CNVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        sv_csvs=[
            f"results/{run}/{sample}/{sample}_SVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
    output:
        snv_tsv="results/combined/combined_snvs.tsv",
        cnv_tsv="results/combined/combined_cnvs.tsv",
        sv_tsv="results/combined/combined_svs.tsv",
    params:
        samplesheet=config["samplesheet"],
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_all.R"
