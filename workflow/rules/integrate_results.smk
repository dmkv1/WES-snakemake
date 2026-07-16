rule combine_results:
    input:
        snv_vcf="results/{run}/{sample}/{sample}.SNV.vcf",
        sv_tsv="work/manta/{run}/{sample}/{sample}.SV.annotated.tsv",
        cnv_cns="work/cnvkit/{run}/{sample}/{sample}.call.cns",
        purity_csv="work/purity/{run}/{sample}/{sample}.purity.csv",
    output:
        xlsx="results/{run}/{sample}/{sample}_results.xlsx",
        snv_csv="work/combine/{run}/{sample}/{sample}_SNVs.csv",
        cnv_csv="work/combine/{run}/{sample}/{sample}_CNVs.csv",
        sv_csv="work/combine/{run}/{sample}/{sample}_SVs.csv",
        qc_csv="work/combine/{run}/{sample}/{sample}_QC.csv",
    params:
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
            f"work/combine/{run}/{sample}/{sample}_SNVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        cnv_csvs=[
            f"work/combine/{run}/{sample}/{sample}_CNVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        sv_csvs=[
            f"work/combine/{run}/{sample}/{sample}_SVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
        qc_csvs=[
            f"work/combine/{run}/{sample}/{sample}_QC.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
    output:
        snv_tsv="results/combined/combined_snvs.tsv",
        cnv_tsv="results/combined/combined_cnvs.tsv",
        sv_tsv="results/combined/combined_svs.tsv",
        qc_tsv="results/combined/combined_qc.tsv",
    params:
        samplesheet=config["samplesheet"],
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_all.R"
