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
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_results.R"


rule combine_all_snvs:
    input:
        snv_csvs=[
            f"results/{run}/{sample}/{sample}_SNVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
    output:
        rds="results/combined/combined_snvs.rds",
    params:
        samplesheet=config["samplesheet"],
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_snvs.R"


rule combine_all_cnvs:
    input:
        cnv_csvs=[
            f"results/{run}/{sample}/{sample}_CNVs.csv"
            for run in runs_dict
            for sample in runs_dict[run]["tumors"]
        ],
    output:
        rds="results/combined/combined_cnvs.rds",
    params:
        samplesheet=config["samplesheet"],
    conda:
        "../envs/r_vcf.yaml"
    script:
        "../scripts/combine_cnvs.R"