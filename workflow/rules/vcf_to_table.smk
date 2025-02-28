rule vcf_to_table:
    input:
        vcf="results/{run}/{sample}/{sample}.snv_indels.vcf",
    output:
        xlsx="results/{run}/{sample}/{sample}.snv_indels.xlsx",
    script:
        "../scripts/vcf_to_table.R"
