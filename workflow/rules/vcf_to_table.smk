rule vcf_to_table:
    input:
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.filtered.vcf",
    output:
        xlsx="tables/{run}/{sample}/{sample}.snv_indels.xlsx",
    script:
        "../scripts/vcf_to_table.R"
