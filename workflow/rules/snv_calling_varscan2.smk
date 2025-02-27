rule samtools_mpileup:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        refg=config["paths"]["refs"]["genome_human"],
    output:
        mpileup=temp("mpileup/{run}/{sample}.mpileup"),
    conda:
        "../envs/varscan.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        samtools mpileup \
        -f {input.refg} \
        -q 1 \
        -B \
        {input.normal} {input.tumor} > {output.mpileup}
        """


rule varscan_somatic:
    input:
        mpileup="mpileup/{run}/{sample}.mpileup",
    output:
        snp=temp("vcf/{run}/{sample}/varscan/{sample}.varscan.snp.vcf"),
        indel="vcf/{run}/{sample}/varscan/{sample}.varscan.indel.vcf",
    params:
        min_coverage=config["params"]["varscan_min_coverage"],
        min_var_freq=config["params"]["varscan_min_var_freq"],
        p_value=config["params"]["varscan_p_value"],
        somatic_p_value=config["params"]["varscan_somatic_p_value"],
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        varscan somatic \
        {input.mpileup} \
        --output-snp {output.snp} \
        --output-indel {output.indel} \
        --mpileup 1 \
        --min-coverage {params.min_coverage} \
        --min-var-freq {params.min_var_freq} \
        --p-value {params.p_value} \
        --somatic-p-value {params.somatic_p_value} \
        --strand-filter 1 \
        --output-vcf 1
        """


rule varscan_filter:
    input:
        vcf_snv="vcf/{run}/{sample}/varscan/{sample}.varscan.snp.vcf",
        vcf_indel="vcf/{run}/{sample}/varscan/{sample}.varscan.indel.vcf",
    output:
        vcf="vcf/{run}/{sample}/varscan/{sample}.varscan.filtered_snp.vcf",
    params:
        min_coverage=config["params"]["varscan_filter_min_coverage"],
        min_reads=config["params"]["varscan_filter_min_reads"],
        min_strands=config["params"]["varscan_filter_min_strands"],
        min_avg_qual=config["params"]["varscan_filter_min_avg_qual"],
        min_var_freq=config["params"]["varscan_filter_min_var_freq"],
        p_value=config["params"]["varscan_filter_p_value"],
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        varscan somaticFilter \
        {input.vcf_snv} \
        --output-file {output.vcf} \
        --min-coverage {params.min_coverage} \
        --min-reads2 {params.min_reads} \
        --min-strands2 {params.min_strands} \
        --min-avg-qual {params.min_avg_qual} \
        --min-var-freq {params.min_var_freq} \
        --p-value {params.p_value} \
        --indel-file {input.vcf_indel}
        """


rule sort_varscan_calls:
    input:
        vcf="vcf/{run}/{sample}/varscan/{sample}.varscan.filtered_snp.vcf",
        bed="refs/regions/regions.bed.gz",
    output:
        vcf_gz=temp("vcf/{run}/{sample}/varscan/{sample}.varscan.filtered.vcf.gz"),
        vcf="vcf/{run}/{sample}/varscan/{sample}.varscan.sorted.vcf",
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        bcftools view -R {input.bed} {output.vcf_gz} | bcftools sort -Ov -o {output.vcf}
        """
