rule samtools_mpileup:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        refg=config["paths"]["refs"]["genome_human"],
    output:
        mpileup="mpileup/{run}/{sample}.mpileup",
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
        snp=temp("vcf/{run}/{sample}/{sample}.varscan.snp.vcf"),
        indel=temp("vcf/{run}/{sample}/{sample}.varscan.indel.vcf"),
    params:
        min_coverage=config["params"]["varscan_min_coverage"],
        min_vaf=config["params"]["varscan_min_VAF"],
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
        --min-var-freq {params.min_vaf} \
        --p-value {params.p_value} \
        --somatic-p-value {params.somatic_p_value} \
        --strand-filter 1 \
        --output-vcf 1
        """

rule varscan_filter:
    input:
        vcf_snv="vcf/{run}/{sample}/{sample}.varscan.snp.vcf",
        vcf_indel="vcf/{run}/{sample}/{sample}.varscan.indel.vcf",
    output:
        vcf=temp("vcf/{run}/{sample}/{sample}.varscan.filtered_snp.vcf"),
    params:
        min_cov=config["params"]["varscan_min_coverage"],
        min_reads=config["params"]["varscan_filter_min_reads"],
        min_strands=config["params"]["varscan_filter_min_strands"],
        min_qual=config["params"]["varscan_filter_min_qual"],
        min_vaf=config["params"]["varscan_filter_min_vaf"],
        p_value=config["params"]["varscan_filter_p_value"],
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        varscan somaticFilter \
        {input.vcf_snv} \
        --output-file {output.vcf} \
        --min-coverage {params.min_cov} \
        --min-reads2 {params.min_reads} \
        --min-strands2 {params.min_strands} \
        --min-avg-qual {params.min_qual} \
        --min-var-freq {params.min_vaf} \
        --p-value {params.p_value} \
        --indel-file {input.vcf_indel}
        """


rule compress_vcf:
    input:
        vcf="vcf/{run}/{sample}/{sample}.varscan.{type}.vcf",
    output:
        vcf=temp("vcf/{run}/{sample}/{sample}.varscan.{type}.vcf.gz"),
        tbi=temp("vcf/{run}/{sample}/{sample}.varscan.{type}.vcf.gz.tbi"),
    wildcard_constraints:
        type="filtered_snp|indel",
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf} &&
        tabix -p vcf {output.vcf}
        """


rule concat_varscan_vcfs:
    input:
        vcf_snv="vcf/{run}/{sample}/{sample}.varscan.filtered_snp.vcf.gz",
        tbi_snv="vcf/{run}/{sample}/{sample}.varscan.filtered_snp.vcf.gz.tbi",
        vcf_indel="vcf/{run}/{sample}/{sample}.varscan.indel.vcf.gz",
        tbi_indel="vcf/{run}/{sample}/{sample}.varscan.indel.vcf.gz.tbi",
    output:
        vcf="vcf/{run}/{sample}/{sample}.varscan.vcf",
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        bcftools concat -a --allow-overlaps {input.vcf_snv} {input.vcf_indel} | \
        bcftools sort -Ov -o {output.vcf}
        """
