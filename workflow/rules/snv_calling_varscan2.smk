rule samtools_mpileup:
    input:
        normal=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor=lambda w: f"bam/{w.run}/{w.sample}.bam",
        refg=config["paths"]["refs"]["genome_human"],
    output:
        "mpileup/{run}/{sample}.mpileup",
    conda:
        "../envs/varscan.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        samtools mpileup \
        -f {input.refg} \
        -q 1 \
        -B \
        {input.normal} {input.tumor} > {output}
        """


rule varscan_somatic:
    input:
        mpileup="mpileup/{run}/{sample}.mpileup",
    output:
        snp=temp("vcf/{run}/{sample}.varscan.snp.vcf"),
        indel=temp("vcf/{run}/{sample}.varscan.indel.vcf"),
    params:
        min_coverage=8,
        min_var_freq=0.1,
        p_value=0.05,
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        varscan somatic \
        {input.mpileup} \
        --output-snp {output.snp} \
        --output-indel {output.indel} \
        --min-coverage {params.min_coverage} \
        --min-var-freq {params.min_var_freq} \
        --p-value {params.p_value} \
        --somatic-p-value {params.p_value} \
        --strand-filter 1 \
        --output-vcf 1 \
        --mpileup 1
        """


rule rename_vcf_samples:
    input:
        "vcf/{run}/{sample}.varscan.{type}.vcf",
    output:
        temp("vcf/{run}/{sample}.varscan.{type}.renamed.vcf"),
    wildcard_constraints:
        type="snp|indel",
    params:
        normal=lambda w: runs_dict[w.run]["normal"],
    script:
        "../scripts/rename_vcf_sample_columns.py"


rule compress_vcf:
    input:
        "vcf/{run}/{sample}.varscan.{type}.renamed.vcf",
    output:
        vcf=temp("vcf/{run}/{sample}.varscan.{type}.vcf.gz"),
        tbi=temp("vcf/{run}/{sample}.varscan.{type}.vcf.gz.tbi"),
    wildcard_constraints:
        type="snp|indel",
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        bgzip -c {input} > {output.vcf} &&
        tabix -p vcf {output.vcf}
        """


rule concat_varscan_vcfs:
    input:
        vcf_snv="vcf/{run}/{sample}.varscan.snp.vcf.gz",
        tbi_snv="vcf/{run}/{sample}.varscan.snp.vcf.gz.tbi",
        vcf_indel="vcf/{run}/{sample}.varscan.indel.vcf.gz",
        tbi_indel="vcf/{run}/{sample}.varscan.indel.vcf.gz.tbi",
    output:
        vcf="vcf/{run}/{sample}.varscan.vcf.gz",
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        bcftools concat -a --allow-overlaps {input.vcf_snv} {input.vcf_indel} | \
        bcftools sort -Oz -o {output.vcf} &&
        bcftools index -t {output.vcf}
        """


rule merge_varscan_vcfs:
    input:
        vcfs=lambda w: [
            f"vcf/{w.run}/{sample}.varscan.vcf.gz"
            for sample in runs_dict[w.run]["tumors"]
        ],
    output:
        vcf="vcf/{run}/{run}.varscan.vcf",
    conda:
        "../envs/varscan.yaml"
    shell:
        """
        bcftools merge --merge none \
            --force-samples \
            {input.vcfs} | \
        bcftools sort -Ov -o {output.vcf}
        """
