import os


rule run_somaticseq:
    input:
        normal_bam=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor_bam=lambda w: f"bam/{w.run}/{w.sample}.bam",
        vcf_mutect=lambda w: f"vcf/{w.run}/{w.sample}/mutect/{w.sample}.mutect2.sorted.vcf",
        vcf_varscan_snv=lambda w: f"vcf/{w.run}/{w.sample}/varscan/{w.sample}.varscan.sorted.vcf",
        vcf_varscan_indel=lambda w: f"vcf/{w.run}/{w.sample}/varscan/{w.sample}.varscan.indel.vcf",
        vcf_strelka_snv=lambda w: f"vcf/{w.run}/{w.sample}/strelka/results/variants/{w.sample}.strelka.sorted.vcf",
        vcf_strelka_indel=lambda w: f"vcf/{w.run}/{w.sample}/strelka/results/variants/somatic.indels.vcf.gz",
        refg=config["paths"]["refs"]["genome_human"],
        regions="refs/regions/regions.bed",
    output:
        snv_vcf=temp("vcf/{run}/{sample}/somaticseq/Consensus.sSNV.vcf"),
        indel_vcf=temp("vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.vcf"),
        table_snv="vcf/{run}/{sample}/somaticseq/Ensemble.sSNV.tsv",
        table_indel="vcf/{run}/{sample}/somaticseq/Ensemble.sINDEL.tsv",
    params:
        outdir=lambda w, output: os.path.dirname(output.snv_vcf),
        ref_path=config["paths"]["refs"]["path"],
        image=config["tools"]["somaticseq"]["image"],
        image_ver=config["tools"]["somaticseq"]["version"],
        pass_threshold=config["params"]["somaticseq"]["pass_threshold"],
        lowqual_threshold=config["params"]["somaticseq"]["lowqual_threshold"],
        homozygous_threshold=config["params"]["somaticseq"]["homozygous_threshold"],
        heterozygous_threshold=config["params"]["somaticseq"]["heterozygous_threshold"],
        minimum_mapping_quality=config["params"]["somaticseq"][
            "minimum_mapping_quality"
        ],
        minimum_base_quality=config["params"]["somaticseq"]["minimum_base_quality"],
        minimum_num_callers=config["params"]["somaticseq"]["minimum_num_callers"],
    resources:
        threads=config["resources"]["threads"],
    log:
        "logs/{run}/{sample}/somaticseq_parallel.log",
    shell:
        """
        docker run --rm \
        -v {params.ref_path}:{params.ref_path} \
        -v $PWD:$PWD -w $PWD \
        --user $(id -u):$(id -g) \
        {params.image}:{params.image_ver} \
        somaticseq_parallel.py \
            -outdir {params.outdir} \
            -ref {input.refg} \
            --pass-threshold {params.pass_threshold} \
            --lowqual-threshold {params.lowqual_threshold} \
            --homozygous-threshold {params.homozygous_threshold} \
            --heterozygous-threshold {params.heterozygous_threshold} \
            --minimum-mapping-quality {params.minimum_mapping_quality} \
            --minimum-base-quality {params.minimum_base_quality} \
            --minimum-num-callers {params.minimum_num_callers} \
            --inclusion-region {input.regions} \
            --threads {resources.threads} \
            paired \
            --tumor-bam-file {input.tumor_bam} \
            --normal-bam-file {input.normal_bam} \
            --mutect2-vcf {input.vcf_mutect} \
            --varscan-snv {input.vcf_varscan_snv} \
            --varscan-indel {input.vcf_varscan_indel} \
            --strelka-snv {input.vcf_strelka_snv} \
            --strelka-indel {input.vcf_strelka_indel} \
            > {log} 2>&1
        """


rule mutect_extract_contigs:
    input:
        vcf_mutect=lambda w: f"vcf/{w.run}/{w.sample}/mutect/{w.sample}.mutect2.sorted.vcf",
        chrom="refs/chrom.txt",
    output:
        contigs=temp("vcf/{run}/{sample}/mutect/{sample}.filtered.contigs.txt"),
    shell:
        """
        while read chr; do grep "^##contig=<ID=${{chr}}," {input.vcf_mutect} >> {output.contigs}; done < {input.chrom}
        """


rule somaticseq_add_contigs:
    input:
        snv_vcf="vcf/{run}/{sample}/somaticseq/Consensus.sSNV.vcf",
        indel_vcf="vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.vcf",
        contigs="vcf/{run}/{sample}/mutect/{sample}.filtered.contigs.txt",
    output:
        snv_vcf=temp("vcf/{run}/{sample}/somaticseq/Consensus.sSNV.contigs.vcf"),
        indel_vcf=temp("vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.contigs.vcf"),
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bcftools reheader --header <(cat <(grep "^##" {input.snv_vcf} | grep -v "^#CHROM") {input.contigs} <(grep "^#CHROM" {input.snv_vcf})) {input.snv_vcf} > {output.snv_vcf}
        bcftools reheader --header <(cat <(grep "^##" {input.indel_vcf} | grep -v "^#CHROM") {input.contigs} <(grep "^#CHROM" {input.indel_vcf})) {input.indel_vcf} > {output.indel_vcf}
        """


rule somaticseq_compress:
    input:
        snv_vcf="vcf/{run}/{sample}/somaticseq/Consensus.sSNV.contigs.vcf",
        indel_vcf="vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.contigs.vcf",
    output:
        snv_vcf_gz=temp("vcf/{run}/{sample}/somaticseq/Consensus.sSNV.contigs.vcf.gz"),
        indel_vcf_gz=temp(
            "vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.contigs.vcf.gz"
        ),
        snv_tbi=temp("vcf/{run}/{sample}/somaticseq/Consensus.sSNV.contigs.vcf.gz.tbi"),
        indel_tbi=temp(
            "vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.contigs.vcf.gz.tbi"
        ),
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bcftools sort -Ov {input.snv_vcf} | bgzip > {output.snv_vcf_gz}
        tabix -p vcf {output.snv_vcf_gz}
        bcftools sort -Ov {input.indel_vcf} | bgzip > {output.indel_vcf_gz}
        tabix -p vcf {output.indel_vcf_gz}
        """


rule somaticseq_concat:
    input:
        snv_vcf_gz="vcf/{run}/{sample}/somaticseq/Consensus.sSNV.contigs.vcf.gz",
        indel_vcf_gz="vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.contigs.vcf.gz",
        snv_tbi="vcf/{run}/{sample}/somaticseq/Consensus.sSNV.contigs.vcf.gz.tbi",
        indel_tbi="vcf/{run}/{sample}/somaticseq/Consensus.sINDEL.contigs.vcf.gz.tbi",
    output:
        merged="vcf/{run}/{sample}/somaticseq/{sample}.consensus.merged.vcf",
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bcftools concat {input.snv_vcf_gz} {input.indel_vcf_gz} -a -Ov -o {output.merged}
        """
