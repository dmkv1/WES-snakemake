rule export_bed:
    input:
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.merged.vcf",
    output:
        bed=temp("vcf/{run}/{sample}/depth/{sample}.snv_indels.bed"),
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bcftools query -f'%CHROM\t%POS\n' {input.vcf} > {output.bed}
        """


rule calculate_depth:
    input:
        bed=lambda w: f"vcf/{w.run}/{w.sample}/depth/{w.sample}.snv_indels.bed",
        normal_bam=lambda w: f"bam/{w.run}/{runs_dict[w.run]['normal']}.bam",
        tumor_bam=lambda w: f"bam/{w.run}/{w.sample}.bam",
    output:
        depth_table="vcf/{run}/{sample}/depth/{sample}.depth.tsv"
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        samtools depth -b {input.bed} {input.normal_bam} {input.tumor_bam} > {output.depth_table}
        """


rule add_depth_to_vcf:
    input:
        depth_table="vcf/{run}/{sample}/depth/{sample}.depth.tsv",
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.merged.vcf",
    output:
        depth_table_gz=temp("vcf/{run}/{sample}/depth/{sample}.depth.tsv.gz"),
        depth_table_tbi=temp("vcf/{run}/{sample}/depth/{sample}.depth.tsv.gz.tbi"),
        dp_header=temp("tmp/{run}_{sample}.dp_header.txt"),
        vcf="vcf/{run}/{sample}/somaticseq/{sample}.consensus.merged.dp.vcf",
    conda:
        "../envs/sam_vcf_tools.yaml"
    shell:
        """
        bgzip -c {input.depth_table} > {output.depth_table_gz}
        tabix -s1 -b2 -e2 {output.depth_table_gz}
        
        echo -e '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (samtools depth)">' > {output.dp_header}

        bcftools annotate -s NORMAL,TUMOR -a {output.depth_table_gz} -h {output.dp_header} -c CHROM,POS,FORMAT/DP {input.vcf} -Ov -o {output.vcf}
        """
