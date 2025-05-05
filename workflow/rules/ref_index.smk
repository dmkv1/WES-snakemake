rule index_bed:
    output:
        bed="reference/regions/regions.bed",
        gz="reference/regions/regions.bed.gz",
        tbi="reference/regions/regions.bed.gz.tbi",
    params:
        bed=get_ref_path(config["refs"]["regions_bedfile"], use_container=False),
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        grep -v '^browser\|^track' {params.bed} > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip > {output.gz}
        tabix -p bed {output.gz}
        """


rule get_chromosomes:
    input:
        bed="reference/regions/regions.bed",
    output:
        chrom="reference/chrom.txt",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        cut -f1 {input.bed} | uniq | sort -k1,1V -k2,2n > {output.chrom}
        """
