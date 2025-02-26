rule index_bed:
    input:
        bed=config["paths"]["refs"]["regions_bedfile"],
    output:
        bed="refs/regions/regions.bed",
        gz="refs/regions/regions.bed.gz",
        tbi="refs/regions/regions.bed.gz.tbi",
    shell:
        r"""
        grep -v '^browser\|^track' {input.bed} > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip > {output.gz}
        tabix -p bed {output.gz}
        """

rule get_chromosomes:
    input:
        bed="refs/regions/regions.bed",
    output:
        chrom="refs/chrom.txt"
    shell:
        """
        cut -f1 {input.bed} | uniq | sort -k1,1V -k2,2n > {output.chrom}
        """
