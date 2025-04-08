rule index_bed:
    input:
        refdir=config["paths"]["refs"]["path"],
        bed=config["paths"]["refs"]["regions_bedfile"],
    output:
        bed="reference/regions/regions.bed",
        gz="reference/regions/regions.bed.gz",
        tbi="reference/regions/regions.bed.gz.tbi",
    singularity:
        "docker://ubuntu:24.04"
    shell:
        r"""
        grep -v '^browser\|^track' {input.refdir}/{input.bed} > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip > {output.gz}
        tabix -p bed {output.gz}
        """


rule get_chromosomes:
    input:
        bed="refs/regions/regions.bed",
    output:
        chrom="refs/chrom.txt",
    singularity:
        "docker://ubuntu:20.04"
    shell:
        """
        cut -f1 {input.bed} | uniq | sort -k1,1V -k2,2n > {output.chrom}
        """
