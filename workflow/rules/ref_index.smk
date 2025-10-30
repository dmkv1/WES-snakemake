rule index_bed:
    input:
        bed=lambda wildcards: config["probe_configs"][wildcards.probe_version][
            "regions_bedfile"
        ],
    output:
        bed="work/refs/regions/{probe_version}/regions.bed",
        gz="work/refs/regions/{probe_version}/regions.bed.gz",
        tbi="work/refs/regions/{probe_version}/regions.bed.gz.tbi",
    conda:
        "../envs/bwamem.yaml"
    shell:
        r"""
        grep -v '^browser\|^track' {input.bed} > {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip > {output.gz}
        tabix -p bed {output.gz}
        """


rule get_chromosomes:
    input:
        bed="work/refs/regions/{probe_version}/regions.bed",
    output:
        chrom="work/refs/regions/{probe_version}/chrom.txt",
    conda:
        "../envs/bwamem.yaml"
    shell:
        """
        cut -f1 {input.bed} | uniq | sort -k1,1V -k2,2n > {output.chrom}
        """
