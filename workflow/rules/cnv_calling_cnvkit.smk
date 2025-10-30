def get_cnvkit_reference(wildcards):
    probe_version = probe_dict[wildcards.run][wildcards.sample]
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    chr_sex = sample_row["Chr_sex"].iloc[0]

    ref_key = "cnvkit_ref_m" if chr_sex == "XY" else "cnvkit_ref_f"
    return config["probe_configs"][probe_version][ref_key]


def get_normal_sample_name(wildcards):
    return runs_dict[wildcards.run]["normal"]


def get_sample_sex(wildcards):
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    chr_sex = sample_row["Chr_sex"].iloc[0]
    return "male" if chr_sex == "XY" else "female"


def is_male_reference(wildcards):
    ref_path = get_cnvkit_reference(wildcards)
    return "--male-reference" if "_male.cnn" in ref_path or "ref_m" in ref_path else ""


rule cnvkit_coverage_and_fix:
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        targets=lambda w: config["probe_configs"][get_probe_version(w)][
            "cnvkit_targets"
        ],
        antitargets=lambda w: config["probe_configs"][get_probe_version(w)][
            "cnvkit_antitargets"
        ],
        reference=get_cnvkit_reference,
    output:
        cnr="work/cnvkit/{run}/{sample}/{sample}.cnr",
    params:
        prefix="work/cnvkit/{run}/{sample}/{sample}",
    threads: config["resources"]["threads"]
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_coverage_fix_{run}_{sample}.log",
    shell:
        """
        cnvkit.py coverage {input.bam} {input.targets} \
            -o {params.prefix}.targetcoverage.cnn \
            -p {threads} > {log} 2>&1
        
        cnvkit.py coverage {input.bam} {input.antitargets} \
            -o {params.prefix}.antitargetcoverage.cnn \
            -p {threads} >> {log} 2>&1
        
        cnvkit.py fix {params.prefix}.targetcoverage.cnn \
            {params.prefix}.antitargetcoverage.cnn \
            {input.reference} \
            -o {output.cnr} >> {log} 2>&1
        
        rm {params.prefix}.targetcoverage.cnn {params.prefix}.antitargetcoverage.cnn
        """


rule gatk_haplotypecaller:
    input:
        bam="results/{run}/{sample}/bam/{sample}.bam",
        bai="results/{run}/{sample}/bam/{sample}.bai",
        refg=config["refs"]["genome_human"],
        regions=lambda w: f"work/refs/regions/{probe_dict[w.run][w.sample]}/regions.bed",
    output:
        vcf="work/haplotypecaller/{run}/{sample}/{sample}.germline.vcf",
        idx="work/haplotypecaller/{run}/{sample}/{sample}.germline.vcf.idx",
    params:
        ref_path=config["refs"]["path"],
    threads: config["resources"]["threads"]
    resources:
        java_max_gb=config["resources"]["java_max_gb"],
        java_min_gb=config["resources"]["java_min_gb"],
    log:
        "work/logs/HaplotypeCaller_{run}_{sample}.log",
    container:
        "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        gatk \
        --java-options "-Xms{resources.java_min_gb}G -Xmx{resources.java_max_gb}G" \
        HaplotypeCaller \
        -R {input.refg} \
        -I {input.bam} \
        -L {input.regions} \
        -O {output.vcf} \
        --native-pair-hmm-threads {threads} \
        > {log} 2>&1
        """


rule cnvkit_merge_germline_and_filter_hetsnp:
    input:
        normal=lambda w: f"work/haplotypecaller/{w.run}/{runs_dict[w.run]['normal']}/{runs_dict[w.run]['normal']}.germline.vcf",
        tumor="work/haplotypecaller/{run}/{sample}/{sample}.germline.vcf",
    output:
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    params:
        normal=lambda w: runs_dict[w.run]["normal"],
        tumor=lambda w: w.sample,
    conda:
        "../envs/bcftools.yaml"
    log:
        "work/logs/merge_filter_germline_{run}_{sample}.log",
    shell:
        """
        bgzip -c {input.normal} > {input.normal}.gz 2> {log}
        tabix -p vcf {input.normal}.gz 2>> {log}
        
        bgzip -c {input.tumor} > {input.tumor}.gz 2>> {log}
        tabix -p vcf {input.tumor}.gz 2>> {log}
        
        bcftools merge -m none {input.normal}.gz {input.tumor}.gz -Ou 2>> {log} | \
        bcftools view -s {params.normal},{params.tumor} -Ou | \
        bcftools filter -i "
                TYPE='snp' &&
                N_ALT==1 &&
                STRLEN(REF)==1 &&
                STRLEN(ALT)==1 &&
                FORMAT/DP[0]>10 &&
                FORMAT/AD[0:1]>3 &&
                GT[0]=='0/1' &&
                FORMAT/DP[1]>10" -Ov -o {output.vcf} 2>> {log}
        """


rule filter_chromosomes:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.cnr",
    output:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
    params:
        chromosomes=[str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"],
    run:
        def check_chromosome(chrom, allowed_chroms):
            """Check if chromosome (with or without chr prefix) is allowed"""
            return chrom.replace("chr", "") in allowed_chroms


        def filter_file(file_in, file_out, allowed_chroms):
            """Filter file keeping only allowed chromosomes"""
            with open(file_in) as in_handle, open(file_out, "w") as out_handle:
                header = True
                for line in in_handle:
                    if header:
                        out_handle.write(line)
                        header = False
                    else:
                        chrom = line.split("\t")[0]
                        if check_chromosome(chrom, allowed_chroms):
                            out_handle.write(line)


        filter_file(input.cnr, output.cnr, params.chromosomes)


rule cnvkit_segment:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.cns",
    params:
        normal=get_normal_sample_name,
        tumor=lambda w: w.sample,
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_segment_{run}_{sample}.log",
    shell:
        """
        cnvkit.py segment {input.cnr} \
            -o {output.cns} \
            -m cbs \
            -v {input.vcf} \
            -n {params.normal} \
            -i {params.tumor} \
            > {log} 2>&1
        """


rule cnvkit_segmetrics:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns="work/cnvkit/{run}/{sample}/{sample}.cns",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.segmetrics.cns",
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_segmetrics_{run}_{sample}.log",
    shell:
        """
        cnvkit.py segmetrics {input.cnr} \
            -s {input.cns} \
            --ci --alpha 0.5 --smooth-bootstrap \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_filter_ci:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.segmetrics.cns",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.ci_filtered.cns",
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_filter_ci_{run}_{sample}.log",
    shell:
        """
        cnvkit.py call {input.cns} \
            --method none \
            --filter ci \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_add_ttest:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns=lambda w: (
            f"work/cnvkit/{w.run}/{w.sample}/{w.sample}.ci_filtered.cns"
            if config["params"]["cnvkit"]["filter_ci"]
            else f"work/cnvkit/{w.run}/{w.sample}/{w.sample}.segmetrics.cns"
        ),
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.ttest.cns",
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_ttest_{run}_{sample}.log",
    shell:
        """
        cnvkit.py segmetrics {input.cnr} \
            -s {input.cns} \
            --t-test \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_call:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.ttest.cns",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.call.cns",
    params:
        normal=get_normal_sample_name,
        tumor=lambda w: w.sample,
        sample_sex=get_sample_sex,
        male_reference=is_male_reference,
        purity=get_purity,
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_call_{run}_{sample}.log",
    shell:
        """
        cnvkit.py call {input.cns} \
            -v {input.vcf} \
            -n {params.normal} \
            -i {params.tumor} \
            --sample-sex {params.sample_sex} \
            {params.male_reference} \
            --purity {params.purity} \
            -m clonal \
            --center median \
            --drop-low-coverage \
            -o {output.cns} \
            > {log} 2>&1
        """


rule cnvkit_plots:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns="work/cnvkit/{run}/{sample}/{sample}.call.cns",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    output:
        scatter="results/{run}/{sample}/{sample}.scatter.png",
        diagram="results/{run}/{sample}/{sample}.diagram.pdf",
    params:
        normal=get_normal_sample_name,
        tumor=lambda w: w.sample,
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_plots_{run}_{sample}.log",
    shell:
        """
        cnvkit.py scatter {input.cnr} \
            -s {input.cns} \
            -v {input.vcf} \
            -n {params.normal} \
            -i {params.tumor} \
            --title {params.tumor} \
            --y-max 4 --y-min -4 \
            --fig-size 10 5 \
            -o {output.scatter} \
            > {log} 2>&1

        cnvkit.py diagram {input.cnr} \
            -s {input.cns} \
            --no-gene-labels \
            --title {params.tumor} \
            -o {output.diagram} \
            > {log} 2>&1
        """
