def get_sample_sex_purecn(wildcards):
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    chr_sex = sample_row["Chr_sex"].iloc[0]
    return "M" if chr_sex == "XY" else "F"


def get_normal_sample_name(wildcards):
    return runs_dict[wildcards.run]["normal"]


rule cnvkit_export_seg:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.cns",
    output:
        seg="work/purecn/{run}/{sample}/{sample}.seg",
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_export_seg_{run}_{sample}.log",
    shell:
        """
        cnvkit.py export seg {input.cns} \
            --enumerate-chroms \
            -o {output.seg} \
            > {log} 2>&1
        """


rule purecn_run:
    input:
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        seg="work/purecn/{run}/{sample}/{sample}.seg",
    output:
        csv="work/purecn/{run}/{sample}/{sample}.csv",
        rds="work/purecn/{run}/{sample}/{sample}.rds",
        pdf="work/purecn/{run}/{sample}/{sample}.pdf",
        genes="work/purecn/{run}/{sample}/{sample}_genes.csv",
        loh="work/purecn/{run}/{sample}/{sample}_loh.csv",
    params:
        out_dir="work/purecn/{run}/{sample}",
        sample_id=lambda w: w.sample,
        genome=config["refs"]["purecn"]["genome"],
        sex=get_sample_sex_purecn,
    threads: 4
    conda:
        "../envs/purecn.yaml"
    log:
        "work/logs/purecn_{run}_{sample}.log",
    shell:
        """
        Rscript $CONDA_PREFIX/lib/R/library/PureCN/extdata/PureCN.R \
            --out {params.out_dir} \
            --sampleid {params.sample_id} \
            --tumor {input.cnr} \
            --sex {params.sex} \
            --seg-file {input.seg} \
            --vcf {input.vcf} \
            --genome {params.genome} \
            --fun-segmentation Hclust \
            --min-base-quality 20 \
            --post-optimize \
            --cores {threads} \
            --force --seed 42 \
            > {log} 2>&1
        """


rule cnvkit_call_with_purity:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.ttest.cns",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
        purecn_csv="work/purecn/{run}/{sample}/{sample}.csv",
    output:
        cns="work/cnvkit/{run}/{sample}/{sample}.purecn_calibrated.cns",
    params:
        normal=get_normal_sample_name,
        tumor=lambda w: w.sample,
        sample_sex=lambda w: get_sample_sex(w),
        male_reference=lambda w: is_male_reference(w),
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_call_purity_{run}_{sample}.log",
    run:
        import pandas as pd
        
        # Extract purity and ploidy from PureCN
        purecn_df = pd.read_csv(input.purecn_csv)
        purity = purecn_df.loc[0, "Purity"]
        ploidy = round(purecn_df.loc[0, "Ploidy"])
        
        purity_ploidy_arg = f"--purity {purity} --ploidy {ploidy}"
        male_ref_arg = params.male_reference
        
        shell(f"""
        cnvkit.py call {{input.cns}} \
            -v {{input.vcf}} \
            -n {{params.normal}} \
            -i {{params.tumor}} \
            --sample-sex {{params.sample_sex}} \
            {male_ref_arg} \
            {purity_ploidy_arg} \
            -m clonal \
            --center median \
            --drop-low-coverage \
            -o {{output.cns}} \
            > {{log}} 2>&1
        """)


rule cnvkit_plots_calibrated:
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
        cns="work/cnvkit/{run}/{sample}/{sample}.purecn_calibrated.cns",
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
    output:
        scatter="results/{run}/{sample}/{sample}.scatter.calibrated.png",
        diagram="results/{run}/{sample}/{sample}.diagram.calibrated.pdf",
    params:
        normal=get_normal_sample_name,
        tumor=lambda w: w.sample,
    container:
        "docker://etal/cnvkit:0.9.11"
    log:
        "work/logs/cnvkit_plots_calibrated_{run}_{sample}.log",
    shell:
        """
        cnvkit.py scatter {input.cnr} \
            -s {input.cns} \
            -v {input.vcf} \
            -n {params.normal} \
            -i {params.tumor} \
            --title "{params.tumor} (PureCN calibrated)" \
            --y-max 4 --y-min -4 \
            --fig-size 10 5 \
            -o {output.scatter} \
            > {log} 2>&1

        cnvkit.py diagram {input.cnr} \
            -s {input.cns} \
            --no-gene-labels \
            --title "{params.tumor} (PureCN calibrated)" \
            -o {output.diagram} \
            > {log} 2>&1
        """