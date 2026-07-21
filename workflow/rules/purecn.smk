def get_sample_sex_purecn(wildcards):
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    gender = sample_row["gender"].iloc[0]
    return "M" if gender == "m" else "F"


rule cnvkit_export_seg:
    input:
        cns="work/cnvkit/{run}/{sample}/{sample}.cns",
    output:
        seg="work/purecn/{run}/{sample}/{sample}.seg",
    container:
        config["containers"]["cnvkit"]
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
        normaldb=get_purecn_normaldb,
        mapping_bias=get_purecn_mapping_bias,
        snp_blacklist=config["refs"]["purecn"]["snp_blacklist"],
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
        PURECN_SCRIPT=$(Rscript -e 'cat(system.file("extdata", "PureCN.R", package = "PureCN"))')
        Rscript "$PURECN_SCRIPT" \
            --out {params.out_dir} \
            --sampleid {params.sample_id} \
            --tumor {input.cnr} \
            --sex {params.sex} \
            --seg-file {input.seg} \
            --vcf {input.vcf} \
            --genome {params.genome} \
            --normaldb {input.normaldb} \
            --mapping-bias-file {input.mapping_bias} \
            --snp-blacklist {input.snp_blacklist} \
            --fun-segmentation Hclust \
            --min-base-quality 20 \
            --post-optimize \
            --cores {threads} \
            --force --seed 42 \
            > {log} 2>&1
        """


rule resolve_purity_source:
    """Single decision point for tumor purity/ploidy: PureCN (paired tumors,
    when not Flagged and config-enabled) or the samplesheet constant
    (tumor-only, PureCN disabled/failed-to-be-confident, or Flagged).
    cnvkit_call and combine_results both read this sidecar rather than each
    independently deciding, so they can never disagree on which value was
    actually used."""
    input:
        purecn_csv=lambda w: (
            f"work/purecn/{w.run}/{w.sample}/{w.sample}.csv"
            if is_purecn_eligible(w) else []
        ),
    output:
        purity_csv="work/purity/{run}/{sample}/{sample}.purity.csv",
    params:
        eligible=is_purecn_eligible,
        toggle=config["params"]["cnv"]["purity_source"],
        sheet_purity=get_purity,
    log:
        "work/logs/resolve_purity_source_{run}_{sample}.log",
    run:
        import csv
        import os

        purecn_purity = purecn_ploidy = purecn_flagged = ""
        source = "samplesheet"

        if params.eligible:
            with open(input.purecn_csv) as fh:
                row = next(csv.DictReader(fh))
            purecn_purity = row["Purity"]
            purecn_ploidy = row["Ploidy"]
            purecn_flagged = row["Flagged"]

            if params.toggle == "purecn":
                is_flagged = str(purecn_flagged).strip().upper() == "TRUE"
                has_valid_purity = purecn_purity not in ("", "NA", None)
                if not is_flagged and has_valid_purity:
                    source = "purecn"
                else:
                    source = "samplesheet_fallback_flagged"

        purity = purecn_purity if source == "purecn" else params.sheet_purity
        ploidy = purecn_ploidy if source == "purecn" else ""

        os.makedirs(os.path.dirname(output.purity_csv), exist_ok=True)
        with open(output.purity_csv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "run",
                    "sample",
                    "purity",
                    "ploidy",
                    "source",
                    "purecn_purity",
                    "purecn_ploidy",
                    "purecn_flagged",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "run": wildcards.run,
                    "sample": wildcards.sample,
                    "purity": purity,
                    "ploidy": ploidy,
                    "source": source,
                    "purecn_purity": purecn_purity,
                    "purecn_ploidy": purecn_ploidy,
                    "purecn_flagged": purecn_flagged,
                }
            )
