def get_sample_sex_purecn(wildcards):
    sample_row = samples[
        (samples["ID"] == wildcards.run) & (samples["sample"] == wildcards.sample)
    ]
    gender = sample_row["gender"].iloc[0]
    return "M" if gender == "m" else "F"


rule purecn_tumor_coverage:
    # Reformats the tumor .cnr into PureCN's GATK3-style coverage
    # (Target/total_coverage/on_target) so PureCN can run its OWN segmentation.
    # PureCN misreads CNVkit's .cnr `depth` column as per-interval read counts
    # and drops every interval below its 100-read floor; the depth*width scale
    # here matches how the NormalDB coverage was built (WES-PON-smk
    # purecn_coverage), so tumor and normals are comparable. This replaces the
    # old cnvkit_export_seg handoff: feeding PureCN CNVkit's CBS segmentation
    # (--seg-file) inherited its over-segmentation (500-700 segments > PureCN's
    # max.segments=300) and collapsed purity to a flat NON-ABERRANT optimum.
    # Letting PureCN segment (native CBS, noise-calibrated) recovered PDX purity
    # from 0.18-0.32 to 0.86-0.93. CNVkit's own calling arm is unaffected.
    input:
        cnr="work/cnvkit/{run}/{sample}/{sample}.filtered.cnr",
    output:
        cov="work/purecn/{run}/{sample}/{sample}.purecn_cov.txt",
    run:
        cov = pd.read_csv(input.cnr, sep="\t")
        cov["on_target"] = cov["gene"] != "Antitarget"
        cov["Target"] = (
            cov["chromosome"]
            + ":"
            + (cov["start"] + 1).astype(str)
            + "-"
            + cov["end"].astype(str)
        )
        cov["total_coverage"] = cov["depth"] * (cov["end"] - cov["start"])
        cov[["Target", "total_coverage", "on_target"]].to_csv(
            output.cov, sep="\t", index=False
        )


rule purecn_run:
    input:
        vcf="work/cnvkit/{run}/{sample}/{sample}.hetsnp.vcf",
        cov="work/purecn/{run}/{sample}/{sample}.purecn_cov.txt",
        normaldb=get_purecn_normaldb,
        mapping_bias=get_purecn_mapping_bias,
        snp_blacklist=config["refs"]["purecn"]["snp_blacklist"],
    output:
        csv="work/purecn/{run}/{sample}/{sample}.csv",
        rds="work/purecn/{run}/{sample}/{sample}.rds",
        pdf="work/purecn/{run}/{sample}/{sample}.pdf",
        # No _genes.csv: gene-level calls need gene-annotated intervals, which
        # PureCN previously got from the CNVkit .cnr `gene` column. The native
        # coverage (Target/total_coverage/on_target) carries no genes, so PureCN
        # skips gene-level calls ("--intervals does not contain gene symbols").
        # This artifact was unused anyway — the report's gene-level CNVs come
        # from CNVkit's .cns, and nothing reads PureCN _genes.csv. To restore it,
        # pass a PureCN interval file (--intervals) built with gene symbols.
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
            --tumor {input.cov} \
            --sex {params.sex} \
            --vcf {input.vcf} \
            --genome {params.genome} \
            --normaldb {input.normaldb} \
            --mapping-bias-file {input.mapping_bias} \
            --snp-blacklist {input.snp_blacklist} \
            --fun-segmentation CBS \
            --min-base-quality 20 \
            --post-optimize \
            --cores {threads} \
            --force --seed 42 \
            > {log} 2>&1
        """


rule resolve_purity_source:
    """Single decision point for the purity/ploidy that cnvkit_call and
    combine_results consume, in priority order:
      1. 'known'        — samplesheet tumor_fraction (orthogonal ground truth:
                          sorting/cytometry, PDX=1). Overrides PureCN.
      2. 'purecn'       — PureCN's estimate, IF it ran, did not Fail, and its
                          ONLY flag (if any) is POOR GOF. No numeric GOF cutoff;
                          any other flag (NON-ABERRANT, LOW PURITY, NOISY
                          SEGMENTATION, EXCESSIVE LOH, contamination...) rejects.
      3. 'assumed_pure' — no known value and no usable PureCN estimate: purity 1
                          (cnvkit calls unrescaled, as in the pre-PureCN default;
                          e.g. tumor-only samples).
    Both consumers read this one sidecar, so they can never disagree on the
    value actually used."""
    input:
        purecn_csv=lambda w: (
            f"work/purecn/{w.run}/{w.sample}/{w.sample}.csv"
            if is_purecn_eligible(w) else []
        ),
    output:
        purity_csv="work/purity/{run}/{sample}/{sample}.purity.csv",
    params:
        eligible=is_purecn_eligible,
        use_purecn=config["params"]["cnv"]["use_purecn_purity"],
        known=get_known_purity,
    log:
        "work/logs/resolve_purity_source_{run}_{sample}.log",
    run:
        import csv
        import os

        purecn_purity = purecn_ploidy = purecn_flagged = purecn_comment = ""
        purecn_failed = False

        if params.eligible and input.purecn_csv:
            with open(input.purecn_csv) as fh:
                row = next(csv.DictReader(fh))
            purecn_purity = row["Purity"]
            purecn_ploidy = row["Ploidy"]
            purecn_flagged = row["Flagged"]
            purecn_comment = row.get("Comment", "")
            purecn_failed = str(row.get("Failed", "")).strip().upper() == "TRUE"

        def purecn_usable():
            # Accept PureCN only if enabled, it ran, didn't fail, has a purity,
            # and its only flag reason (if flagged) is POOR GOF.
            if not (params.use_purecn and params.eligible and input.purecn_csv):
                return False
            if purecn_failed or purecn_purity in ("", "NA", None):
                return False
            if str(purecn_flagged).strip().upper() != "TRUE":
                return True  # not flagged
            reasons = [r.strip() for r in str(purecn_comment).split(";") if r.strip()]
            return bool(reasons) and all(
                r.upper().startswith("POOR GOF") for r in reasons
            )

        # Integer baseline ploidy for cnvkit --ploidy (which is type=int) and the
        # QC table, resolved independently of the purity source. PureCN reports
        # only a continuous Ploidy, so round it half-up to the nearest integer
        # (min 1); fall back to diploid when PureCN did not run, failed, or
        # produced no usable value.
        def round_ploidy(val):
            try:
                p = float(val)
            except (TypeError, ValueError):
                return None
            return max(1, int(p + 0.5)) if p > 0 else None

        ploidy_int = None if purecn_failed else round_ploidy(purecn_ploidy)
        ploidy = str(ploidy_int) if ploidy_int is not None else "2"

        known = params.known  # float in (0,1] or None
        if known is not None:
            purity, source = str(known), "known"
        elif purecn_usable():
            purity, source = purecn_purity, "purecn"
        else:
            purity, source = "1", "assumed_pure"

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
                    "tumor_fraction",
                    "purecn_purity",
                    "purecn_ploidy",
                    "purecn_flagged",
                    "purecn_comment",
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
                    "tumor_fraction": "" if known is None else known,
                    "purecn_purity": purecn_purity,
                    "purecn_ploidy": purecn_ploidy,
                    "purecn_flagged": purecn_flagged,
                    "purecn_comment": purecn_comment,
                }
            )
