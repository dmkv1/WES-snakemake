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
    benchmark:
        "work/benchmarks/purecn_tumor_coverage/{run}_{sample}.tsv",
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
    resources:
        # R holds the interval coverage and the segmentation for the sample in
        # memory across the whole fit, and peaks near 10 GB on a WES tumour.
        mem_mb=16384,
    benchmark:
        "work/benchmarks/purecn_run/{run}_{sample}.tsv",
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
    combine_results consume.

    Purity, in priority order:
      1. 'known'        — samplesheet tumor_fraction (orthogonal ground truth:
                          sorting/cytometry, PDX=1). Overrides PureCN.
      2. 'purecn'       — PureCN's estimate, IF it ran, did not Fail, and
                          produced a numeric Purity. Flag content (LOW PURITY,
                          NON-ABERRANT, POOR GOF, NOISY SEGMENTATION,
                          EXCESSIVE LOH, ...) no longer rejects the estimate:
                          a flagged fit is still the best purity we have, and
                          forcing purity 1 instead only throws away signal —
                          see purity_confidence below for the QC signal that
                          replaces the old reject. Before 2.3.0 any flag other
                          than POOR GOF rejected the estimate outright, which
                          silently zeroed real CNAs in ~1 in 4 samples cohort-
                          wide (WES-MCL-II notebook, 2026-09-09).
      3. 'assumed_pure' — no known value and no usable PureCN estimate: purity 1
                          (cnvkit calls unrescaled, as in the pre-PureCN default;
                          e.g. tumor-only samples, or PureCN failed/didn't run).

    Ploidy (tracked separately as ploidy_source, resolved independently of the
    purity source above — see the 2026-08-24 notebook entries for why a
    sample can be e.g. purity 'assumed_pure' + ploidy 'purecn' at once), in
    priority order:
      1. 'known'   — samplesheet known_ploidy (orthogonal ground truth:
                     karyotype, flow DNA index). Overrides PureCN. Added after
                     P038 PDCL: karyotype ~42 chr (near-diploid, one wrongly
                     tetraploid subclone) but PureCN fit ploidy 3 anyway (POOR
                     GOF alone doesn't reject a PureCN estimate) — corroborated
                     wrong by an independent BCL2 FISH-vs-WES discordance.
      2. 'purecn'  — PureCN's raw Ploidy, IF it ran and did not Fail (no flag
                     check — flags gate purity usability, not ploidy).
      3. 'default' — no known value and no usable PureCN estimate: diploid (2).

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
        known_ploidy=get_known_ploidy,
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

        # Integer baseline ploidy for cnvkit --ploidy (which is type=int) and the
        # QC table. PureCN reports only a continuous Ploidy, so round it half-up
        # to the nearest integer (min 1); fall back to diploid when PureCN did
        # not run, failed, or produced no usable value. No flag check here (see
        # docstring) — flags gate purity usability, not this rounding step.
        def round_ploidy(val):
            try:
                p = float(val)
            except (TypeError, ValueError):
                return None
            return max(1, int(p + 0.5)) if p > 0 else None

        known = params.known  # float in (0,1] or None
        purecn_available = bool(params.use_purecn and params.eligible and input.purecn_csv)
        purity, source, purity_confidence = resolve_purity(
            known, purecn_purity, purecn_failed, purecn_available
        )

        known_ploidy = params.known_ploidy  # int or None
        purecn_ploidy_int = None if purecn_failed else round_ploidy(purecn_ploidy)
        if known_ploidy is not None:
            ploidy, ploidy_source = str(known_ploidy), "known"
        elif purecn_ploidy_int is not None:
            ploidy, ploidy_source = str(purecn_ploidy_int), "purecn"
        else:
            ploidy, ploidy_source = "2", "default"

        os.makedirs(os.path.dirname(output.purity_csv), exist_ok=True)
        with open(output.purity_csv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "run",
                    "sample",
                    "purity",
                    "purity_confidence",
                    "ploidy",
                    "source",
                    "ploidy_source",
                    "tumor_fraction",
                    "known_ploidy",
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
                    "purity_confidence": purity_confidence,
                    "ploidy": ploidy,
                    "source": source,
                    "ploidy_source": ploidy_source,
                    "tumor_fraction": "" if known is None else known,
                    "known_ploidy": "" if known_ploidy is None else known_ploidy,
                    "purecn_purity": purecn_purity,
                    "purecn_ploidy": purecn_ploidy,
                    "purecn_flagged": purecn_flagged,
                    "purecn_comment": purecn_comment,
                }
            )
