# WES-snakemake Pipeline — Known Issues & TODO

## Blocking / Pressing (open, affect result correctness)

25. **PureCN is not wired into the DAG** — `purecn.smk` exists (`purecn_run`, `cnvkit_call_with_purity`, `cnvkit_plots_calibrated` rules) but is never `include`d in the `Snakefile`, and no `rule all` input references its outputs: it's dead code, not currently run. Decide: wire it in (in which case #12 below matters) or remove it entirely (in which case #12 is moot and `purecn.smk`/`purecn.yaml` should go with it). Also blocks FACETS (#13) being a genuine second, orthogonal estimator — right now there's only zero estimators wired in, not one.

12. **PureCN receives pre-segmented `.cns` + `--fun-segmentation Hclust`** — contradictory; PureCN's own segmentation is better for purity estimation. Should feed `.cnr` ratios directly. *(Only relevant once PureCN is actually wired in — see #25.)*

---

## Good-to-haves (open, enhancements / non-blocking)

8. **xengsort mouse-read stats not captured** — fraction of murine reads is critical QC for PDX samples and is currently lost.

13. **No orthogonal purity/ploidy estimator (FACETS)** — purity is currently a hardcoded samplesheet constant (CTRL=0, tumors=1); PureCN is the only estimator and its quality is gated by the BAF/het-SNP track (#6, #12). Add **FACETS** (WES tumor/normal, allele-specific CN + purity/ploidy from depth-ratio + BAF) as a second, orthogonal estimator so purity can be reported as a "two methods agree" number rather than an assumption. Cheap interim cross-check available with no new tool: `2 × modal VAF of clonal PASS SNVs in CN-neutral CNVkit segments ≈ purity` (Mutect2 VCF + CNVkit `.cns`, both already produced).

19. **No MAF-format SNV output for maftools / cohort-level analysis** — current SNV outputs are per-sample VEP-annotated VCF (`{sample}.SNV.vcf`) and the merged Excel/TSV output. Both are fine for per-sample inspection but unusable as input to `maftools` (oncoplots, mutational signatures, TMB summaries, somatic interactions, lollipop plots, etc.), which is the standard downstream entry point for cohort SNV analysis.

    **STALE — roadmap below predates the Funcotator→VEP migration (#23, done).** It assumed Funcotator's native `--output-file-format MAF` emitter; that path no longer exists. Re-scope to `vcf2maf` (VEP-based) before implementing — the sample-ID wiring and cohort-aggregation steps below still apply, but step 1's Funcotator MAF call must become a `vcf2maf` call over the VEP-annotated VCF.

    **Implementation roadmap (needs re-scoping, see above):**

    1. **Per-sample MAF rule** in `snv_calling_mutect2.smk`. Two options:
       - *(preferred)* Add a sibling rule `funcotate_maf` that consumes the same PASS-filtered, AF-filtered (tumor-only) VCF as the existing `funcotate` rule and runs Funcotator a second time with `--output-file-format MAF`. Keeps the VCF path byte-identical; cheap re-run of just the annotation step.
       - *(alternative)* Add a second `output:` to the existing rule and emit both formats in one Funcotator call — slightly faster but couples the two outputs (any future change to one re-triggers the other). Prefer option 1 for separability.

    2. **Sample-ID wiring.** Funcotator MAF requires identity columns the VCF mode doesn't:
       - `--tumor-id {sample}` (always)
       - `--normal-id {normal}` — paired mode only; pull via the same helper used by `run_mutect2` (e.g. `_get_normal_sample(wildcards)`). Tumor-only: omit the flag; maftools handles single-sample MAFs fine.
       - `--sequence-dictionary` already implicit via the reference fasta `.dict`; no extra wiring.

    3. **Output path.** `results/{run}/{sample}/{sample}.SNV.maf` — sibling of the existing `.SNV.vcf`. Add to the `all` target and to the per-sample output list referenced by `combine_results.R` inputs (the R script itself does not need to read the MAF; just declare it so Snakemake builds it).

    4. **Cohort aggregation rule.** New rule `aggregate_cohort_maf` (top-level `results/cohort.maf`):
       - Inputs: all per-sample `{sample}.SNV.maf` for tumor/PDX samples (exclude CTRL — they have no SNV MAF anyway).
       - Shell: take header (`#version` line + column header) from the first file, then `tail -n +N` from each (skipping both comment and header). ~5 lines of awk or grep -v.
       - One file, no per-run subdirectory — maftools wants a single cohort MAF.

    5. **R verification snippet** (not part of pipeline, just for the verification run):
       ```r
       library(maftools)
       maf <- read.maf("results/cohort.maf")
       plotmafSummary(maf); oncoplot(maf, top = 30)
       ```

    **Caveats to document in the rule comment / README:**
    - MAF's column set is fixed; the rich extra Funcotator annotations currently preserved in the VCF/xlsx (full gnomAD_exome + gnomAD_genome AFs, ClinVar significance, COSMIC details, HGNC fields, etc.) get folded into MAF standard columns or pushed into `Other_Transcripts` — some are lost. MAF is **additional**, not a replacement: keep the VCF and the xlsx workbook unchanged.
    - Funcotator MAF emits one row per *transcript-effect*, not one per variant — `maftools::read.maf` collapses to canonical by default, but be aware when grep-ing the raw MAF.
    - Tumor-only entries will have `Matched_Norm_Sample_Barcode = NONE`; fine for maftools, but flag if downstream code assumes a normal exists.
    - PDX samples: MAF will be emitted on the xenograft-filtered (xengsort graft) BAM, same as the VCF — no special handling needed.

    **Estimated effort:** ~30 min implementation (one rule + one aggregation rule + samplesheet/target wiring) + one end-to-end verification run on an existing patient (e.g. `P005_PT`, already validated for #4/#5/#6) to confirm the MAF loads cleanly in maftools. No new containers, no new reference resources, no DAG restructuring.

21. **Add Delly as a second SV caller, merged with Manta.** WES-smk is Manta-only; BALSAMIC runs Manta + Delly merged via SVDB. **Cohort-specific payoff:** DLBCL is defined by IG translocations (BCL2/BCL6/MYC) and Manta and Delly have non-overlapping breakpoint sensitivity — a consensus/union materially improves translocation recall. This is the highest-value caller addition for this cohort (higher than a second SNV caller).
    - **Plan:** new `run_delly` rule (`delly call -x <exclude.bed>`, paired = tumor+normal with somatic pre-filter `delly filter -f somatic -s samples.tsv`; tumor-only = single-sample); PASS-filter; merge Manta + Delly with SVDB (`--no_intra --bnd_distance 5000 --overlap 0.80`, priority list) or a lighter `bcftools`/survivor merge; feed the merged VCF to the existing `annotate_sv` (AnnotSV) rule in `sv_calling_manta.smk`.
    - **Resources:** Delly exclude-regions BED for hg38 (ships with Delly). No PON required for somatic Delly. (Delly-CNV is a separate optional path — not in scope here; CNV stays CNVkit/PureCN.)
    - **Effort:** medium. New caller + merge step + AnnotSV input swap. Verify on a known-translocation case.

24. **`gatk_haplotypecaller`'s tumor-side force-call is heavier than needed for BAF extraction** — the paired-tumor call added for #6 already restricts `-L` to the normal's het sites, so it's much cheaper than the normal's full-exome call, but HaplotypeCaller still does local re-assembly + HMM realignment per active region even when restricted to single-bp intervals — overkill for what's actually needed (ref/alt depth at *known* positions, not fresh genotyping). GATK's `CollectAllelicCounts` (pileup-based, no assembly) is the tool built for exactly this — it's what GATK's own ModelSegments/allelic-CN workflow uses. The common-SNP resource already in config (`contamination_resource`, used for `GetPileupSummaries`/#5) could seed candidate sites for both normal and tumor, skipping HaplotypeCaller's assembly cost on both sides, not just the tumor's. Not a correctness bug — #6 is verified and working — this is a speed optimization, lower priority than the items above.

---

## Multi-lane (parked — branch `multi-lane`)

15. **Read group `RGID`/`RGPU` parsing fragile** — encodes only instrument:flowcell, not lane; RGIDs can collide for multi-lane samples. *(Multi-lane fix attempted, broke on RG-string parsing/bwa `-R` escaping — parked on branch `multi-lane`, not resolved here.)*

---

## Done

1. ~~**No UMI-aware trimming or deduplication for V8+UTR samples**~~ ✅ — investigated, **not a pipeline bug; no sample requires UMI removal**. The V8+UTR prep (SureSelect XT HS, G9702-90005) does generate a 10-bp i5 MBC, but CLIP/BaseSpace demux delivers FASTQs with the MBC **dropped at FASTQ generation**: i5 written as a constant 8-bp *sample* index, no I2 read, raw binned NextSeq qualities, PCR-duplicate load consistent with un-collapsed data (triangulated across three independent checks). The barcode is physically absent from every delivered FASTQ, so coordinate-based `MarkDuplicates` is the only valid option and is exactly what runs — existing results are **not compromised on the UMI axis**. `fastp` defaults (PE adapter overlap-trim, auto polyG for NextSeq, quality/N/length filter) are adequate for non-UMI WES. Recovery is possible *only* by regenerating FASTQs from BaseSpace with i5/MBC retrieval (facility request, not a code change). Detail: memory `reference_wes_umi_demux.md`.
   - *Optional, non-blocking future-proofing:* add an explicit `umi` samplesheet column (orthogonal to `probes`) gating a dormant fgbio path; set `umi=false` for all current rows (byte-identical behavior, no remap) so correctly-demuxed future runs work without touching existing results.

2. ~~**BWA-MEM `-M` flag**~~ ✅ — fixed: replaced with `-Y` in `bam_mapping_gatk.smk`.

3. ~~**BQSR trained on tumor BAM**~~ ✅ — investigated, **known caveat, not a pipeline bug; per-sample BQSR is the correct catch-all here**. Current pipeline runs `BaseRecalibrator` independently per sample (`bam_mapping_gatk.smk:144–205`), exactly as GATK4 Best Practices prescribes. The GATK3-era "apply matched-normal's recal table to the tumor" advice is (a) **negligible at this mutation load** — coding TMB ≈ 4.5–12 mut/Mb (150–400 CDS SNVs / ~33 Mb); somatic loci miscounted as errors are ~1e-5 of all tabulated bases, far below the BQSR signal floor; the concern only bites for *ultra-hypermutators* (≥100 mut/Mb, POLE/MMR-deficient, ~10⁴ mutations), which this cohort is not; and (b) **unsafe for this dataset** — most patients (e.g. P005) have the CTRL sequenced years before the models, with the V8+UTR resequenced samples on a different run *and* capture chemistry (251125/V8 vs 230405/V6); applying the normal's run/chemistry-specific table to such a tumor corrects the *wrong* systematic-error profile, worse than per-sample BQSR. Per-sample BQSR is left as the deliberate catch-all; existing results are **not compromised on the BQSR axis**.
   - *Optional, non-blocking future-proofing:* if a genuinely ultra-hypermutated sample ever enters the cohort, the correct hardening is **not** normal-table-sharing but adding matched-normal germline (or eventual somatic) calls to the tumor's `--known-sites` mask — over-engineering for the current cohort, noted only for completeness.

4. ~~**No orientation bias modeling**~~ ✅ — fixed: added `--f1r2-tar-gz` output to `run_mutect2`, new `learn_read_orientation_model` rule, and `--orientation-bias-artifact-priors` to `filter_mutect2_calls` in `snv_calling_mutect2.smk`. Verified end-to-end on `P005_PT`: model trains (SUCCESS, all 32 contexts converged), `FilterMutectCalls` consumes the priors, and an `orientation` row now appears in the filtering stats with the filter applied in the VCF — at zero measurable sensitivity cost (FNR 0.0). Applied unconditionally (GATK default best-practice); old runs branched, not retro-applied.

5. ~~**No contamination estimation**~~ ✅ — fixed: added `contamination_resource` (`small_exac_common_3.hg38.vcf.gz`) to config, a generic per-sample `get_pileup_summaries` rule (capture-region ∩ common-SNP intersection, runs for tumour and matched normal), a paired/tumour-only `calculate_contamination` rule, and `--contamination-table` + `--tumor-segmentation` into `filter_mutect2_calls` (alongside #4's orientation priors). Verified end-to-end on `P005_PT` (paired, V6 CTRL / V6 PT): all rules `SUCCESS`, contamination 3.3e-4 ± 1.6e-4 (pure tumour, correctly filters nothing), a `contamination` row now present in the filtering stats. Applied unconditionally (GATK best-practice); cross-panel pairs (V6 normal / V8 tumour) tolerated by GATK with coarser segmentation, as noted for #3/#5.

6. ~~**HaplotypeCaller run on tumor BAM for het-SNPs**~~ ✅ — fixed (Design A) in `cnv_calling_cnvkit.smk`: het-SNP *discovery* is now normal-only via the new `extract_normal_het_sites` rule; the paired-tumor `gatk_haplotypecaller` force-calls those sites with `--alleles` (GATK 4.6 HaplotypeCaller — `--alleles` force-calls "regardless of evidence"; there is no `--force-call`, that is Mutect2-only) and restricts `-L` to them (also faster). Normal / tumour-only paths keep plain exome calling; `_is_paired_tumor` short-circuits tumour-only safely; merge/filter rule unchanged. Verified end-to-end on `P005_PT`: 48,543/49,347 normal het sites retained (98.4%; the ~1.6% lost are genuine low tumour coverage, `DP[1]≤10`, not silent CNV-region loss), and skew is preserved at imbalanced sites (e.g. `chr6:7288933` normal `0/1:55,50` → tumour `1/1:8,114`, BAF≈0.93 — exactly the LoH site class previously dropped). BAF track no longer blanked in aneuploid/LoH regions; gates PureCN (#12) and FACETS (#13) quality. *(Perf follow-up: #24.)*

7. ~~**`fastp` adapter removal is overlap-only (`--detect_adapter_for_pe` not set)**~~ ✅ — added config-gated `params.fastp.detect_adapter_for_pe` toggling the `--detect_adapter_for_pe` flag in the `fastp_trim` rule. Sequence-based detection now layers on top of (does not replace) the default R1/R2 overlap analysis — verified additive on `P005_CTRL`. Shipped default **`true`** for new runs; existing/old runs to be branched and left as-is (the flag alters trimmed output slightly, so it is not retroactively applied to already-processed results).

9. ~~**`MarkDuplicates` no-index intermediate, no `--intervals` in `BaseRecalibrator`**~~ ✅ — fixed in `bam_mapping_gatk.smk`: `mark_duplicates` now emits its `.bai` inline (`--CREATE_INDEX true`), and `create_base_recalibration` restricts `BaseRecalibrator` to the kit's capture BED (`covered_bedfile`) with `--interval-padding` (new `params.bqsr.interval_padding: 100`) instead of scanning the whole genome.

10. ~~**FastQC runs on post-BQSR BAM**~~ ✅ — fixed: `fastqc` moved to new `qc.smk` and now runs on the trimmed FASTQs (`{sample}_R1/R2.fq.gz`, output of `fastp_trim`) instead of the recalibrated BAM, so its chemistry plots reflect actual sequencing quality rather than GATK-rewritten scores.

11. ~~**Population AF filtering uses Mutect2's germline prior VCF**~~ ✅ — investigated, **sound and authority-aligned, not a defect**. Tumour-only filtering = PON (`run_mutect2`) + population-AF hard filter (`filter_population_af`, AF>0.001), the standard strategy used by the Chapuy group. The explicit filter is *additive, not redundant*: without a matched normal, `FilterMutectCalls`' germline filter is weak. The resource is gnomAD-derived `af-only-gnomad` rather than the cited 1000G — a deliberate, defensible substitution (gnomAD is far deeper; ≥ as thorough at the 0.001 cutoff). Known marginal limitation, consciously not fixed: gnomAD is blood-derived, so `AF>0.001` can drop bona-fide CH/recurrent-somatic drivers (`DNMT3A`, `JAK2 V617F`, …) in tumour-only mode. A cancerhotspots/COSMIC allele-exemption would address it but requires an off-DAG, hg19→hg38-liftover resource the user must hand-build — friction not worth the marginal gain; documented here rather than implemented.

14. ~~**`wildcard_constraints` blocks underscores in run IDs**~~ ✅ — the `{run}` constraint still forbids `_`/`.`/`/` (they collide with the pervasive `{run}_{sample}` filename scheme; full support would mean renaming that everywhere for no current benefit — no ID uses them). Converted the silent breakage into a loud early error: the Snakefile now validates every `ID` and raises with the offending value(s) before any DAG is built.

16. ~~**`CollectHsMetrics` missing**~~ ✅ — fixed: new `collect_hs_metrics` rule in `qc.smk` runs on the final recalibrated BAM against bait/target interval lists (new `bed_to_interval_list` rule, converting the kit's `covered_bedfile`/`target_regions_bedfile` split added to `probe_configs`). `multiqc` now depends on the parsed metric files (fastp json, fastqc zip, Picard metrics, mosdepth summary, hs_metrics) rather than rendered HTML, driven by a new `multiqc_config.yaml`.

17. ~~**CNVkit and Mutect2 PON reference construction not in Snakemake**~~ ✅ — done, realized as the standalone **`WES-PON-smk`** pipeline. Both reference-building workflows (which consume the same input — all normal BAMs, per capture kit) are unified there: `alignment`/`bqsr`/`coverage`/`qc` on the normals, then (1) `cnvkit_access` + `cnvkit_autobin` (with refFlat annotation, see #18), (2) `cnvkit_coverage` per normal per kit, (3) `cnvkit_reference` per kit/sex → replaces CNVkit-nf; (4) `genomicsdb_import` + `create_somatic_pon` → replaces the old Mutect2 PON pipeline. Adding a normal triggers a single re-run of `WES-PON-smk` to update all references. (Note: the old concern that `.cnn` gene labels carry the ENST bug is moot — they are now built with `--annotate refFlat`; and even before, `.cnn` labels are metadata only, `cnvkit.py fix` uses depth values.)

18. ~~**CNVkit targets BED has mixed gene annotation (ENST IDs at some loci)**~~ ✅ — fixed in `../WES-PON-smk`: targets/antitargets are now built fresh by `cnvkit_autobin` with `cnvkit.py ... --annotate {refflat}` (`cnvkit.smk:54`, refFlat at `config.yaml:14`) instead of reusing the legacy CNVkit-nf BEDs, so HGNC symbols are assigned from refFlat rather than inheriting Agilent `ens|ENST` tags. The original-analysis detail below is retained for context. (Capture-gap genes like **CDKN2A**/**STAT6** that have no probe coverage are still absent — annotation can't recover what isn't captured.)

    **Original issue:** the `targets.bed` files in the CNVkit-nf reference directories (`reference_V6/targets.bed`, `reference_V8/targets.bed`) were produced by the legacy CNVkit-nf Nextflow pipeline and carry raw Ensembl transcript IDs (`ENST...`) at a subset of loci instead of HGNC gene symbols. This happens because the Agilent `Regions.bed` files used as input contain `ens|ENST...` tags that propagate when `cnvkit.py target` is run without a gene-symbol annotation source (refFlat). Most intervals are fine, but loci where the Agilent BED has only `ens|` entries (no `ref|` gene symbol) end up with transcript IDs as their probe label. Those IDs are carried verbatim through the entire cnvkit chain (`.cnr` → `.cns` → `_CNVs.csv` → `combined_cnvs.rds`), causing affected genes to be silently absent from the callset by name even though their loci are covered. Known affected genes in the DLBCL driver list: **CDKN2A** (probe labeled `ENST00000404796` = MTAP, an adjacent gene), **STAT6** (probes labeled with neighboring gene transcripts — capture gap, not fixable by annotation alone).

    **Root cause:** `cnvkit.py target` was called without `--annotate refFlat.txt`. The Agilent BED mixed format was passed through as-is.

    **Fix:** Add a `cnvkit_annotate_targets` rule in `ref_index.smk` that runs before `cnvkit_coverage_and_fix`:

    ```python
    rule cnvkit_annotate_targets:
        input:
            targets=lambda w: config["probe_configs"][w.probe_version]["cnvkit_targets_raw"],
            refflat=config["refs"]["refflat"],
        output:
            bed="work/refs/cnvkit/{probe_version}/targets_annotated.bed",
        container:
            "docker://etal/cnvkit:0.9.11"
        log:
            "work/logs/cnvkit_annotate_targets_{probe_version}.log",
        shell:
            """
            cnvkit.py target {input.targets} \
                --annotate {input.refflat} \
                -o {output.bed} \
                > {log} 2>&1
            """
    ```

    Config changes required:
    - Rename `cnvkit_targets` → `cnvkit_targets_raw` in `probe_configs` (keeps the original file path)
    - Add `refflat: "/mnt/data/NGS/refs/UCSC/refFlat.txt"` under `refs:`
    - Update `cnvkit_coverage_and_fix` to use `work/refs/cnvkit/{probe_version}/targets_annotated.bed` as its targets input

    This makes annotation reproducible, version-controlled, and a tracked Snakemake dependency — updating refFlat automatically re-triggers the full cnvkit chain. The raw Agilent-derived files are never modified.

    **Partial workaround applied in `pdxpdcl_dlbcl` branch:** `cnvkit.py target --annotate refFlat.txt` was run manually on both kits, producing `targets_annotated.bed` alongside the originals. `config.yaml` on that branch points to the annotated files. A post-hoc ENST→symbol remapping was also applied in `combine_cnvs.R` (recovers PRDM1, SOCS1; cannot recover CDKN2A/STAT6 which were mislabeled at the BED level). See `docs/cnvkit_gene_annotation_issue.md` in the `pdxpdcl_dlbcl` branch for full details.

20. ~~**Add somalier — genotype-based sample-swap and sex QC**~~ ✅ — fixed: new `somalier_extract` (per-sample, off the final recalibrated BAM) / `somalier_ped` (per-run PED encoding samplesheet `Chr_sex`) / `somalier_relate` (per-run gather) chain in `workflow/rules/somalier.smk`, wired into `rule all` and into `multiqc`'s inputs so MultiQC's native somalier module picks up the report. Runs via conda (`somalier=0.3.2`), not a container — the `brentp/somalier:v0.3.2`/`v0.3.2-1` Docker images both ship a broken `/bin/bash` (readline symbol relocation failure), so bioconda was used instead, matching how this repo already handles other small CLI tools (bcftools, mosdepth, fastqc). Report-only, no DAG gating: mismatches surface in the HTML/MultiQC report for human review. Verified end-to-end on `P005` (CTRL/PT pair): both samples' inferred sex matches `Chr_sex` (XY → male), and CTRL/PT relatedness is 0.997 with 1.000 concordance, confirming correct tumor/normal pairing.

22. ~~**Manta correctness fixes**~~ ✅ — fixed in `sv_calling_manta.smk`: new `cap_manta_bam_quality` rule caps base qualities at 70 (pysam, new `envs/pysam.yaml`) on a Manta-only copy of each sample's recalibrated BAM (`work/manta/{run}/{sample}/{sample}.bqcap.bam`) — the BAM used by every other rule is untouched. `run_manta` now points `--tumorBam`/`--normalBam` at the capped BAMs, and its shell block pipes the raw `tumorSV`/`somaticSV` VCF through the bundled `convertInversion.py` (added `samtools=1.18` to `envs/manta.yaml`, which also supplies `bgzip`/`tabix` — `1.21`+ was tried and doesn't solve alongside `python=2.7`'s legacy zlib/ncurses pins, so `1.18` is the practical ceiling here) before re-compressing/indexing to the rule's declared output — folding paired-BND inversion records into proper `INV` calls ahead of `filter_manta_variants`/`annotate_sv`. Not yet re-verified end-to-end on a known-translocation case; do that before trusting new INV calls in a report.

23. ~~**Migrate SNV annotation from Funcotator to VEP**~~ ✅ — fixed: `funcotator` rule replaced by a `vep` rule (`snv_calling_mutect2.smk:241`, `docker://ensemblorg/ensembl-vep:release_116.0`), and `combine_results.R` rewired to parse the VEP CSQ field instead of Funcotator's FUNCOTATION field. Landed in the same commit that should have, but did not, update this file — caught retroactively while auditing commits against TODO.md. Leaves #19 (MAF output) stale, since it assumed Funcotator's native MAF emitter — see the note under #19.

26. ~~**Conda envs under-pinned — reproducibility risk**~~ ✅ — fixed: `envs/pysam.yaml` now pins `python=3.12` alongside `pysam=0.24.0` (previously no Python version was pinned, so the solver had picked a free-threaded `python-3.14.6` build that forced the GIL back on at import). `envs/r_vcf.yaml` now pins `r-base=4.3.3`, `bioconductor-variantannotation=1.48.1`, `r-tidyverse=2.0.0`, `r-openxlsx=4.2.8` (previously only `r-base=4.3` was pinned, the other three were unpinned). Verified via a full wet-run on host. `envs/purecn.yaml`'s unpinned `r-optparse` left as-is (moot while `purecn.smk` is dead code, #25). Conda build-string pinning remains an unresolved deeper gap, tracked separately in `portability-todo.md`.
