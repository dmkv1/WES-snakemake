# WES-snakemake Pipeline — Known Issues & TODO

## Critical (affect result correctness)

1. ~~**No UMI-aware trimming or deduplication for V8+UTR samples**~~ ✅ — investigated, **not a pipeline bug; no sample requires UMI removal**. The V8+UTR prep (SureSelect XT HS, G9702-90005) does generate a 10-bp i5 MBC, but CLIP/BaseSpace demux delivers FASTQs with the MBC **dropped at FASTQ generation**: i5 written as a constant 8-bp *sample* index, no I2 read, raw binned NextSeq qualities, PCR-duplicate load consistent with un-collapsed data (triangulated across three independent checks). The barcode is physically absent from every delivered FASTQ, so coordinate-based `MarkDuplicates` is the only valid option and is exactly what runs — existing results are **not compromised on the UMI axis**. `fastp` defaults (PE adapter overlap-trim, auto polyG for NextSeq, quality/N/length filter) are adequate for non-UMI WES. Recovery is possible *only* by regenerating FASTQs from BaseSpace with i5/MBC retrieval (facility request, not a code change). Detail: memory `reference_wes_umi_demux.md`.
   - *Optional, non-blocking future-proofing:* add an explicit `umi` samplesheet column (orthogonal to `probes`) gating a dormant fgbio path; set `umi=false` for all current rows (byte-identical behavior, no remap) so correctly-demuxed future runs work without touching existing results.

2. ~~**BWA-MEM `-M` flag**~~ ✅ — fixed: replaced with `-Y` in `bam_mapping_gatk.smk`.

3. ~~**BQSR trained on tumor BAM**~~ ✅ — investigated, **known caveat, not a pipeline bug; per-sample BQSR is the correct catch-all here**. Current pipeline runs `BaseRecalibrator` independently per sample (`bam_mapping_gatk.smk:144–205`), exactly as GATK4 Best Practices prescribes. The GATK3-era "apply matched-normal's recal table to the tumor" advice is (a) **negligible at this mutation load** — coding TMB ≈ 4.5–12 mut/Mb (150–400 CDS SNVs / ~33 Mb); somatic loci miscounted as errors are ~1e-5 of all tabulated bases, far below the BQSR signal floor; the concern only bites for *ultra-hypermutators* (≥100 mut/Mb, POLE/MMR-deficient, ~10⁴ mutations), which this cohort is not; and (b) **unsafe for this dataset** — most patients (e.g. P005) have the CTRL sequenced years before the models, with the V8+UTR resequenced samples on a different run *and* capture chemistry (251125/V8 vs 230405/V6); applying the normal's run/chemistry-specific table to such a tumor corrects the *wrong* systematic-error profile, worse than per-sample BQSR. Per-sample BQSR is left as the deliberate catch-all; existing results are **not compromised on the BQSR axis**.
   - *Optional, non-blocking future-proofing:* if a genuinely ultra-hypermutated sample ever enters the cohort, the correct hardening is **not** normal-table-sharing but adding matched-normal germline (or eventual somatic) calls to the tumor's `--known-sites` mask — over-engineering for the current cohort, noted only for completeness.

4. ~~**No orientation bias modeling**~~ ✅ — fixed: added `--f1r2-tar-gz` output to `run_mutect2`, new `learn_read_orientation_model` rule, and `--orientation-bias-artifact-priors` to `filter_mutect2_calls` in `snv_calling_mutect2.smk`. Verified end-to-end on `P005_PT`: model trains (SUCCESS, all 32 contexts converged), `FilterMutectCalls` consumes the priors, and an `orientation` row now appears in the filtering stats with the filter applied in the VCF — at zero measurable sensitivity cost (FNR 0.0). Applied unconditionally (GATK default best-practice); old runs branched, not retro-applied.

5. ~~**No contamination estimation**~~ ✅ — fixed: added `contamination_resource` (`small_exac_common_3.hg38.vcf.gz`) to config, a generic per-sample `get_pileup_summaries` rule (capture-region ∩ common-SNP intersection, runs for tumour and matched normal), a paired/tumour-only `calculate_contamination` rule, and `--contamination-table` + `--tumor-segmentation` into `filter_mutect2_calls` (alongside #4's orientation priors). Verified end-to-end on `P005_PT` (paired, V6 CTRL / V6 PT): all rules `SUCCESS`, contamination 3.3e-4 ± 1.6e-4 (pure tumour, correctly filters nothing), a `contamination` row now present in the filtering stats. Applied unconditionally (GATK best-practice); cross-panel pairs (V6 normal / V8 tumour) tolerated by GATK with coarser segmentation, as noted for #3/#5.

6. ~~**HaplotypeCaller run on tumor BAM for het-SNPs**~~ ✅ — fixed (Design A) in `cnv_calling_cnvkit.smk`: het-SNP *discovery* is now normal-only via the new `extract_normal_het_sites` rule; the paired-tumor `gatk_haplotypecaller` force-calls those sites with `--alleles` (GATK 4.6 HaplotypeCaller — `--alleles` force-calls "regardless of evidence"; there is no `--force-call`, that is Mutect2-only) and restricts `-L` to them (also faster). Normal / tumour-only paths keep plain exome calling; `_is_paired_tumor` short-circuits tumour-only safely; merge/filter rule unchanged. Verified end-to-end on `P005_PT`: 48,543/49,347 normal het sites retained (98.4%; the ~1.6% lost are genuine low tumour coverage, `DP[1]≤10`, not silent CNV-region loss), and skew is preserved at imbalanced sites (e.g. `chr6:7288933` normal `0/1:55,50` → tumour `1/1:8,114`, BAF≈0.93 — exactly the LoH site class previously dropped). BAF track no longer blanked in aneuploid/LoH regions; gates PureCN (#12) and FACETS (#13) quality.

---

## Significant (suboptimal but not wrong)

7. ~~**`fastp` adapter removal is overlap-only (`--detect_adapter_for_pe` not set)**~~ ✅ — added config-gated `params.fastp.detect_adapter_for_pe` toggling the `--detect_adapter_for_pe` flag in the `fastp_trim` rule. Sequence-based detection now layers on top of (does not replace) the default R1/R2 overlap analysis — verified additive on `P005_CTRL`. Shipped default **`true`** for new runs; existing/old runs to be branched and left as-is (the flag alters trimmed output slightly, so it is not retroactively applied to already-processed results).

8. **xengsort mouse-read stats not captured** — fraction of murine reads is critical QC for PDX samples and is currently lost.

9. **`MarkDuplicates` no-index intermediate, no `--intervals` in `BaseRecalibrator`** — minor sort-order inefficiency in the BWA→GATK chain.

10. **FastQC runs on post-BQSR BAM** — recalibrated quality scores make FastQC plots meaningless for detecting sequencing chemistry problems. Should run on raw or trimmed FASTQs.

11. **Population AF filtering uses Mutect2's germline prior VCF** (`af-only-gnomad`) — modified resource, not standard gnomAD; the 0.001 AF threshold is largely redundant with what `FilterMutectCalls` already does internally.

12. **PureCN receives pre-segmented `.cns` + `--fun-segmentation Hclust`** — contradictory; PureCN's own segmentation is better for purity estimation. Should feed `.cnr` ratios directly.

---

## Notable (correctness / maintainability)

13. **No orthogonal purity/ploidy estimator (FACETS)** — purity is currently a hardcoded samplesheet constant (CTRL=0, tumors=1); PureCN is the only estimator and its quality is gated by the BAF/het-SNP track (#6, #12). Add **FACETS** (WES tumor/normal, allele-specific CN + purity/ploidy from depth-ratio + BAF) as a second, orthogonal estimator so purity can be reported as a "two methods agree" number rather than an assumption. Cheap interim cross-check available with no new tool: `2 × modal VAF of clonal PASS SNVs in CN-neutral CNVkit segments ≈ purity` (Mutect2 VCF + CNVkit `.cns`, both already produced).

14. **`wildcard_constraints` blocks underscores in run IDs** — silent breakage if any samplesheet `ID` contains `_`.

15. **Read group `RGID`/`RGPU` parsing fragile** — encodes only instrument:flowcell, not lane; RGIDs can collide for multi-lane samples.

16. **`CollectHsMetrics` missing** — no capture efficiency, on-target rate, or fold-enrichment in MultiQC.
