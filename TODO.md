# WES-snakemake Pipeline — Known Issues & TODO

## Critical (affect result correctness)

1. ~~**No UMI-aware trimming or deduplication for V8+UTR samples**~~ ✅ — investigated, **not a pipeline bug; no sample requires UMI removal**. The V8+UTR prep (SureSelect XT HS, G9702-90005) does generate a 10-bp i5 MBC, but CLIP/BaseSpace demux delivers FASTQs with the MBC **dropped at FASTQ generation**: i5 written as a constant 8-bp *sample* index, no I2 read, raw binned NextSeq qualities, PCR-duplicate load consistent with un-collapsed data (triangulated across three independent checks). The barcode is physically absent from every delivered FASTQ, so coordinate-based `MarkDuplicates` is the only valid option and is exactly what runs — existing results are **not compromised on the UMI axis**. `fastp` defaults (PE adapter overlap-trim, auto polyG for NextSeq, quality/N/length filter) are adequate for non-UMI WES. Recovery is possible *only* by regenerating FASTQs from BaseSpace with i5/MBC retrieval (facility request, not a code change). Detail: memory `reference_wes_umi_demux.md`.
   - *Optional, non-blocking future-proofing:* add an explicit `umi` samplesheet column (orthogonal to `probes`) gating a dormant fgbio path; set `umi=false` for all current rows (byte-identical behavior, no remap) so correctly-demuxed future runs work without touching existing results.

2. ~~**BWA-MEM `-M` flag**~~ ✅ — fixed: replaced with `-Y` in `bam_mapping_gatk.smk`.

3. **BQSR trained on tumor BAM** — known GATK caveat; for hypermutated lymphoma samples this is non-trivial. Matched normal's recal table should be applied to the tumor.

4. **No orientation bias modeling** — `LearnReadOrientationModel` not run; no `--f1r2-tar-gz` output from Mutect2; no `--orientation-bias-artifact-priors` in `FilterMutectCalls`. FFPE and oxidative damage artifacts (C>T) pass as PASS somatic variants.

5. **No contamination estimation** — `GetPileupSummaries` + `CalculateContamination` missing before `FilterMutectCalls`. Without this, Mutect2 filtering runs on incorrect priors; contaminated samples pass more false positives.

6. **HaplotypeCaller run on tumor BAM for het-SNPs** — diploid caller on an aneuploid tumor miscalls/misses het SNPs in CNV regions, degrading the BAF track precisely where it matters most. Should use normal BAM only.

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

13. **`wildcard_constraints` blocks underscores in run IDs** — silent breakage if any samplesheet `ID` contains `_`.

14. **Read group `RGID`/`RGPU` parsing fragile** — encodes only instrument:flowcell, not lane; RGIDs can collide for multi-lane samples.

15. **`CollectHsMetrics` missing** — no capture efficiency, on-target rate, or fold-enrichment in MultiQC.
