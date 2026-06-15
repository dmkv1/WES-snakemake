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

11. ~~**Population AF filtering uses Mutect2's germline prior VCF**~~ ✅ — investigated, **sound and authority-aligned, not a defect**. Tumour-only filtering = PON (`run_mutect2`) + population-AF hard filter (`filter_population_af`, AF>0.001), the standard strategy used by the Chapuy group. The explicit filter is *additive, not redundant*: without a matched normal, `FilterMutectCalls`' germline filter is weak. The resource is gnomAD-derived `af-only-gnomad` rather than the cited 1000G — a deliberate, defensible substitution (gnomAD is far deeper; ≥ as thorough at the 0.001 cutoff). Known marginal limitation, consciously not fixed: gnomAD is blood-derived, so `AF>0.001` can drop bona-fide CH/recurrent-somatic drivers (`DNMT3A`, `JAK2 V617F`, …) in tumour-only mode. A cancerhotspots/COSMIC allele-exemption would address it but requires an off-DAG, hg19→hg38-liftover resource the user must hand-build — friction not worth the marginal gain; documented here rather than implemented.

12. **PureCN receives pre-segmented `.cns` + `--fun-segmentation Hclust`** — contradictory; PureCN's own segmentation is better for purity estimation. Should feed `.cnr` ratios directly.

---

## Notable (correctness / maintainability)

13. **No orthogonal purity/ploidy estimator (FACETS)** — purity is currently a hardcoded samplesheet constant (CTRL=0, tumors=1); PureCN is the only estimator and its quality is gated by the BAF/het-SNP track (#6, #12). Add **FACETS** (WES tumor/normal, allele-specific CN + purity/ploidy from depth-ratio + BAF) as a second, orthogonal estimator so purity can be reported as a "two methods agree" number rather than an assumption. Cheap interim cross-check available with no new tool: `2 × modal VAF of clonal PASS SNVs in CN-neutral CNVkit segments ≈ purity` (Mutect2 VCF + CNVkit `.cns`, both already produced).

14. **`wildcard_constraints` blocks underscores in run IDs** — silent breakage if any samplesheet `ID` contains `_`.

15. **Read group `RGID`/`RGPU` parsing fragile** — encodes only instrument:flowcell, not lane; RGIDs can collide for multi-lane samples.

16. **`CollectHsMetrics` missing** — no capture efficiency, on-target rate, or fold-enrichment in MultiQC.

17. ~~**CNVkit and Mutect2 PON reference construction not in Snakemake**~~ ✅ — done, realized as the standalone **`../WES-PON-smk`** pipeline. Both reference-building workflows (which consume the same input — all normal BAMs, per capture kit) are unified there: `alignment`/`bqsr`/`coverage`/`qc` on the normals, then (1) `cnvkit_access` + `cnvkit_autobin` (with refFlat annotation, see #18), (2) `cnvkit_coverage` per normal per kit, (3) `cnvkit_reference` per kit/sex → replaces CNVkit-nf; (4) `genomicsdb_import` + `create_somatic_pon` → replaces the old Mutect2 PON pipeline. Adding a normal triggers a single re-run of `WES-PON-smk` to update all references. (Note: the old concern that `.cnn` gene labels carry the ENST bug is moot — they are now built with `--annotate refFlat`; and even before, `.cnn` labels are metadata only, `cnvkit.py fix` uses depth values.)

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

19. **No MAF-format SNV output for maftools / cohort-level analysis** — current SNV outputs are per-sample Funcotator-annotated VCF (`{sample}.SNV.vcf`) and the merged Excel workbook. Both are fine for per-sample inspection but unusable as input to `maftools` (oncoplots, mutational signatures, TMB summaries, somatic interactions, lollipop plots, etc.), which is the standard downstream entry point for cohort SNV analysis. Funcotator natively emits MAF (`--output-file-format MAF`), so this is a low-effort additive output — not a replacement for the existing VCF/xlsx.

    **Implementation roadmap:**

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

---

## Planned enhancements (from somatic-pipeline comparison — see `../comparison.md`)

Four features confirmed for implementation after benchmarking WES-smk against BALSAMIC-TGA and exome_somatic. **None of these touch the Mutect2 PON or the CNVkit reference build** — they sit downstream of, or orthogonal to, PON construction (see note at the end of this section). WGS is explicitly out of scope.

20. **Add somalier — genotype-based sample-swap and sex QC.** Current pairing and sex are taken from the samplesheet on trust: tumor↔normal pairing is never verified against the data, and `Chr_sex` is consumed blindly to select the CNVkit reference (`reference_male/female.cnn`) — a wrong `Chr_sex` silently picks the wrong reference and corrupts CNV calls with no error. BALSAMIC runs somalier relatedness (pair check) + sex prediction; WES-smk has no equivalent. Sample swaps are the most common wet-lab error and this is the cheapest control against them.
    - **Plan:** `somalier extract` per BAM (post-BQSR `{sample}.bam`) against a sites VCF → per-sample `.somalier`; `somalier relate` per run over a run's samples to confirm matched tumor/normal are the same individual; use somalier's inferred sex to validate the samplesheet `Chr_sex` (warn/gate on mismatch, since it drives CNVkit reference selection).
    - **Resources:** one sites VCF (`somalier` ships hg38 sites, e.g. `sites.hg38.vcf.gz`) — add under `refs:`. No PON dependency.
    - **Effort:** low. One extract rule (scatter per sample) + one relate rule (gather per run) + a small check feeding MultiQC (somalier has a native MultiQC module). No DAG restructuring.

21. **Add Delly as a second SV caller, merged with Manta.** WES-smk is Manta-only; BALSAMIC runs Manta + Delly merged via SVDB. **Cohort-specific payoff:** DLBCL is defined by IG translocations (BCL2/BCL6/MYC) and Manta and Delly have non-overlapping breakpoint sensitivity — a consensus/union materially improves translocation recall. This is the highest-value caller addition for this cohort (higher than a second SNV caller).
    - **Plan:** new `run_delly` rule (`delly call -x <exclude.bed>`, paired = tumor+normal with somatic pre-filter `delly filter -f somatic -s samples.tsv`; tumor-only = single-sample); PASS-filter; merge Manta + Delly with SVDB (`--no_intra --bnd_distance 5000 --overlap 0.80`, priority list) or a lighter `bcftools`/survivor merge; feed the merged VCF to the existing `annotate_sv` (AnnotSV) rule in `sv_calling_manta.smk`.
    - **Resources:** Delly exclude-regions BED for hg38 (ships with Delly). No PON required for somatic Delly. (Delly-CNV is a separate optional path — not in scope here; CNV stays CNVkit/PureCN.)
    - **Effort:** medium. New caller + merge step + AnnotSV input swap. Verify on a known-translocation case.

22. **Manta correctness fixes.** Two small items BALSAMIC applies that WES-smk's `run_manta` omits:
    - **`convertInversion.py`** — Manta represents inversions as paired BND records; running the bundled `convertInversion.py` post-call folds them into proper `INV` records, giving correct SV-type and cleaner AnnotSV interpretation. Slots between `run_manta` and `filter_manta_variants`.
    - **Base-quality cap (→70) on the Manta input BAM** — BALSAMIC caps base qualities for Manta. Directly relevant here given the binned NextSeq qualities documented in #1; uncapped binned/odd qualities can skew Manta's scoring. A small pysam pass on the BAM feeding `run_manta` (does not affect the BAM used elsewhere).
    - **Effort:** low for both. Pure pre/post-processing around the existing Manta rule; no new resources.

23. **Migrate SNV annotation from Funcotator to VEP — important for DLBCL.** WES-smk uses GATK Funcotator (`funcotator` rule, `snv_calling_mutect2.smk:241`); both comparison pipelines use Ensembl VEP. For DLBCL specifically VEP is the better engine: richer/standard consequence annotation (`--everything`, HGVS, canonical/MANE), and a clean path to the layered annotation DLBCL interpretation needs — COSMIC, cancerhotspots, and the CH-driver exemption noted in #11 — via `--custom`/vcfanno overlays that Funcotator can't easily express. VEP is also the lingua franca for downstream DLBCL tooling.
    - **Plan:** replace the `funcotator` rule with a `vep` rule (`vep --offline --cache --everything --hgvsg --variant_class --vcf --assembly GRCh38 --fork N`, FASTA-based) → `{sample}.SNV.vep.vcf.gz`. Optionally layer `--custom COSMIC` + a cancerhotspots overlay to action #11's CH-driver rescue at the same time.
    - **Interaction with #19 (MAF):** #19 assumed Funcotator's native MAF emitter. With VEP, the cohort MAF comes from `vcf2maf` (VEP-based) instead — re-scope #19 to `vcf2maf` once VEP lands, or implement #19 first on Funcotator and migrate the MAF step with this item. Decide ordering before starting either.
    - **Resources:** VEP cache (merged, GRCh38) + reference FASTA (already present); optional COSMIC VCF + cancerhotspots list. **Replaces** the Funcotator data-sources dependency; no PON impact.
    - **Effort:** medium. New rule + cache install + re-wire `combine_results.R` / `integrate_results.smk` to the VEP VCF's CSQ field instead of Funcotator's FUNCOTATION field (the downstream parser change is the real work, not the VEP call itself). Verify the result tables reproduce the expected fields on `P005_PT`.

### PON impact of items 20–23

**None of these four change the PON pipeline.** The Mutect2 PON (`GenomicsDBImport` + `CreateSomaticPanelOfNormals`) and the CNVkit reference (`cnvkit.py coverage`/`reference` on the normals) are built **upstream** of everything 20–23 touches:
- somalier (#20) reads the finished normal BAMs for QC; it does not feed PON construction.
- Delly (#21) and the Manta fixes (#22) are SV-side and use the BAMs directly, no PON.
- VEP (#23) is pure annotation, downstream of calling.

The normal BAMs that feed the PON come out of the alignment→dedup→BQSR chain (`bam_mapping_gatk.smk`), and **none of 20–23 alter that chain** — so a PON built now stays byte-for-byte valid after these features land. **Safe to build the PON on all current normals now.**

No outstanding reference-build caveat either: the unified reference pipeline (#17) already exists as **`../WES-PON-smk`** and incorporates the refFlat CNVkit fix (#18) — `cnvkit_autobin` runs `cnvkit.py ... --annotate {refflat}` (`cnvkit.smk:54`), so a single run rebuilds **both** a correct Mutect2 PON (`genomicsdb_import` + `create_somatic_pon`) and a correctly-annotated CNVkit reference (`cnvkit_reference`, per probes/sex). **Safe to run `WES-PON-smk` on all current normals now** — nothing in 20–23 invalidates its outputs.
