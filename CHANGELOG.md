# Changelog

This project uses [Semantic Versioning](https://semver.org), adapted for a pipeline.
The public contract is the samplesheet schema, the `config.yaml` keys and the output
paths. See [Versioning](README.md#versioning) in the README.

Each release states whether it changes results for the same input data. A version that
changes results needs a re-run before you compare old and new cohorts.

## [2.0.0] - 2026-08-20

**This release breaks existing setups.** The `resources` block of `config.yaml` is
restructured: three keys are removed and the GATK rules read a tier instead, so a 1.3.0
`config.yaml` raises a `KeyError` while Snakemake parses the rules, before any job runs.
Read the Breaking changes section before you upgrade. No samplesheet column, no
reference file and no output path changes.

**This release changes `results/combined/cohort.maf`.** Calls do not change: the same
input data gives the same VCFs, and every other `combined_*` table is byte-identical.
The MAF gains three columns and re-classifies a handful of consequence terms, so re-run
`cohort_maf` before you compare a MAF against a 1.3.0 cohort.

### Breaking changes

`config.yaml` keys:

| Old | New |
|---|---|
| `resources.java_max_gb`, `resources.java_min_gb`, `resources.mem_mb` | `resources.gatk.<tier>.java_max_gb`, `.java_min_gb` and `.mem_mb`, for the tiers `light`, `medium`, `markdup` and `heavy` |

New required `config.yaml` keys, which have no defaults: `resources.gatk.*`,
`resources.mutect2_threads`, `resources.haplotypecaller_threads`,
`resources.manta_threads`, `resources.vep_threads` and `resources.delly_mem_mb`.

Copy the `resources` block from `config.yaml.example` into an existing `config.yaml`
rather than patching key by key.

### Added

* `resources.gatk.light`, `.medium`, `.markdup` and `.heavy` in `config.yaml`. Each is a
  `java_min_gb` / `java_max_gb` / `mem_mb` triple. Every GATK rule names a tier instead
  of the one global heap, so `ApplyBQSR` no longer reserves what `Mutect2` reserves.
  Ceilings are about 2.5x the `max_rss` the rule's benchmark recorded on a 9 GB WES BAM.
* `resources.mutect2_threads`, `resources.haplotypecaller_threads`,
  `resources.manta_threads`, `resources.vep_threads` and `resources.delly_mem_mb`. They
  cap the tools that stop scaling before `resources.threads`; Snakemake lowers each to
  the cores the host has, so they need no editing for a smaller machine.
* A memory reservation on the five non-GATK rules that had none: `bwa_map` 12 GB,
  `purecn_run` 16 GB, `cnvkit_coverage_and_fix` 8 GB, `vep` 4 GB and
  `cap_manta_bam_quality` 1 GB. These are set in the rules, not in `config.yaml`, because
  they do not scale with the host. The scheduler had been booking each at the profile's
  4 GB `default-resources` figure.
* An `io_heavy` resource, capped at 4 in `profiles/default/config.yaml`. `bwa_map`,
  `mark_duplicates`, `sort_bam`, `apply_base_recalibration` and `cap_manta_bam_quality`
  each declare 1. Their memory tiers alone would let many run at once, and on a
  spinning-disk array that many concurrent whole-BAM rewrites seek-thrash rather than
  deliver throughput. Raise it on NVMe.
* `benchmark:` on 28 rules, writing `work/benchmarks/<rule>/<run>_<sample>.tsv`, one file
  per unit for the rules that run per unit and a single file for `multiqc` and
  `build_xengsort_index`. `benchmark-extended: true` in the profile puts the declared
  resources and the `max_rss` they were measured against in the same row.
* `cohort.maf` gains `Gene`, the Ensembl gene ID, and `vcf_position` / `vcf_variant`.
  MAF coordinates are re-derived to MAF convention and are not a valid join key back to
  `combined_snvs.tsv`; join on `Tumor_Sample_Barcode` plus the two `vcf_*` columns.

### Changed

* `resources.manta_max_gb` 32 to 6 and `resources.manta_mem_mb` 40960 to 8192. Manta
  parcels the call out to its own scheduler under `--jobs` and holds under 1 GB in the
  parent, so `--memGb` is a ceiling for that scheduler, not a heap it fills.
* `resources.xengsort_index_mem_mb` 65536 to 32768. The index build measured 21.7 GB.
* `run_delly` reserves `resources.delly_mem_mb`, 4 GB, instead of the global
  `resources.mem_mb`, 40 GB. It measured 1.1 GB.
* `sort_bam` reserves 30% over its per-thread buffers instead of 15%. The buffers
  themselves, the merge pass and the BGZF output measured about 11% over on a WES BAM;
  the rest is margin.
* `profiles/default/config.yaml` moves to 32 cores and 245760 MB, half of a 64-thread,
  504 GB host. Scale it with the `resources` block, as before.

### Fixed

* `combined_to_maf.R` maps ten SO terms it previously did not: `start_retained_variant`,
  `splice_donor_5th_base_variant`, `splice_donor_region_variant`,
  `splice_polypyrimidine_tract_variant`, `TFBS_ablation`, `TFBS_amplification`,
  `regulatory_region_ablation`, `regulatory_region_amplification`,
  `feature_elongation` and `feature_truncation`. A variant whose only consequence was
  one of these was written as `Targeted_Region`; it now gets `Silent`, `Splice_Region`
  or `IGR`. In a multi-term consequence list the term also ranked last instead of at its
  vcf2maf priority. The SnpEff-only synonyms vcf2maf carries are still left out, and the
  header comment now says so, so the next re-port does not read the gap as an oversight.
* An SO term with no entry in either table is now collected and reported once per run.
  It was warned about per variant, and R collapses anything past 50 warnings into "There
  were 50 or more warnings", which loses the term names.
* `HGVSp_Short` substitutes only whole three-letter codons. It replaced any occurrence of
  a codon name anywhere in the string.
* The read-group warnings print once per Snakemake invocation. Rules with a `run:`
  directive re-import the Snakefile in a worker process, so the list was repeated ten
  times in the main log.

### Documentation

* The README describes the profile the repository ships, not the 16-thread, 64 GB host
  it was first sized for.
* `config.yaml.example` states where each memory figure comes from and which benchmark
  to re-read after a cohort at different coverage.

## [1.3.0] - 2026-08-19

**This release does not change results.** The same input data gives the same calls as
1.2.0. It changes what the scheduler believes each rule costs, so a run on a large host
no longer over-subscribes CPU or memory.

### Added

* `resources.xengsort_mem_mb` and `resources.xengsort_index_mem_mb` in `config.yaml`.
  Both are required. `run_xengsort` and `build_xengsort_index` declare them, so the
  scheduler stops falling back to the profile's 4 GB default. xengsort holds the whole
  human and mouse hash in the process and shares nothing between concurrent jobs, so the
  footprint multiplies by the number of slots.

### Changed

* `run_manta` and `run_delly` declare their CPUs with `threads:` instead of a
  `resources: threads=` entry. `resources.threads` is not Snakemake's CPU reservation,
  so both rules were booked as one core while running many. `run_delly` now asks for one
  thread per input BAM, two for a pair and one tumour-only, which is all its OpenMP loop
  over the BAMs can use.
* The quality-capped BAMs that Manta reads, `work/manta/{run}/{sample}/{sample}.bqcap.bam`,
  are `temp()`. A run's normal is held until the last tumour of that run is called.
  Re-deriving one costs a single-threaded pysam pass over the full BAM, around 30 to 40
  minutes for 10 GB, so pass `--notemp` when a Manta rerun over already-called samples is
  expected.

### Fixed

* `run_manta` declares the `.bai` of its tumour and normal inputs. `configManta.py`
  resolves the indices itself, so they were undeclared, and `temp()` reaped them before
  the calling job ran.

## [1.2.0] - 2026-08-10

**This release changes results for tumour-normal pairs whose two samples use different
capture kits.** Every other pair and every tumour-only sample is unaffected. If the
cohort mixes kits within a pair, re-run Mutect2 before you compare against a 1.1.0
cohort.

### Changed

* Mutect2 in paired mode now calls on the intersection of the tumour's and the matched
  normal's capture kits, instead of the tumour's kit alone. Outside the intersection the
  normal has no coverage and cannot veto a call, and the Panel of Normals is not applied
  in paired mode, so germline sites passed unopposed. A pair on one kit passes the same
  BED twice, which is a no-op.

### Fixed

* The rules that write to `tmp` now create it. `mark_duplicates`, `sort_bam`,
  `create_base_recalibration`, `apply_base_recalibration`, `bed_to_interval_list` and
  `collect_hs_metrics` failed on a fresh checkout, where the directory does not exist
  yet. `mark_duplicates` also passed a hardcoded `tmp` to `--TMP_DIR` instead of its own
  `tmp_dir` parameter.
* `mosdepth` runs with `--no-per-base`. It was writing a per-base depth track that no
  rule declares and no downstream step reads, which cost runtime and disk per sample.
  The three declared outputs, the summary, the region distribution and the thresholds
  BED, are unchanged.

### Removed

* `workflow/scripts/fastq_header.py` is no longer a vendored copy of `wesingest/header.py`
  from the lab's WES repository. The parser is self-contained, so the workflow stays
  clone-and-run without that repository present. `tests/test_vendored_parser.py`, which
  compared the two sources object by object, is deleted along with the coupling. Parsing
  behaviour does not change: the commit touches docstrings, comments and the removed
  test only.

### Documentation

* The DELLY exclusion template carries both `chr`-prefixed and unprefixed contig names.
  The README claimed it was `chr`-prefixed only.
* The README links `WES-PON-smk` to its repository.
* `units.py` and `fastq_header.py` describe the samplesheet shapes they accept without
  reference to `wesingest`. A bare `sample,fq1,fq2` sheet and a fully specified one are
  both first-class, as before.

## [1.1.0] - 2026-08-07

**This release does not change results.**

Host configuration moves out of version control.

### Changed

* The repository tracks `config.yaml.example` and ignores `config.yaml`. A new checkout
  needs one extra step before the first run:

  ```bash
  cp config.yaml.example config.yaml
  ```

  Then repoint the `/path/to/...` placeholders. An existing `config.yaml` is untouched.
  It becomes an ignored file, so host paths no longer appear in a commit.

### Note

A checkout of `v1.0.0` or earlier overwrites the working `config.yaml`, because those
commits track the file. Keep a copy of the host configuration outside the repository
before you check out an older version.

## [1.0.0] - 2026-08-07

The first tagged release. It adds a second SV caller, a purity model, sample identity QC,
multi-lane support and a test suite. It replaces Funcotator with VEP.

**This release changes results.** The alignment, the duplicate marking and the annotation
all changed. Do not compare its output against a 0.1.0 cohort. Re-run instead.

**This release breaks existing setups.** The samplesheet columns, the `config.yaml` keys
and several output paths changed. Read the Breaking changes section before you upgrade.

### Breaking changes

Samplesheet columns:

| Old | New | Note |
|---|---|---|
| `Chr_sex`, values `XX` and `XY` | `gender`, values `f` and `m` | The name and the vocabulary both changed |
| `probes` | `capture_kit` | The values are still `probe_configs` keys |
| `purity` | `tumor_fraction` | The value must be in `(0,1]`. Blank or `NA` means unknown. `0` is no longer valid |

A run `ID` must not contain `_`, `.` or `/`. The pipeline stops at startup if it does,
because output files use the `{ID}_{sample}` pattern.

`config.yaml` keys:

| Old | New |
|---|---|
| `tumor_only.pon.<kit>` | `panel_of_normals.mutect2.<kit>` |
| `probe_configs.<kit>.regions_bedfile` | `probe_configs.<kit>.covered_bedfile`, which is the Covered BED, not the Regions BED |
| `probe_configs.<kit>.cnvkit_ref_m` and `cnvkit_ref_f` | `panel_of_normals.cnvkit.cnvkit_ref_m.<kit>` and `cnvkit_ref_f.<kit>`, with the sex key first |
| `probe_configs.<kit>.library_prep` | Removed. The `LB` tag comes from the `library` column or the sample name |
| `refs.funcotator_data_sources.*` | `refs.vep_cache.*` |

New required `config.yaml` keys, which have no defaults: `containers.*`,
`panel_of_normals.path`, `probe_configs.<kit>.target_regions_bedfile`,
`refs.contamination_resource`, `refs.vep_cache.*`, `refs.purecn.snp_blacklist`,
`refs.somalier_sites`, `params.cnv.use_purecn_purity`, `params.fastp.*`, `params.bqsr.*`,
`params.delly.*`, `params.somalier.*`, `resources.mem_mb`, `resources.sort_mem` and
`resources.manta_mem_mb`.

New required reference files: `small_exac_common_3.hg38.vcf.gz`, the VEP offline cache
for GRCh38 v116, `sites.hg38.vcf.gz`, `human.hg38.excl.tsv`, `hg38_simpleRepeats.bed`,
and the Covered BED per capture kit.

New required Panel of Normals artifacts from `WES-PON-smk`: `normalDB_<kit>_{m,f}_hg38.rds`
and `mapping_bias_<kit>_{m,f}_hg38.rds`. The CNVkit references are renamed from
`reference_male.cnn` and `reference_female.cnn` to `reference_m.cnn` and `reference_f.cnn`.

Output paths:

| Old | New |
|---|---|
| `results/combined/combined_data.rds` | `results/combined/combined_{snvs,cnvs,svs,qc}.tsv` |
| `results/{run}/{sample}/{sample}_{SNVs,CNVs,SVs}.csv` | `work/combine/{run}/{sample}/`, no longer deliverables |
| `results/metrics/{run}_{sample}.mosdepth.*` | `results/metrics/{run}/{sample}.mosdepth.*` |
| `results/qc/fastp/{run}/{sample}_fastp.*` | `results/qc/fastp/{run}/{sample}.{unit}_fastp.*` |
| `results/qc/fastqc/{run}/{sample}_fastqc.html` | `results/qc/fastqc/{run}/{sample}.{unit}_R{1,2}_fastqc.*` |
| `work/manta/{run}/{sample}/{sample}.SV.annotated.tsv` | `work/sv_merge/{run}/{sample}/{sample}.SV.annotated.tsv` |

`results/{run}/{sample}/{sample}.SNV.vcf.idx` is no longer produced, because VEP writes
no index.

The conda environment is renamed from `snakemake-wes` to `snakemake`.

### Results change

These changes alter the output for the same input data:

* bwa mem uses `-Y -K 100000000` in place of `-M`. Supplementary alignments are now
  soft-clipped, and the result no longer depends on the thread count. Every BAM changes.
* MarkDuplicates runs on query-grouped input, so it marks supplementary and secondary
  reads differently.
* BQSR is restricted to the Covered BED with `--interval-padding 100`.
* Mutect2 filtering now uses orientation-bias priors and a contamination estimate, so the
  PASS set changes.
* VEP replaces Funcotator, so every annotation column changes. The `Gencode_*`,
  `ClinVar_VCF_*`, `Cosmic_*`, `HGNC_*` and `dbSNP_*` columns are replaced by VEP CSQ
  fields such as `SYMBOL`, `Consequence`, `IMPACT`, `HGVSc`, `HGVSp`, `CLIN_SIG`, `SIFT`,
  `PolyPhen` and `gnomADe_*`.
* `cnvkit call` takes the purity and the ploidy from the purity resolver, in place of the
  raw samplesheet `purity` with no ploidy.
* Manta reads a base-quality-capped BAM, and its paired breakend records become
  inversions.
* SVs are now a Manta and DELLY consensus, not Manta alone.

### Added

* **DELLY v2.5.1** as a second SV caller, with the author's short-read settings.
* **SURVIVOR 1.0.7** merges the Manta and DELLY calls into a consensus set. `SUPP_VEC`
  records the caller support. AnnotSV now annotates the merged set.
* **PureCN** estimates the tumor purity and the ploidy for paired runs. A resolver rule
  chooses between the samplesheet value, the PureCN estimate and an assumed pure sample,
  and writes one sidecar file that both `cnvkit call` and the report read.
* **somalier 0.3.2** sample identity QC. It relates the whole cohort, every sample
  against every other sample, and writes a relatedness table and a heatmap. It reports
  only and does not stop the run.
* **Multi-lane support.** A samplesheet row expands into one alignment unit per FASTQ
  pair. Each unit carries a derived and validated read group, and MarkDuplicates gathers
  the units. New optional columns `flowcell`, `lane`, `barcode` and `library` override
  the derived values. `results/metadata/units.tsv` records the provenance.
* **Mutect2 orientation-bias and contamination modelling** through
  `LearnReadOrientationModel`, `GetPileupSummaries` and `CalculateContamination`.
* **GATK CollectHsMetrics** capture metrics.
* **Cohort outputs**: `results/combined/cohort.maf` and four cohort TSV tables.
* **A test suite** with 71 tests over the samplesheet parser and the read-group logic,
  plus a GitHub Actions workflow and a coverage badge.
* Base quality capping for Manta, and inversion conversion for its breakend records.
* `results/metadata/units.tsv` and `results/metadata/samples.tsv`.
* Startup validation for run and sample names, for the gender value, for the tumor
  fraction range, and for the presence of a Mutect2 PON on a tumor-only run.

### Changed

* FastQC now reads the trimmed FASTQ files per unit, not the final BAM, because BQSR
  makes the FastQC chemistry plots meaningless.
* MultiQC gained a row-naming configuration so a multi-lane sample reads correctly.
* The container images moved from hardcoded strings into `config["containers"]`.
* The Singularity bind mounts moved from the profile into the Snakefile, which builds
  them from `refs.path` and `panel_of_normals.path`.
* `launch.sh` records the Snakemake process ID, and `stop.sh` reports a clear error when
  that file is absent.
* The conda environments pin their package versions.
* GATK 4.6.1.0 to 4.6.2.0. CNVkit 0.9.11 to 0.9.14. MultiQC 1.31 to 1.33.
* The shipped `config.yaml` holds placeholder paths, not the paths of one host.
* The resource defaults assume a 16-thread host with 64 GB of RAM, in place of the
  64-thread host they assumed before. Raise them for a larger machine.
* The README was rewritten against the code.

### Removed

* **GATK Funcotator**, replaced by VEP.
* **GATK AddOrReplaceReadGroups** and **FixMateInformation**. bwa writes the read group
  inline, and `samtools fixmate -m` does the mate fixing.
* `cnvkit export seg`, and the unused `cnvkit_call_with_purity` and
  `cnvkit_plots_calibrated` rules.
* The `bundle_data.R`, `combine_snvs.R`, `combine_svs.R` and `combine_cnvs.R` scripts.

### Fixed

* The CNVkit male reference detection read the reference filename. It now uses the sample
  sex directly.
* The xengsort and mosdepth rows did not merge in the MultiQC report.
* The tumor B-allele frequency came from a HaplotypeCaller force call, which needs
  assembly. It now comes from `CollectAllelicCounts`, which is a pileup.
* PureCN over-segmented when it received the CNVkit segments. It now runs its own
  segmentation.
* The FASTQ index parsing read the wrong field of the read header.

## [0.1.0] - 2026-04-09

The state of the pipeline before this changelog started. It is tagged retroactively at
the tip of `main`.

Mutect2 SNV calling, CNVkit copy number calling, Manta SV calling, Funcotator
annotation, xengsort host read filtering for PDX samples, and an Excel report per tumor
sample.

[2.0.0]: https://github.com/dmkv1/WES-snakemake/compare/v1.3.0...v2.0.0
[1.3.0]: https://github.com/dmkv1/WES-snakemake/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/dmkv1/WES-snakemake/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/dmkv1/WES-snakemake/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/dmkv1/WES-snakemake/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/dmkv1/WES-snakemake/releases/tag/v0.1.0
