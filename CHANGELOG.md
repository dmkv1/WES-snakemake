# Changelog

This project uses [Semantic Versioning](https://semver.org), adapted for a pipeline.
The public contract is the samplesheet schema, the `config.yaml` keys and the output
paths. See [Versioning](README.md#versioning) in the README.

Each release states whether it changes results for the same input data. A version that
changes results needs a re-run before you compare old and new cohorts.

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

[1.0.0]: https://github.com/dmkv1/WES-snakemake/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/dmkv1/WES-snakemake/releases/tag/v0.1.0
