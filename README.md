# WES-snakemake: Somatic Variant Calling Pipeline

[![version](https://img.shields.io/github/v/tag/dmkv1/WES-snakemake?label=version&sort=semver)](https://github.com/dmkv1/WES-snakemake/blob/development/CHANGELOG.md)
[![tests](https://github.com/dmkv1/WES-snakemake/actions/workflows/tests.yml/badge.svg)](https://github.com/dmkv1/WES-snakemake/actions/workflows/tests.yml)
![coverage](.github/badges/coverage.svg)

A Snakemake pipeline that calls somatic short variants (SNV and indel), copy number
variants (CNV) and structural variants (SV) from whole-exome sequencing data.
It runs in tumor-normal paired mode and in tumor-only mode with a Panel of Normals.
It also supports Patient-Derived Xenograft (PDX) samples.

## Requirements

The pipeline needs two toolchains, because `profiles/default/config.yaml` sets both
`use-conda` and `use-singularity`:

* A conda solver: conda, mamba or micromamba.
* Apptainer or Singularity.

Create the driver environment from `environment.yml`. For micromamba:

```bash
micromamba env create -f environment.yml
micromamba activate snakemake
```

Snakemake builds all other environments per rule, and pulls the container images on the
first run.

## Reference data

The pipeline does not download reference data. Get each resource before the first run.

| Resource | Config key | Source |
|---|---|---|
| hg38 FASTA, plus `.fai`, `.dict` and the bwa index | `refs.genome_human` | [GATK resource bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0) (`Homo_sapiens_assembly38.fasta`) |
| Mouse mm39 FASTA, for PDX samples | `refs.genome_host` | UCSC or Ensembl |
| BQSR known sites (3 VCFs) | `refs.known_sites` | GATK bucket: `Homo_sapiens_assembly38.dbsnp138.vcf`, `Homo_sapiens_assembly38.known_indels.vcf.gz`, `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` |
| gnomAD germline resource | `refs.germline_resource` | GATK somatic bundle: `af-only-gnomad.hg38.vcf.gz` |
| Common variants for contamination | `refs.contamination_resource` | GATK somatic bundle: `small_exac_common_3.hg38.vcf.gz` |
| VEP offline cache, GRCh38 | `refs.vep_cache.dir` | VEP `INSTALL.pl --CACHEDIR` |
| AnnotSV annotations | `refs.annotsv_annotations.path` | [INSTALL_annotations.sh](https://github.com/lgmgeo/AnnotSV/blob/master/bin/INSTALL_annotations.sh) |
| PureCN SNP blacklist | `refs.purecn.snp_blacklist` | UCSC Simple Repeats, see below |
| somalier sites | `refs.somalier_sites` | somalier releases: `sites.hg38.vcf.gz` |
| DELLY exclude regions | `params.delly.exclude` | DELLY exclude template, see below |
| Capture kit BED files | `probe_configs.<kit>.covered_bedfile` and `.target_regions_bedfile` | The exome kit manufacturer |

The pipeline builds no index. Build the bwa index once with `bwa index`.

**VEP cache version.** `refs.vep_cache.cache_version` must agree with the tag in
`containers.vep`. The repository pins 116 and `ensemblorg/ensembl-vep:release_116.0`.

**DELLY exclude regions.** The file holds telomeres, centromeres and unplaced contigs,
as the DELLY author recommends for short reads. Download the hg38 template:

```bash
wget -O /path/to/refs/delly/human.hg38.excl.tsv \
  https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv
```

The template uses `chr`-prefixed contig names, which agree with
`Homo_sapiens_assembly38.fasta`. An empty value removes the `-x` option from DELLY.

**PureCN SNP blacklist.** PureCN excludes repeat-region variants from the purity and
ploidy fit. Build the BED from the UCSC Tandem Repeats Finder track:

```bash
wget -O simpleRepeat.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
zcat simpleRepeat.txt.gz | cut -f2-4 > hg38_simpleRepeats.bed   # chrom, chromStart, chromEnd, 0-based
```

### Panel of Normals

The sibling pipeline `WES-PON-smk` builds four sets of panel files. Each set has one
entry per capture kit. Two sets also have one entry per sex.

| Artifact | Config key | Used for |
|---|---|---|
| `targets.bed` and `antitargets.bed` | `probe_configs.<kit>.cnvkit_targets` and `.cnvkit_antitargets` | CNVkit coverage |
| `reference_m.cnn` and `reference_f.cnn` | `panel_of_normals.cnvkit.cnvkit_ref_m` and `cnvkit_ref_f` | CNVkit bias correction, all runs |
| `pon.vcf.gz` | `panel_of_normals.mutect2.<kit>` | Mutect2 artifact filtering, tumor-only runs |
| `normalDB_*.rds` and `mapping_bias_*.rds` | `panel_of_normals.purecn.normaldb_{m,f}` and `mapping_bias_{m,f}` | PureCN purity fit, paired runs |

The pipeline stops at startup if a tumor-only run uses a capture kit that has no
Mutect2 PON entry.

## Configuration

`config.yaml` holds the run settings. Repoint paths to actual reference and PON locations.

| Key group | Contents |
|---|---|
| `samplesheet` | Path to the samplesheet CSV, relative to the repository root |
| `containers` | Pinned image tags for GATK, VEP, CNVkit and DELLY |
| `probe_configs.<kit>` | Per capture kit: `covered_bedfile`, `target_regions_bedfile`, `cnvkit_targets`, `cnvkit_antitargets` |
| `refs` | The reference data paths from the table above, plus `refs.path` |
| `params` | Tool settings: `rg.strict`, `xengsort`, `cnvkit.filter_ci`, `cnv.use_purecn_purity`, `fastp`, `bqsr`, `delly`, `somalier` |
| `panel_of_normals` | The panel paths from the table above, plus `panel_of_normals.path` |
| `tumor_only.af_threshold` | The gnomAD allele frequency cutoff for tumor-only runs. Default 0.001 |
| `resources` | Threads, Java heap limits and memory for the scheduler |

**Bind roots.** The Snakefile sets `APPTAINER_BIND` and `SINGULARITY_BIND` from exactly
two values: `refs.path` and `panel_of_normals.path`. Every reference and panel path must
be inside one of these two directories. A path outside them fails inside the container
with a "file not found" message. Do not add bind options to the profile, because the
Snakefile controls the binds.

**Capture kit keys.** The keys of `probe_configs` are the values that the samplesheet
`capture_kit` column uses. To add a kit, add a `probe_configs` block and the matching
`panel_of_normals` entries.

**Notifications.** `telegram_bot_env` points to a file that exports `TELEGRAM_BOT_TOKEN`
and `TELEGRAM_CHAT_ID`. The pipeline then sends a message when the run starts, stops or
fails. `run_id` is a label for those messages. An empty `telegram_bot_env` disables the
notifications.

## Samplesheet

Create a comma-separated table with these columns:

* **ID**: Sample group identifier, usually a patient ID. Samples that share an ID are
  processed together. The value must not contain `_`, `.` or `/`, because output files
  use the `{ID}_{sample}` pattern.
* **sample**: Sample identifier, unique in the given ID.
* **sample_type**: `CTRL` for the matched normal, `PDX` for a xenograft. Any other value
  is a tumor sample.
* **gender**: `f` for female, `m` for male. It selects the sex-matched CNVkit and PureCN
  panel files.
* **capture_kit**: The exome kit version. It must be a key in `probe_configs`.
* **tumor_fraction**: The known tumor cell fraction, in `(0,1]`. Blank or `NA` means
  unknown, and the purity then comes from PureCN. `0` is not valid.
* **fq1**: Path to the R1 FASTQ, or a glob that matches several. See below.
* **fq2**: Path to the R2 FASTQ, which must match `fq1`.

Example:

| ID   | sample   | sample_type | gender | capture_kit | tumor_fraction | fq1 | fq2 |
|------|----------|-------------|--------|-------------|----------------|-----|-----|
| Pt01 | Normal01 | CTRL        | m      | V8+UTR      |                | /path/to/Normal01_R1.fastq.gz | /path/to/Normal01_R2.fastq.gz |
| Pt01 | Tumor01  | Tumor       | m      | V8+UTR      | 0.7            | /path/to/Tumor01_R1.fastq.gz  | /path/to/Tumor01_R2.fastq.gz  |
| Pt01 | Tumor02  | Tumor       | m      | V8+UTR      | 0.3            | /path/to/Tumor02_R1.fastq.gz  | /path/to/Tumor02_R2.fastq.gz  |
| Pt01 | Model01  | PDX         | m      | V8+UTR      | 1              | /path/to/Model01_R1.fastq.gz  | /path/to/Model01_R2.fastq.gz  |

### Multi-lane samples and read groups

The pipeline aligns a sample once per lane. Each alignment gets its own read group.
`MarkDuplicates` then gathers the per-lane BAMs. There are two ways to write this:

* **One row with globbed paths**, for example
  `fq1: /path/SAMPLE_S3_L*_R1_001.fastq.gz`. The glob expands to one *unit* per matched
  pair. This is what the lab ingest emits. A plain path is a glob that matches one file,
  so a single-file row is the one-unit case.
* **One row per FASTQ pair**. Repeat the sample across rows with the same `ID` and
  `sample`. The per-sample columns must agree on every row of a sample. If they
  disagree, the run stops and the message names the conflict.

The pipeline derives the read groups from the data. Four optional columns override this
per row. Add them only when the data cannot supply the value:

* **flowcell**, **lane** and **barcode**. These values are asserted, not derived. A row
  that carries an explicit `lane` must name exactly one FASTQ pair.
* **library**. This is the physical DNA library, and it sets the `LB` tag. It defaults
  to the sample name, which assumes one library. Set it when a sample has two
  independent preparations. `MarkDuplicates` finds duplicates only inside one library.
  If two preparations share one `LB`, it discards real independent evidence as PCR
  duplicates. A resequenced library is not a new library.

`results/metadata/units.tsv` records the result for each unit, together with the exact
`@RG` line. The MultiQC report shows the same data in a **Read groups** table. The
`rg_source` column gives the provenance:

| `rg_source` | Meaning |
|---|---|
| `sheet` | The samplesheet asserts the flowcell and the lane |
| `header+lane` | The flowcell comes from the reads, the lane from the filename. This is the clean per-lane case |
| `header_nolane` | The flowcell is known, the lane is not. **The file can be lane-merged**, so its read group applies to the file, not to a lane |
| `filename` | The read header does not parse. The lane comes from the filename |
| `positional` | There is no provenance. The pipeline numbers the units in order |

The pipeline never takes the lane from a read header. A merged FASTQ can start on one
lane and end on another, so its first record says nothing about the remainder of the
file. The pipeline still reads the header lane. If that lane disagrees with the
filename, the pipeline reports a warning. The disagreement is the signal that a file
spans lanes or that somebody renamed it.

`params.rg.strict` in `config.yaml` controls the response. `true` turns those warnings
into errors, which suits a curated samplesheet. `false` keeps them as warnings, which
suits an unattended run, where a degraded result is better than a stopped one.

## Analysis modes

### Tumor-normal paired mode

The pipeline uses paired mode for a group that has a `CTRL` sample. It compares every
tumor and PDX sample in the group against the single normal.

### Tumor-only mode

The pipeline uses tumor-only mode for a group that has no `CTRL` sample. In this mode:

* Mutect2 uses a Panel of Normals for artifact filtering, in place of a matched normal.
* A population allele frequency filter runs after the PASS filter. It removes variants
  with a gnomAD AF above `tumor_only.af_threshold`.
* CNVkit runs without B-allele frequency data, and HaplotypeCaller does not run.
* PureCN does not run, so the purity comes from the samplesheet or defaults to 1.
* Manta runs in tumor-only mode.
* DELLY calls the tumor alone and applies no somatic filter, because that filter needs a
  matched normal.

Each capture kit in a tumor-only run needs a Mutect2 PON:

```yaml
panel_of_normals:
  path: /path/to/WES-PON-smk
  mutect2:
    "V8+UTR": /path/to/WES-PON-smk/results/PON/mutect2/SureSelectV8UTR/pon.vcf.gz
```

The pipeline validates this at startup.

### PDX samples

A row with `sample_type: PDX` gets an extra host-read filter. The pipeline builds the
xengsort index on the first run that has a PDX sample. The index build is slow, and it
happens once.

## Run the pipeline

Run the commands from the repository root, because `config.yaml`, `profiles/default`
and `multiqc_config.yaml` are relative paths.

```bash
./launch.sh          # full run
./launch.sh -n       # dry run, any extra arguments pass through to snakemake
./stop.sh            # stop a running pipeline
```

`launch.sh` runs in the foreground and writes the log to `snakemake.log`. It records the
Snakemake process ID in `snakemake.pid` and removes that file when the run ends.
`stop.sh` reads `snakemake.pid` and sends `SIGTERM`.

`profiles/default/config.yaml` holds the run settings: 48 cores and a 450 GB memory
limit. Lower both values for a smaller host, and lower `resources.threads` and
`resources.java_max_gb` in `config.yaml` as well. The profile has no retry setting, so a
failed job does not run again.

### Tests

```bash
pytest
```

The suite covers the samplesheet parser and the read-group logic. It is pure Python, so
it runs in seconds and needs no container. Run it after you edit a samplesheet, because
it catches a bad row before a long run starts.

To see the coverage report:

```bash
pytest --cov=workflow/scripts --cov-report=term-missing
```

GitHub Actions runs the same suite on every push to `development` and on every pull
request. It regenerates the coverage badge from the result. The badge covers the Python
helper modules only. It does not measure the R scripts or the Snakemake rules, and an
end-to-end run of the workflow is the only test for those.

## Pipeline steps

This command regenerates the rule graph:

```bash
snakemake --rulegraph --profile profiles/default | dot -Tpng -o rulegraph.png
```

![Pipeline rulegraph](rulegraph.png)

### 1. Quality control and trimming

**fastp** removes adapter sequences and low-quality bases. It runs once per unit. The
HTML and JSON reports go to the MultiQC report.

### 2. Host read filter (PDX samples only)

**xengsort** classifies the trimmed reads against the human (hg38) and mouse (mm39)
genomes. Only the reads that it classifies as "graft", which is human, continue to the
alignment.

### 3. Alignment and BAM preprocessing

The GATK steps use the `broadinstitute/gatk:4.6.2.0` container.

* **bwa mem**, once per unit: `bwa mem -Y -K 100000000 -R '@RG...' | samtools fixmate -m`.
  bwa writes the read group inline, so no later rule adds it. The `-K` option fixes the
  chunk size, so the result does not depend on the thread count. The output is
  query-grouped, not coordinate-sorted.
* **GATK MarkDuplicates** gathers every unit of a sample with `--ASSUME_SORT_ORDER
  queryname`. This is the point where the lanes of a sample join. The metrics go to
  `results/metrics/`.
* **samtools sort** performs the single coordinate sort, after the duplicate marking.
  This order follows GATK Best Practices, so the tool marks supplementary reads
  correctly.
* **GATK BaseRecalibrator and ApplyBQSR** recalibrate the base qualities against dbSNP
  and the known indel sites. BQSR runs for every sample and has no off switch.

`results/{run}/{sample}/bam/{sample}.bam` is the final analysis BAM. All intermediate
BAMs are temporary, and Snakemake deletes them.

### 4. Quality metrics

* **mosdepth** calculates the depth of coverage against the capture BED file.
* **FastQC** reads the trimmed FASTQ files, not the BAM. BQSR changes the base
  qualities, which makes the FastQC chemistry plots meaningless.
* **GATK CollectHsMetrics** reports the capture efficiency, the on-target fraction and
  the fold enrichment.
* **MultiQC** collects fastp, FastQC, mosdepth, duplicate metrics, BQSR tables,
  CollectHsMetrics, the xengsort logs, the somalier tables and two custom tables. The
  custom tables give the read groups and the row type. The report is
  `results/qc/multiqc_report.html`.

The MultiQC report uses three row levels. FastQC rows are per unit and per read. fastp
rows are per unit. Every metric from the finished BAM has one row per sample.

### 5. SNV and indel calling

* **GATK Mutect2** calls somatic variants inside the capture regions.
  * *Paired mode*: tumor and matched normal BAMs, with the gnomAD germline resource.
  * *Tumor-only mode*: the PON and the gnomAD germline resource, with no normal BAM.
* **GATK LearnReadOrientationModel** builds the orientation-bias priors from the F1R2
  counts. These priors cover FFPE and OxoG artifacts.
* **GATK GetPileupSummaries and CalculateContamination** estimate the cross-sample
  contamination. In paired mode, the calculation uses the matched normal.
* **GATK FilterMutectCalls** applies the statistical filter model, with the orientation
  priors and the contamination estimate.
* **bcftools** keeps the PASS variants and sorts them.
* **Population AF filter** (tumor-only only) removes the variants with a gnomAD AF above
  `tumor_only.af_threshold`.
* **VEP** annotates the variants offline from the local cache, with `--everything` and
  `--pick`. The container is `ensemblorg/ensembl-vep:release_116.0`.

### 6. Copy number variant calling

**CNVkit** runs in the `etal/cnvkit:0.9.14` container.

* CNVkit calculates the coverage at the target and antitarget regions. It then corrects
  the coverage against the sex-matched panel reference for the capture kit.
* The pipeline keeps the standard chromosomes only: 1 to 22, X, Y and M.
* *Paired mode*: **GATK HaplotypeCaller** calls the germline variants of the normal
  only. The pipeline extracts the heterozygous SNPs of the normal. **GATK
  CollectAllelicCounts** then counts the alleles of the tumor at those same sites. Both
  results merge into one B-allele frequency VCF, which informs the CBS segmentation.
* *Tumor-only mode*: the segmentation runs without B-allele frequency data.
* `cnvkit segmetrics` calculates the segment confidence intervals by bootstrap.
* An optional confidence-interval filter runs if you set `params.cnvkit.filter_ci: true`.
* `cnvkit segmetrics --t-test` adds the t-test statistics.
* `cnvkit call` makes the copy number calls with the clonal model and median centering.
  It reads the purity and the ploidy from the sidecar file of step 7.
* CNVkit draws a scatter plot (PNG) and a chromosome diagram (PDF) per sample.

### 7. Purity and ploidy

**PureCN** estimates the tumor purity and the ploidy. It runs for the tumors of a paired
run only, because it needs the B-allele frequency data and the panel normal database.
PureCN performs its own CBS segmentation on the reformatted CNVkit coverage.

One rule resolves the purity, and it uses this order:

1. The samplesheet `tumor_fraction`, if the value is known.
2. The PureCN estimate, if `params.cnv.use_purecn_purity` is true, PureCN ran, PureCN
   did not report `Fail`, and its only flag is `POOR GOF`.
3. Purity 1 and ploidy 2.

The result goes to `work/purity/{run}/{sample}/{sample}.purity.csv`. Both `cnvkit call`
and the results report read this one file, so the two cannot disagree. PureCN runs for
every paired tumor even when `params.cnv.use_purecn_purity` is false, so its numbers
stay available for comparison.

### 8. Sample identity

**somalier** extracts a genotype sketch per sample. It then relates the whole cohort,
every sample against every other sample. A cohort-wide comparison makes a sample swap
between patients visible. An R script adds the expected relationship and a
PASS/WARN/FAIL verdict, and draws a heatmap. The thresholds are in `params.somalier`.

This step reports only. A FAIL verdict does not stop the pipeline.

### 9. Structural variant calling

Two independent callers run, and one merge step reconciles them.

* **Manta** detects SVs inside the capture regions. Manta rejects high base qualities,
  so it reads a copy of the BAM with the base qualities capped at 70. No other caller
  uses that copy.
  * *Paired mode*: somatic SVs. *Tumor-only mode*: tumor SVs. The pipeline converts the
    paired breakend records into inversions and renames the output to one consistent
    name, so the later rules need no branch.
* **DELLY** is the second caller. It reads the final analysis BAM directly, with no
  base-quality cap. It runs with the DELLY author's short-read settings, which
  `params.delly` controls. These are the exclude regions (`-x`), the minimum mapping
  quality `-q 20`, the insert-size MAD cutoff `-s 15` and the minimum clique size
  `-z 5`. DELLY has no exome mode, so the pipeline restricts its calls to the capture
  BED after the call step.
  * *Paired mode*: `delly sr` on the tumor and the matched normal, then
    `delly filter -f somatic`. *Tumor-only mode*: `delly sr` on the tumor alone, with no
    somatic filter.
* **bcftools** keeps the PASS variants of both callers and sorts them.
* **SURVIVOR** merges the two PASS VCF files into one consensus set. It records the
  caller support in `SUPP_VEC`, where bit 1 is Manta and bit 2 is DELLY. The breakpoint
  distance and the minimum size are in `params.delly`.
* A size filter removes the intra-chromosomal SVs longer than `params.delly.max_sv_size`.
  It always keeps the breakend records.
* **AnnotSV** annotates the consensus set against the GRCh38 database. If no variant
  passes the filters, the pipeline writes an empty file.

### 10. Results integration

`combine_results.R` merges the results of one tumor sample into one Excel workbook.

**Output: `{sample}_results.xlsx`**, four sheets:

| Sheet | Contents |
|---|---|
| **SNVs** | VEP-annotated PASS SNVs and indels, with AF, GT, DP and AD per sample, the copy number context from CNVkit (total CN, major and minor allele CN) and the Cancer Cell Fraction |
| **CNVs** | The CNVkit segment-level copy number calls |
| **SVs** | The AnnotSV-annotated consensus SVs, with an `SV_callers` column that names the supporting callers |
| **Sample_QC** | The purity, the ploidy, the purity source, and the PureCN estimate with its flags |

**Cancer Cell Fraction.** The CCF estimates the fraction of cancer cells that carry the
mutation:

```
CCF = (AF × tumor_CN) / (purity × expected_mutant_copies)
```

`expected_mutant_copies` is `tumor_CN` for a homozygous (1/1) call, and 1 for a
heterozygous call. The CCF has a cap of 1.0. The purity comes from the sidecar file of
step 7. Variants outside a CNVkit segment get the normal copy number: 2 for the
autosomes, and 1 for chrX and chrY in a male sample.

Two more steps run after the per-sample workbooks:

* `combine_all.R` joins every sample into four cohort tables in `results/combined/`.
* `combined_to_maf.R` converts the cohort SNV table into `results/combined/cohort.maf`.
  The MAF step is separate, so a failure there leaves the cohort tables intact.

## Outputs

`results/` holds the deliverables.

| Path | Description |
|---|---|
| `results/{run}/{sample}/{sample}_results.xlsx` | Combined SNV, CNV and SV workbook |
| `results/{run}/{sample}/{sample}.SNV.vcf` | VEP-annotated SNV and indel VCF |
| `results/{run}/{sample}/{sample}.scatter.png` | CNVkit genome-wide scatter plot |
| `results/{run}/{sample}/{sample}.diagram.pdf` | CNVkit per-chromosome diagram |
| `results/{run}/{sample}/bam/{sample}.bam` and `.bai` | Final recalibrated analysis BAM |
| `results/combined/combined_snvs.tsv` | Cohort SNV table |
| `results/combined/combined_cnvs.tsv` | Cohort CNV table |
| `results/combined/combined_svs.tsv` | Cohort SV table |
| `results/combined/combined_qc.tsv` | Cohort QC table, with the purity and the somalier verdict |
| `results/combined/cohort.maf` | Cohort MAF, for maftools and similar tools |
| `results/combined/combined_relatedness.tsv` | somalier pair table with the verdicts |
| `results/combined/relatedness_heatmap.png` | somalier relatedness heatmap |
| `results/metadata/units.tsv` | One row per alignment unit, with the exact `@RG` line |
| `results/metadata/samples.tsv` | One validated row per sample |
| `results/qc/multiqc_report.html` | Aggregated QC report |
| `results/qc/somalier/cohort/cohort.html` | somalier interactive report |
| `results/qc/fastp/{run}/{sample}.{unit}_fastp.html` | fastp report per unit |
| `results/qc/fastqc/{run}/{sample}.{unit}_R{1,2}_fastqc.html` | FastQC report per unit and read |
| `results/metrics/dupl_metrics_{run}_{sample}.txt` | MarkDuplicates metrics |
| `results/metrics/{run}_{sample}.recal_data.table` | BQSR recalibration table |
| `results/metrics/{run}_{sample}.hs_metrics.txt` | CollectHsMetrics capture metrics |
| `results/metrics/mutect2_filteringStats_{run}_{sample}.tsv` | Mutect2 filter statistics |
| `results/metrics/{run}/{sample}.mosdepth.summary.txt` | Coverage summary by region |
| `results/metrics/{run}/{sample}.thresholds.bed.gz` | Coverage thresholds (10x, 20x, 30x, 50x) |

`work/` holds the intermediate files. Snakemake deletes the temporary files
automatically. Three groups of results stay in `work/` permanently, and they reach you
through the workbook and the cohort tables:

* `work/cnvkit/{run}/{sample}/{sample}.call.cns`, the CNVkit segment calls.
* `work/sv_merge/{run}/{sample}/{sample}.SV.annotated.tsv`, the AnnotSV output.
* `work/purecn/{run}/{sample}/`, the PureCN report, and
  `work/purity/{run}/{sample}/{sample}.purity.csv`, the resolved purity.

## Versioning

The pipeline uses [Semantic Versioning](https://semver.org). A git tag is the source of
truth, and the `VERSION` file holds the same number. `CHANGELOG.md` records every
release.

The public contract is the samplesheet schema, the `config.yaml` keys and the output
paths. A version part changes for these reasons:

| Part | Reason |
|---|---|
| MAJOR | A samplesheet column, a `config.yaml` key or an output path changes. A reference file or a Panel of Normals artifact becomes necessary. An output is removed |
| MINOR | A tool, a rule or an output is added, and an existing configuration still runs |
| PATCH | A bug fix, a documentation change, or a pinned version bump that does not change results |

A pipeline has one rule that a software library does not need. **Any release that changes
results for the same input data must say so at the top of its changelog entry.** A tool
that produces different numbers is a scientific event, even when no key and no path
changes. Do not compare a cohort across such a release. Run the old data again.

To make a release:

1. Update `CHANGELOG.md` and `VERSION`.
2. Commit, then tag: `git tag -a v1.0.0 -m "v1.0.0"`.
3. Push both: `git push origin development --follow-tags`.

The version badge reads the newest tag, so it updates by itself.

## Known limitations

**The tumor-only population filter can remove true drivers.** In tumor-only mode, the
pipeline removes the variants with a gnomAD AF above `tumor_only.af_threshold`, which
defaults to 0.001. Some recurrent somatic mutations exceed that frequency in gnomAD, for
example DNMT3A variants and JAK2 V617F. The filter therefore removes a small number of
true somatic drivers. Read the tumor-only SNV results with this in mind. A paired run
does not use this filter.

**The preprocessing differs from the panel baseline.** `WES-PON-smk` built the panels of
normals with an older preprocessing chain. That chain used AddOrReplaceReadGroups,
FixMateInformation and coordinate-sorted duplicate marking. This pipeline uses per-unit
read groups, `samtools fixmate` and query-grouped duplicate marking. The two chains mark
supplementary reads differently, so the coverage is not identical to the panel baseline.
The lab accepts this difference. `TODO.md` item 27 tracks the correction, which needs a
rebuild of both panels.

## Acknowledgments

This pipeline uses the following open-source tools:

#### GATK (Broad Institute)

```
Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
```

#### Mutect2

```
Cibulskis, K., Lawrence, M., Carter, S. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnol 31, 213-219 (2013). https://doi.org/10.1038/nbt.2514
```

#### Ensembl VEP

```
McLaren W, Gil L, Hunt SE et al. The Ensembl Variant Effect Predictor. Genome Biology 17:122 (2016). https://doi.org/10.1186/s13059-016-0974-4
```

#### CNVkit

```
Talevich E, Shain AH, Botton T, Bastian BC. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLOS Computational Biology 12(4) (2016). https://doi.org/10.1371/journal.pcbi.1004873
```

#### PureCN

```
Riester M, Singh AP, Brannon AR et al. PureCN: copy number calling and SNV classification using targeted short read sequencing. Source Code for Biology and Medicine 11:13 (2016). https://doi.org/10.1186/s13029-016-0060-z
```

#### Manta

```
Chen X et al. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics 32(8):1220-1222 (2016). https://doi.org/10.1093/bioinformatics/btv710
```

#### DELLY

```
Rausch T et al. DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics 28(18):i333-i339 (2012). https://doi.org/10.1093/bioinformatics/bts378
```

#### SURVIVOR

```
Jeffares DC et al. Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast. Nature Communications 8:14061 (2017). https://doi.org/10.1038/ncomms14061
```

#### AnnotSV

```
Geoffroy V et al. AnnotSV: an integrated tool for structural variations annotation. Bioinformatics 34(20):3572-3574 (2018). https://doi.org/10.1093/bioinformatics/bty304
```

#### somalier

```
Pedersen BS, Bhetariya PJ, Brown J et al. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches. Genome Medicine 12:62 (2020). https://doi.org/10.1186/s13073-020-00761-2
```

#### Xengsort (TU Dortmund University)

```
Zentgraf, J., Rahmann, S. Fast lightweight accurate xenograft sorting. Algorithms Mol Biol 16, 2 (2021). https://doi.org/10.1186/s13015-021-00181-w
```

#### BWA

```
Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler Aligner. Bioinformatics 25(14):1754-1760 (2009). https://doi.org/10.1093/bioinformatics/btp324
```

#### SAMtools and BCFtools

```
Danecek P, Bonfield JK, Liddle J et al. Twelve years of SAMtools and BCFtools. GigaScience 10(2):giab008 (2021). https://doi.org/10.1093/gigascience/giab008
```

#### fastp

```
Chen S et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17):i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560
```

#### FastQC

```
Andrews S. FastQC: a quality control tool for high throughput sequence data (2010). https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
```

#### mosdepth

```
Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics 34(5):867-868 (2018). https://doi.org/10.1093/bioinformatics/btx699
```

#### MultiQC

```
Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19):3047-3048 (2016). https://doi.org/10.1093/bioinformatics/btw354
```

#### Snakemake

```
Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 3; peer review: 2 approved]. F1000Research 2025, 10:33 (https://doi.org/10.12688/f1000research.29032.3)
```
