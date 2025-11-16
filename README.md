# WES-snakemake: Somatic Variant Calling Pipeline

A Snakemake-based pipeline for calling short somatic nucleotide variants (SNV) in tumor-normal pairs, with support for Patient-Derived Xenograft (PDX) samples.

## Installation and Dependencies

Use your python environment manager (conda, etc.) to install environment from `environment.yml `. E.g., for micromamba:

```bash
micromamba env create -f environment.yml
```

Then activate `snakemake-wes` environment.

### Configuration

Modify `config.yaml` to specify:

* Paths to reference files
* Computational resources (threads, memory)

### Reference files

Resources which can be downloaded from the [GATK resource bucket](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0):
* Human reference genome hg38 in fasta format (`Homo_sapiens_assembly38.fasta`)
* Known sites for base recalibration - can be downloaded from , files `Homo_sapiens_assembly38.dbsnp138.vcf`, `Homo_sapiens_assembly38.known_indels.vcf.gz` and `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
* Germline resource for Mutect2 - `af-only-gnomad.hg38.vcf.gz`
* Funcotator data sources - [`funcotator_dataSources.v1.8.hg38.20230908s`](https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s)
* Host reference genome in fasta format for PDX samples (mm39)

Other resources:
* Exome capture regions in BED format should be provided by the exome library preparation kit manufacturer.
* AnnotSV annotations could be downloaded using [INSTALL_annotations shell script](https://github.com/lgmgeo/AnnotSV/blob/master/bin/INSTALL_annotations.sh).

## Running the Pipeline

### Samplesheet Preparation

Create comma-separated table `samplesheet.csv` with columns:

* **ID**: Sample group identifier for grouping samples from the same origin. Usually a patient ID.
* **sample**: Sample identifier unique for the given ID.
* **sample_type**: Sample type - `CTRL` for normal samples, `PDX` for xenograft samples, any other value would be interpreted as a tumor sample.
* **Chr_sex**: Chromosomal sex of each sample, `XX` for female and `XY` for male.
* **probes**: Exome library kit version.
* **purity**: Known cancer/normal cell ratio, from 0 to 1. Used in CNVkit `call` and CCF calculation.
* **fq1**: Full path to R1 FASTQ file
* **fq2**: Full path to R2 FASTQ file

Example format:

| ID   |   sample   | sample_type | Chr_sex | probes |  purity |  fq1 | fq2 |
|------|------------|-------------|---------|--------|---------|------|-----|
| Pt01 |  Normal01  |     CTRL    |    XY   | V8+UTR |    0    | /path/to/Normal01_R1.fastq.gz | /path/to/Normal01_R2.fastq.gz |
| Pt01 |  Tumor01   |     Tumor   |    XY   | V8+UTR |    0.7  | /path/to/Tumor01_R1.fastq.gz  | /path/to/Tumor01_R2.fastq.gz  |
| Pt01 |  Tumor02   |     Tumor   |    XY   | V8+UTR |    0.3  | /path/to/Tumor02_R1.fastq.gz  | /path/to/Tumor02_R2.fastq.gz  |
| Pt01 |  Model01   |     PDX     |    XY   | V8+UTR |    1    | /path/to/Model01_R1.fastq.gz  | /path/to/Model01_R2.fastq.gz  |

Each run must have exactly one normal/germline sample to which all other samples would be compared.

### Launching the Pipeline

Before launching the pipeline, check the profile config in `profiles/default/config.yaml`. The workflow requires both conda and singularity.

**For containers to work you have to bind the reference directory!** I.e., all required reference sources should be put in this directory. Edit the `singularity-args` parameter in the profile config:

```yaml
singularity-args: "-B /path/to/refs:/path/to/refs"
```

Use `snakemake --profile profiles/default -n` to test the pipeline ("dry run").

Use `launch.sh` which would attempt to launch snakemake in the detached mode in the background, or run in the foreground using `snakemake --profile profiles/default`.

`launch.sh` records the process ID which could be used to stop the execution using `stop.sh`.

## Outputs

TBD

## Acknowledgments

This pipeline uses several open-source tools:

#### GATK (Broad Institute)

```
Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
```

#### Mutect

```
Cibulskis, K., Lawrence, M., Carter, S. et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnol 31, 213–219 (2013). https://doi.org/10.1038/nbt.2514
```

#### Xengsort (TU Dortmund University)

```
Zentgraf, J., Rahmann, S. Fast lightweight accurate xenograft sorting. Algorithms Mol Biol 16, 2 (2021). https://doi.org/10.1186/s13015-021-00181-w
```

#### Snakemake

```
Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.2) 
```
