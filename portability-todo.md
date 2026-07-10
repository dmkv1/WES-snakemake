# WES-snakemake Portability Plan

Goal: a consumer with Docker/Apptainer (no conda/mamba, no `/mnt/data/NGS/...`
layout) can point the pipeline at their own reference tree and samplesheet and
run it. Two independent problems, tackled in order: **(1) hardcoded paths**
(blocks everyone regardless of container strategy) and **(2) mixed
conda+container deployment** (blocks anyone without a conda/mamba toolchain
and undermines the reproducibility the pinned containers are supposed to
give).

Reference/PON data itself (hundreds of GB) is never bundled — always an
external mount. Not in scope here.

---

## 1. Config/path portability (do first, independent of containers)

1. **Split `config.yaml` into template + local override.**
   - Add `config.yaml.example` (committed) with placeholder paths and inline
     comments describing the expected directory layout.
   - Add `config.yaml` to `.gitignore` (verify it isn't already tracked with
     real paths committed — it currently is, per `git status`; decide whether
     to `git rm --cached` it or just stop tracking future edits).
   - Document in README: copy the example, edit paths, done.
   - **Effort:** low.

2. **Consolidate reference layout under one configurable root.** Right now
   paths are hardcoded piecemeal: `refs.path`, `refs.genome_human`,
   `refs.genome_host`, `refs.known_sites[]`, `refs.germline_resource`,
   `refs.contamination_resource`, `refs.vep_cache.dir`,
   `refs.annotsv_annotations.path`, `refs.somalier_sites`, plus
   `probe_configs.*.covered_bedfile`/`target_regions_bedfile`, plus
   `panel_of_normals.path` and its per-kit `cnvkit_ref_m/f`/`mutect2` maps.
   All of these currently happen to live under `/mnt/data/NGS/...` but nothing
   enforces that — a consumer can easily end up with a bind mount that misses
   one of them, surfacing as an opaque "file not found" *inside the
   container* rather than a clear error at startup.
   - Keep the existing key structure (it's fine — don't force a rigid
     subdirectory convention that fights how consumers already have refs
     laid out), but require every path in `config.yaml` to resolve under one
     of a small, explicit list of root directories declared once, e.g.:
     ```yaml
     bind_roots:
       - /path/to/refs
       - /path/to/pon
     ```
   - `Snakefile`'s `SINGULARITY_BIND` construction (currently
     `[refs.path, panel_of_normals.path]`, `Snakefile:8-10`) should build from
     `bind_roots` instead, so adding a new external resource (e.g. a
     consumer's refs live across two unrelated mounts) is a one-line config
     change, not a code change.
   - **Effort:** low-medium.

3. **Startup path validation, extended from the existing fastq check.**
   `Snakefile:34-44` already validates fastq paths exist and fails fast with
   a clear message before Snakemake builds the DAG. Extend the same pattern
   to every reference/PON/BED path pulled from `config.yaml` (genome fastas +
   `.dict`/`.fai`, known-sites VCFs, germline/contamination resources, VEP
   cache dir, AnnotSV annotations dir, somalier sites, per-probe BED files,
   PON files for probe kits actually used in the samplesheet). One
   consolidated error listing everything missing, not a chain of
   trial-and-error re-runs.
   - **Effort:** low. Mechanical extension of code that already exists.

4. **Fix the `SINGULARITY_BIND` vs. README double-documentation.** README
   currently tells the user to *also* hand-edit `singularity-args:
   "-B /path:/path"` in the profile (`README.md:110-116`), duplicating what
   `Snakefile:8-10` already does automatically via `SINGULARITY_BIND`. Pick
   one:
   - **Recommended:** keep the automatic env-var approach (item 2 above,
     driven by `bind_roots`), delete the manual-bind instructions from the
     README entirely. One mechanism, not two that can silently disagree.
   - **Effort:** low (mostly deleting stale docs once item 2 lands).

5. **Gate the Telegram webhook harder.** `_tg_notify` (`Snakefile:139-156`)
   is already config-gated (`telegram_bot_env` empty ⇒ no-op) so it's not a
   blocker, but a fresh consumer's `config.yaml.example` should ship with it
   commented out / empty by default rather than relying on them noticing.
   - **Effort:** trivial, bundle with item 1.

---

## 2. Eliminate the conda/container split

Currently `use-conda: true` **and** `use-singularity: true`
(`profiles/default/config.yaml`) — a consumer needs both a conda/mamba
solver *and* Apptainer/Singularity working. Converting every remaining
`conda:` rule to `container:` drops the requirement to "Apptainer or Docker
only" and closes the reproducibility gap between the 3 version-pinned
containers (GATK/CNVkit/VEP) and the conda envs, which pin package version
but not build string and can drift on re-solve.

Per-rule mapping below. "Direct" = an existing single-tool biocontainers
image is a drop-in replacement, just change `conda:` → `container:`. "Custom"
= no single upstream image covers the rule's tool combination; build a small
project-maintained Dockerfile instead. Tags below were checked live against
`quay.io/biocontainers` — pin exact tags in `config["containers"]` the same
way GATK/CNVkit/VEP already are, don't float `:latest`.

6. **`fastp_trim`** (`bam_mapping_gatk.smk:1`) — `workflow/envs/fastp.yaml` →
   **Direct**: `quay.io/biocontainers/fastp:0.23.4--h125f33a_5` (verify exact
   build tag at swap time). Effort: trivial.

7. **`bwa_map`** (`bam_mapping_gatk.smk:29`) — pipes `bwa mem | samtools
   view` in one shell command, needs both tools in one container.
   `workflow/envs/bwamem.yaml` has no single-tool upstream equivalent.
   - **Option A (preferred): split the rule** into `bwa_mem` (SAM output,
     `container: biocontainers/bwa:0.7.19--h577a1d6_1`) → `sam_to_bam`
     (`container: biocontainers/samtools:1.21--...`). Costs one intermediate
     SAM on disk (`temp()`, cleaned up automatically) per sample; simplest,
     no custom image to maintain.
   - **Option B: custom image** — tiny Dockerfile installing
     `bwa=0.7.19 samtools=1.21` (mirrors `bwamem.yaml` exactly), built and
     pushed once. Keeps the current single-rule pipe (no intermediate SAM),
     more maintenance (rebuild/repush on version bump).
   - Also used by `ref_index.smk`'s `index_bed`/`get_chromosomes` — but
     those only need `bgzip`/`tabix`/`grep`/`sort`/`cut`, not bwa or samtools
     at all. Point them at **`quay.io/biocontainers/htslib:1.24--ha79157c_0`**
     instead (the actually-correct minimal tool for those two rules) rather
     than dragging in the bwa/samtools image.
   - **Effort:** low (Option A) or low-medium (Option B, one-time image build).

8. **`mosdepth`, `fastqc`, `multiqc`** (`bam_mapping_gatk.smk:217`,
   `qc.smk:53`, `qc.smk:102`) — currently share one fat `workflow/envs/qc.yaml`
   (mosdepth+fastqc+multiqc) even though each rule only uses one of the three
   tools. **Direct**, one image per rule:
   - `quay.io/biocontainers/mosdepth:0.3.11--h0ec343a_0` (or `_1`, match
     pinned version exactly)
   - `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`
   - `quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0`
   - Delete `qc.yaml` once all three rules are converted.
   - **Effort:** low.

9. **`somalier_extract`, `somalier_relate`** (`somalier.smk:5,46`) —
   **Direct**, but verify carefully:
   `quay.io/biocontainers/somalier:0.3.2--h5205c93_0`. TODO.md #20 notes the
   *vendor's* `brentp/somalier` Docker Hub images (v0.3.2/v0.3.2-1) ship a
   broken `/bin/bash` and that's why this rule is on conda today. The
   biocontainers build is a separate build pipeline (bioconda recipe →
   BioContainers CI, not brentp's own Dockerfile) so it may well not share
   the bug — but **must be smoke-tested before switching**, don't assume.
   If it's also broken, keep this one rule on conda as a documented
   exception (or build a custom minimal image).
   - **Effort:** low, plus one verification run.

10. **`combine_results`, `merge_results`** (`integrate_results.smk:1,21`) —
    `workflow/envs/r_vcf.yaml` (r-base=4.3 + bioconductor-variantannotation +
    r-tidyverse + r-openxlsx). No single upstream image bundles this
    combination. **Custom**: small Dockerfile,
    `FROM bioconductor/bioconductor_docker:RELEASE_3_18` (or the closest
    matching Bioconductor release for R 4.3) + `install.packages`/
    `BiocManager::install` for the four packages. Maintained in-repo
    (`containers/r-vcf/Dockerfile`), built once, pinned by tag.
    - **Effort:** medium (one-time image build + pin), low ongoing.

11. **`build_xengsort_index`, `run_xengsort`** (`host_read_filter.smk:1,28`)
    — `workflow/envs/xengsort.yaml` → **Direct**:
    `quay.io/biocontainers/xengsort:2.0.8--pyhdfd78af_0`.
    - **Effort:** low.

12. **`run_manta`** (`sv_calling_manta.smk:1`) — `workflow/envs/manta.yaml`
    (python2.7 + manta=1.6.0) → **Direct**:
    `quay.io/biocontainers/manta:1.6.0--py27h9948957_6` (or `_5`/`_4`, pin
    to whatever build matches the currently-solved env if reproducibility
    of existing results matters — check `.snakemake/conda/*.yaml` for the
    exact resolved build before picking).
    - **Effort:** low.

13. **`filter_manta_variants`** (`sv_calling_manta.smk:53`) — 
    `workflow/envs/bcftools.yaml` → **Direct**:
    `quay.io/biocontainers/bcftools:1.21--h8b25389_0` (or `_1`).
    - **Effort:** trivial.

14. **`annotate_sv`** (`sv_calling_manta.smk:67`) —
    `workflow/envs/annotsv.yaml` (annotsv=3.5.3 + bcftools=1.21 — bcftools is
    only used for a variant-count check to skip AnnotSV on empty VCFs).
    Cheapest fix: **replace the `bcftools view -H | wc -l` variant-count
    check with a pure-Python line count** (skip `#`-prefixed lines) inside
    the same `shell:` block or move it to a tiny `run:` block — removes the
    bcftools dependency from this rule entirely. Then it's **Direct**:
    `quay.io/biocontainers/annotsv:3.5.3--py313hdfd78af_0` (match pinned
    version).
    - **Effort:** low.

15. **`purecn_run`** (`purecn.smk:31`) — currently dead code, not `include`d
    in `Snakefile` (TODO.md #25). No action needed until #25 is resolved;
    when it is, `workflow/envs/purecn.yaml` (bioconductor-purecn=2.12.0 +
    r-optparse) → **Direct**-ish:
    `quay.io/biocontainers/bioconductor-purecn:2.12.0--r44hdfd78af_0`, but
    the shell command hardcodes `$CONDA_PREFIX/lib/R/library/PureCN/extdata/PureCN.R`
    (`purecn.smk:54`) — that path doesn't exist in the container (no
    `$CONDA_PREFIX`); needs updating to wherever the biocontainer installs R
    libraries (typically `/usr/local/lib/R/library/...` — verify at
    conversion time).
    - **Effort:** low, parked behind #25.

16. **Drop `environment.yml` / micromamba as the way to get Snakemake
    itself**, once items 6-15 land and `use-conda: true` is no longer
    needed anywhere. Replace with either:
    - a `container:` on the top-level driver only needing Apptainer/Docker
      installed natively, or
    - keep `environment.yml` as one *documented* option for users who prefer
      a native conda install, but make it optional, not required.
    - Flip `profiles/default/config.yaml` to `use-conda: false`.
    - **Effort:** low, last step once everything above is converted.

---

## 3. Driver packaging (thin, not monolithic)

17. **Ship a thin driver image**, not a fat one containing GATK+CNVkit+VEP+
    everything baked in. The per-rule `container:` model already externalizes
    the heavy tools into their own pinned images — collapsing them into one
    15-20GB blob would fight that and gain nothing.
    - `containers/driver/Dockerfile`: minimal base + Snakemake + pandas
      (mirrors `environment.yml`'s pins) + Apptainer client (or documented as
      "run natively with Apptainer installed on host, this image only
      wraps Snakemake+pandas").
    - Consumer flow becomes: `docker run` (or native Snakemake) driver
      → driver invokes `snakemake --sdm apptainer` → Snakemake pulls each
      rule's pinned tool image on demand. No conda solve, ever.
    - **Effort:** medium.

18. **`--software-deployment-method` note (Snakemake ≥8 syntax).** Current
    profile uses the older `use-conda`/`use-singularity` flags
    (`profiles/default/config.yaml:2-3`); `environment.yml` already pins
    `snakemake>=8.29.3`, so once item 16 lands, migrate the profile to
    `software-deployment-method: apptainer` (the flags still work but the
    newer syntax is the one Snakemake 8+ documents going forward).
    - **Effort:** trivial, bundle with item 16.

---

## Suggested order

Do **section 1 in full first** — it's required no matter what container
strategy is chosen, and is what actually blocks a consumer today. Section 2
items are independent of each other and can land one rule-group at a time
(each is a small, verifiable diff — convert, dry-run, run on `P005` end-to-end,
compare output byte-for-byte against the conda version before moving to the
next). Section 3 last, once nothing depends on conda.

## Out of scope

- Bundling reference/PON data — always external, documented as a directory
  contract, never shipped in an image.
- Rewriting the per-rule container model into a single monolithic image —
  actively worse for maintenance and image size, rejected above (item 17).
