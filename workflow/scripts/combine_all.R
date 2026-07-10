#!/usr/bin/env Rscript
# combine_all.R -- row-bind every sample's per-type CSVs into three
# cross-cohort TSVs (SNVs, CNVs, SVs). Called by snakemake rule merge_results.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

samplesheet <- snakemake@params[["samplesheet"]]

# Load samplesheet for metadata (drop FASTQ paths — not needed). Rename
# 'probes' up front: CNVkit's per-segment CSV also has a 'probes' column
# (probe count), so joining unrenamed would collide.
meta <- read_csv(samplesheet, show_col_types = FALSE) %>%
  select(-any_of(c("fq1", "fq2"))) %>%
  rename(sID = sample, capture_kit = probes)

# AnnotSV places VCF sample columns (named after the actual sample IDs,
# which differ per row) between FORMAT (pos 14) and Annotation_mode.
# Rename them to generic names before binding.
annotsv_fixed_cols <- c(
  "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length",
  "SV_type", "Samples_ID", "ID", "REF", "ALT", "QUAL",
  "FILTER", "INFO", "FORMAT"
)

rename_sv_sample_cols <- function(df) {
  annot_mode_pos <- match("Annotation_mode", names(df))
  if (is.na(annot_mode_pos)) return(df)
  vcf_sample_cols <- setdiff(names(df)[seq_len(annot_mode_pos - 1L)], annotsv_fixed_cols)
  for (i in seq_along(vcf_sample_cols)) {
    new_name <- if (i == 1L) "FORMAT_vcf" else paste0("FORMAT_vcf_", i)
    names(df)[names(df) == vcf_sample_cols[i]] <- new_name
  }
  df
}

# Row-bind one data type's per-sample CSVs, tag with sample/run, and join
# samplesheet metadata (moved to the front of the result).
combine_csvs <- function(csv_files, rename_fn = identity) {
  parts <- lapply(csv_files, function(f) {
    # Path shape: results/{run}/{sample}/{sample}_{TYPE}.csv
    sample_name <- basename(dirname(f))
    run_name    <- basename(dirname(dirname(f)))

    df <- tryCatch(
      read_csv(f, col_types = cols(.default = col_character())),
      error = function(e) {
        message("Warning: could not read ", f, ": ", e$message)
        NULL
      }
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df <- rename_fn(df)
    df$sID <- sample_name
    df$ID  <- run_name
    df
  })

  combined <- bind_rows(parts)
  combined <- type_convert(combined, col_types = cols())
  combined <- left_join(combined, meta, by = c("sID", "ID"))

  meta_cols <- intersect(names(meta), names(combined))
  other_cols <- setdiff(names(combined), meta_cols)
  combined[, c(meta_cols, other_cols)]
}

combined_snvs <- combine_csvs(snakemake@input[["snv_csvs"]])
combined_cnvs <- combine_csvs(snakemake@input[["cnv_csvs"]])
combined_svs  <- combine_csvs(snakemake@input[["sv_csvs"]], rename_sv_sample_cols)

write_tsv(combined_snvs, snakemake@output[["snv_tsv"]])
write_tsv(combined_cnvs, snakemake@output[["cnv_tsv"]])
write_tsv(combined_svs,  snakemake@output[["sv_tsv"]])

message(sprintf(
  "[OK] Combined %d SNV / %d CNV / %d SV rows across %d samples",
  nrow(combined_snvs), nrow(combined_cnvs), nrow(combined_svs),
  length(unique(combined_snvs$sID))
))
