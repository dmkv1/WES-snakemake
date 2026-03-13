#!/usr/bin/env Rscript
# combine_svs.R — Row-bind all per-sample SV CSVs into one combined RDS
# Called by snakemake rule combine_all_svs
#
# AnnotSV places VCF sample columns between FORMAT (pos 14) and
# Annotation_mode, named after the actual sample IDs. These must be
# renamed to generic names (FORMAT_vcf, FORMAT_vcf_2) before binding.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

csv_files   <- snakemake@input[["sv_csvs"]]
samplesheet <- snakemake@params[["samplesheet"]]
output_rds  <- snakemake@output[["rds"]]

# Load samplesheet for metadata (drop FASTQ paths)
meta <- read_csv(samplesheet, show_col_types = FALSE) %>%
  select(-any_of(c("fq1", "fq2"))) %>%
  rename(sID = sample, capture_kit = probes)

# Fixed AnnotSV column names up to and including FORMAT (positions 1-14).
# Everything between FORMAT and Annotation_mode is a VCF sample column.
annotsv_fixed_cols <- c(
  "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length",
  "SV_type", "Samples_ID", "ID", "REF", "ALT", "QUAL",
  "FILTER", "INFO", "FORMAT"
)

# Read and bind all SV CSVs
sv_list <- lapply(csv_files, function(f) {
  sample_name <- basename(dirname(f))
  run_name    <- basename(dirname(dirname(f)))

  df <- tryCatch(
    read_csv(f, col_types = cols(.default = col_character())),
    error = function(e) {
      message("Warning: could not read ", f, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Rename VCF sample columns to generic FORMAT_vcf[_N]
  annot_mode_pos <- match("Annotation_mode", names(df))
  if (!is.na(annot_mode_pos)) {
    vcf_sample_cols <- setdiff(
      names(df)[seq_len(annot_mode_pos - 1L)],
      annotsv_fixed_cols
    )
    if (length(vcf_sample_cols) > 0) {
      for (i in seq_along(vcf_sample_cols)) {
        new_name <- if (i == 1L) "FORMAT_vcf" else paste0("FORMAT_vcf_", i)
        names(df)[names(df) == vcf_sample_cols[i]] <- new_name
      }
    }
  }

  df$sID <- sample_name
  df$ID  <- run_name
  df
})

combined <- bind_rows(sv_list)
combined <- type_convert(combined, col_types = cols())

# Attach samplesheet metadata
combined <- left_join(combined, meta, by = c("sID", "ID"))

# Move metadata columns to the front
meta_cols <- intersect(names(meta), names(combined))
other_cols <- setdiff(names(combined), meta_cols)
combined <- combined[, c(meta_cols, other_cols)]

saveRDS(combined, output_rds)
message(sprintf("[OK] Combined %d SV rows from %d samples -> %s",
                nrow(combined), length(unique(combined$sID)), output_rds))
