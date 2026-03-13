#!/usr/bin/env Rscript
# combine_snvs.R — Row-bind all per-sample SNV CSVs into one combined RDS
# Called by snakemake rule combine_all_snvs

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

csv_files   <- snakemake@input[["snv_csvs"]]
samplesheet <- snakemake@params[["samplesheet"]]
output_rds  <- snakemake@output[["rds"]]

# Load samplesheet for metadata (drop FASTQ paths — not needed)
meta <- read_csv(samplesheet, show_col_types = FALSE) %>%
  select(-any_of(c("fq1", "fq2"))) %>%
  rename(sID = sample)

# Read and bind all SNV CSVs
snv_list <- lapply(csv_files, function(f) {
  # Extract sample name from path: results/{run}/{sample}/{sample}_SNVs.csv
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

  df$sID <- sample_name
  df$ID  <- run_name
  df
})

combined <- bind_rows(snv_list)

# Re-infer column types now that all samples are combined
# (individual CSVs were read as all-character to avoid type conflicts)
combined <- type_convert(combined, col_types = cols())

# Rename samplesheet 'probes' to avoid collision with CNVkit's 'probes' column
meta <- meta %>% rename(capture_kit = probes)

# Attach samplesheet metadata
combined <- left_join(combined, meta, by = c("sID", "ID"))

# Move metadata columns to the front
meta_cols <- intersect(names(meta), names(combined))
other_cols <- setdiff(names(combined), meta_cols)
combined <- combined[, c(meta_cols, other_cols)]

saveRDS(combined, output_rds)
message(sprintf("[OK] Combined %d SNV rows from %d samples → %s",
                nrow(combined), length(unique(combined$sID)), output_rds))
