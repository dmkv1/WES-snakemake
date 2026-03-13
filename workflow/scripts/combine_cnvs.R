#!/usr/bin/env Rscript
# combine_cnvs.R — Row-bind all per-sample CNV CSVs into one combined RDS
# Called by snakemake rule combine_all_cnvs

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

csv_files   <- snakemake@input[["cnv_csvs"]]
samplesheet <- snakemake@params[["samplesheet"]]
output_rds  <- snakemake@output[["rds"]]

# Load samplesheet for metadata (drop FASTQ paths)
meta <- read_csv(samplesheet, show_col_types = FALSE) %>%
  select(-any_of(c("fq1", "fq2"))) %>%
  rename(sID = sample)

# Read and bind all CNV CSVs
cnv_list <- lapply(csv_files, function(f) {
  # Extract sample name from path: results/{run}/{sample}/{sample}_CNVs.csv
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

combined <- bind_rows(cnv_list)
combined <- type_convert(combined, col_types = cols())

# Attach samplesheet metadata
combined <- left_join(combined, meta, by = c("sID", "ID"))

# Move metadata columns to the front
meta_cols <- intersect(names(meta), names(combined))
other_cols <- setdiff(names(combined), meta_cols)
combined <- combined[, c(meta_cols, other_cols)]

saveRDS(combined, output_rds)
message(sprintf("[OK] Combined %d CNV segments from %d samples -> %s",
                nrow(combined), length(unique(combined$sID)), output_rds))
