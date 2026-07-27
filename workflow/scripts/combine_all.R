#!/usr/bin/env Rscript
# combine_all.R -- row-bind every sample's per-type CSVs into three
# cross-cohort TSVs (SNVs, CNVs, SVs). Called by snakemake rule merge_results.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

samplesheet <- snakemake@params[["samplesheet"]]

# Load samplesheet for metadata (drop FASTQ paths — not needed). The samplesheet
# capture-kit column is already named 'capture_kit' (distinct from CNVkit's
# per-segment 'probes' probe-count column, so the join won't collide).
meta <- read_csv(samplesheet, show_col_types = FALSE) %>%
  select(-any_of(c("fq1", "fq2"))) %>%
  rename(sID = sample)

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

  # Drop any non-key columns the per-sample CSVs share with the samplesheet
  # (e.g. tumor_fraction) so the join doesn't emit duplicate .x/.y pairs;
  # the samplesheet copy is authoritative.
  dup_cols <- setdiff(intersect(names(combined), names(meta)), c("sID", "ID"))
  combined <- combined[, setdiff(names(combined), dup_cols)]
  combined <- left_join(combined, meta, by = c("sID", "ID"))

  meta_cols <- intersect(names(meta), names(combined))
  other_cols <- setdiff(names(combined), meta_cols)
  combined[, c(meta_cols, other_cols)]
}

combined_snvs <- combine_csvs(snakemake@input[["snv_csvs"]])
combined_cnvs <- combine_csvs(snakemake@input[["cnv_csvs"]])
combined_svs  <- combine_csvs(snakemake@input[["sv_csvs"]], rename_sv_sample_cols)
combined_qc   <- combine_csvs(snakemake@input[["qc_csvs"]])

# Collapse the annotated cohort pairs table to one row per sample and join the
# somalier relatedness summary onto combined_qc (keyed by sID = tumor sample):
#   rel_normal        relatedness to this sample's matched normal ('normal' col)
#   rel_min_patient   lowest relatedness among within-patient pairs (LOH-aware)
#   rel_worst_partner the partner sample giving that minimum
#   rel_max_unrelated highest relatedness to any cross-patient sample (swap flag)
#   somalier_flag     worst verdict (FAIL > WARN > PASS) over pairs touching it
rel <- read_tsv(snakemake@input[["relatedness_tsv"]], show_col_types = FALSE)

# Long form: one row per (sample, partner) so every pair is seen from both ends.
rel_long <- bind_rows(
  rel %>% transmute(sID = sample_a, partner = sample_b, same_patient,
                    relatedness, verdict),
  rel %>% transmute(sID = sample_b, partner = sample_a, same_patient,
                    relatedness, verdict)
)

flag_rank <- c(PASS = 1L, WARN = 2L, FAIL = 3L)
rel_summary <- rel_long %>%
  group_by(sID) %>%
  summarise(
    rel_min_patient   = suppressWarnings(min(relatedness[same_patient], na.rm = TRUE)),
    rel_worst_partner = partner[same_patient][which.min(relatedness[same_patient])][1],
    rel_max_unrelated = suppressWarnings(max(relatedness[!same_patient], na.rm = TRUE)),
    somalier_flag     = names(flag_rank)[max(flag_rank[verdict])],
    .groups = "drop"
  ) %>%
  mutate(
    rel_min_patient   = ifelse(is.finite(rel_min_patient), rel_min_patient, NA_real_),
    rel_max_unrelated = ifelse(is.finite(rel_max_unrelated), rel_max_unrelated, NA_real_)
  )

# rel_normal: relatedness of each tumor to the normal named in its QC row.
if ("normal" %in% names(combined_qc)) {
  rel_pair <- rel_long %>% select(sID, partner, relatedness)
  combined_qc <- combined_qc %>%
    left_join(rel_pair, by = c("sID", "normal" = "partner")) %>%
    rename(rel_normal = relatedness)
}

combined_qc <- left_join(combined_qc, rel_summary, by = "sID")

write_tsv(combined_snvs, snakemake@output[["snv_tsv"]])
write_tsv(combined_cnvs, snakemake@output[["cnv_tsv"]])
write_tsv(combined_svs,  snakemake@output[["sv_tsv"]])
write_tsv(combined_qc,   snakemake@output[["qc_tsv"]])

message(sprintf(
  "[OK] Combined %d SNV / %d CNV / %d SV / %d QC rows across %d samples",
  nrow(combined_snvs), nrow(combined_cnvs), nrow(combined_svs), nrow(combined_qc),
  length(unique(combined_snvs$sID))
))
