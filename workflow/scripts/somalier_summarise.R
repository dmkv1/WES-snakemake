#!/usr/bin/env Rscript
# somalier_summarise.R -- annotate the cohort all-vs-all pairs table with
# samplesheet-derived expectations and PASS/WARN/FAIL verdicts, and render a
# full sample x sample relatedness heatmap. Called by rule somalier_summarise.
#
# Flag logic (report-only, thresholds from config params.somalier):
#   within-patient : FAIL if ibs0/n > ibs0_frac_max (true genotype conflict ->
#                    likely swap); else WARN if relatedness < relatedness_warn
#                    (LOH/aneuploidy can legitimately depress this); else PASS.
#   cross-patient  : FAIL if relatedness > unrelated_max (unexpected kinship /
#                    swap); else PASS.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

pairs_file  <- snakemake@input[["pairs"]]
samplesheet <- snakemake@input[["samplesheet"]]
ibs0_frac_max    <- as.numeric(snakemake@params[["ibs0_frac_max"]])
relatedness_warn <- as.numeric(snakemake@params[["relatedness_warn"]])
unrelated_max    <- as.numeric(snakemake@params[["unrelated_max"]])

# sample -> patient (samplesheet ID column drives the same_patient expectation)
meta <- read_csv(samplesheet, show_col_types = FALSE)
patient_of <- setNames(as.character(meta$ID), meta$sample)

pairs <- read_tsv(pairs_file, show_col_types = FALSE)
names(pairs)[1] <- sub("^#", "", names(pairs)[1])  # strip leading '#' on sample_a

pairs <- pairs %>%
  mutate(
    patient_a    = patient_of[sample_a],
    patient_b    = patient_of[sample_b],
    same_patient = patient_a == patient_b,
    expected     = if_else(same_patient, 1L, 0L),
    ibs0_frac    = if_else(n > 0, ibs0 / n, NA_real_),
    verdict = case_when(
      same_patient  & ibs0_frac > ibs0_frac_max      ~ "FAIL",
      same_patient  & relatedness < relatedness_warn ~ "WARN",
      same_patient                                   ~ "PASS",
      !same_patient & relatedness > unrelated_max    ~ "FAIL",
      TRUE                                           ~ "PASS"
    )
  )

out <- pairs %>%
  select(
    sample_a, sample_b, patient_a, patient_b, same_patient, expected,
    relatedness, ibs0, ibs2, n, ibs0_frac, concordance, verdict
  ) %>%
  arrange(factor(verdict, levels = c("FAIL", "WARN", "PASS")),
          desc(same_patient), patient_a, sample_a, sample_b)

write_tsv(out, snakemake@output[["tsv"]])

flagged <- out %>% filter(verdict != "PASS")
if (nrow(flagged) > 0) {
  message(sprintf("[somalier] %d pair(s) flagged:", nrow(flagged)))
  for (i in seq_len(nrow(flagged))) {
    r <- flagged[i, ]
    message(sprintf("  [%s] %s <-> %s  rel=%.3f ibs0_frac=%.4g same_patient=%s",
                    r$verdict, r$sample_a, r$sample_b, r$relatedness,
                    r$ibs0_frac, r$same_patient))
  }
} else {
  message("[somalier] all pairs PASS")
}

# --- Heatmap: full sample x sample matrix, ordered so patient blocks sit on
# the diagonal. Off-diagonal high relatedness = cross-patient swap signal.
samps <- union(pairs$sample_a, pairs$sample_b)
order_df <- tibble(sample = samps, patient = patient_of[samps]) %>%
  arrange(patient, sample)
lvl <- order_df$sample

# Symmetric long form + self-relatedness (1) on the diagonal.
sym <- bind_rows(
  pairs %>% transmute(s1 = sample_a, s2 = sample_b, relatedness, verdict),
  pairs %>% transmute(s1 = sample_b, s2 = sample_a, relatedness, verdict),
  tibble(s1 = lvl, s2 = lvl, relatedness = 1, verdict = "PASS")
) %>%
  mutate(s1 = factor(s1, levels = lvl), s2 = factor(s2, levels = rev(lvl)))

n_s <- length(lvl)
show_text <- n_s <= 40

p <- ggplot(sym, aes(s1, s2, fill = relatedness)) +
  geom_tile(color = "grey90") +
  geom_tile(
    data = filter(sym, verdict == "FAIL"),
    color = "red", linewidth = 0.7, fill = NA
  ) +
  scale_fill_viridis_c(limits = c(0, 1), oob = scales::squish, name = "relatedness") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )

# Numeric labels with a black halo (offset black copies behind white text) so
# they stay legible across the whole viridis range without an extra dependency.
if (show_text) {
  halo <- 0.03
  dirs <- expand.grid(dx = c(-1, 0, 1), dy = c(-1, 0, 1))
  dirs <- dirs[!(dirs$dx == 0 & dirs$dy == 0), ]
  for (k in seq_len(nrow(dirs))) {
    p <- p + geom_text(aes(label = sprintf("%.2f", relatedness)),
                       nudge_x = dirs$dx[k] * halo, nudge_y = dirs$dy[k] * halo,
                       size = 2.5, color = "black")
  }
  p <- p + geom_text(aes(label = sprintf("%.2f", relatedness)),
                     size = 2.5, color = "white")
}

dim_in <- max(4, 0.35 * n_s + 2)
ggsave(snakemake@output[["heatmap_png"]], p, width = dim_in, height = dim_in,
       dpi = 300, bg = "white", limitsize = FALSE)

message(sprintf("[OK] Summarised %d pairs across %d samples", nrow(pairs), n_s))
