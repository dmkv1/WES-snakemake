#!/usr/bin/env Rscript

# Parse SNV vcf from Mutect2+VEP
suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(VariantAnnotation)
  library(GenomicRanges)
})

# VEP's --pick keeps exactly one (MANE Select-preferred) transcript
# annotation per allele in CSQ, so no per-transcript fan-out here.
extractCSQ <-
  function(vcf,
           fields = c()) {
    annotation_colnames <-
      info(vcf@metadata$header)["CSQ", "Description"] %>%
      stringr::str_split(pattern = "Format: ") %>%
      purrr::pluck(1, 2) %>%
      stringr::str_split(pattern = "\\|") %>%
      unlist()

    if (length(fields) == 0) {
      fields <- annotation_colnames
    }

    csq_first <- vapply(
      info(vcf)[["CSQ"]],
      function(x) if (length(x) == 0) NA_character_ else x[1],
      character(1)
    )

    tibble(CSQ = csq_first) %>%
      tidyr::separate(CSQ,
                      into = annotation_colnames,
                      sep = "\\|",
                      fill = "right") %>%
      dplyr::select(all_of(fields))
  }

# VEP percent-encodes reserved delimiter characters (e.g. "%3D" for "=")
# inside free-text CSQ subfields.
decode_url <- function(x) {
  if (is.character(x)) {
    return(vapply(x, function(v) if (is.na(v)) v else utils::URLdecode(v), character(1)))
  } else {
    return(x)
  }
}

# Parse input files paths
input_tsv_sv <- snakemake@input[["sv_tsv"]]
input_file_vcf <- snakemake@input[["snv_vcf"]]
input_cns_cnv <- snakemake@input[["cnv_cns"]]
sample_sex <- snakemake@params[["sample_sex"]]
normal_name <- snakemake@params[["normal"]]

# Single source of truth for purity, resolved by resolve_purity_source (see
# purecn.smk) — the same value cnvkit_call used for --purity, so CCF math
# here and the copy number baked into input_cns_cnv can never disagree.
purity_info <- read_csv(snakemake@input[["purity_csv"]], show_col_types = FALSE)
purity <- purity_info$purity[1]

# Fixed paths for debugging
# input_file_vcf <- "/home/dmitryk/Projects/WES_analysis/WES-snakemake/results/P005/P005_PT/P005_PT.SNV.vcf"
# input_cns_cnv <- "/home/dmitryk/Projects/WES_analysis/WES-snakemake/work/cnvkit/P005/P005_PT/P005_PT.call.cns"
# input_tsv_sv <- "/home/dmitryk/Projects/WES_analysis/WES-snakemake/work/manta/P005/P005_PT/P005_PT.SV.annotated.tsv"
# purity = 1
# sample_sex = "male"

# --- Parse Mutect2 vcf ----
vcf <- suppressWarnings(VariantAnnotation::readVcf(input_file_vcf, "hg38"))
vcf <- vcf[rowRanges(vcf)$FILTER == "PASS"]

allowed_chrs <- c(paste0("chr", (1:22)), "chrX", "chrY", "chrM", "chrMT")
vcf <- vcf[word(names(vcf), 1, sep = ":") %in% allowed_chrs, ]

if (length(vcf) == 0) {
  # No PASS variants — create empty result with expected columns
  result_snv <- tibble(
    Position = character(), Variant = character(), FILTER = character(),
    AF_tumor = numeric(), GT_tumor = character(),
    DP_tumor = integer(), AD_REF_tumor = integer(), AD_ALT_tumor = integer(),
    normal_cn = integer(), tumor_cn = integer(),
    tumor_cn1 = integer(), tumor_cn2 = integer(),
    expected_mutant_copies = integer(), CCF = numeric()
  )
} else {
  filter.df <- rowRanges(vcf) %>%
    as.data.frame() %>%
    dplyr::select("FILTER")

  # Determine generic role names for genotype columns
  # Mutect2 orders sample columns alphabetically, NOT always [normal, tumor].
  # Use the known normal sample name to assign roles correctly.
  n_geno_samples <- ncol(geno(vcf)[["AF"]])
  if (n_geno_samples == 2) {
    vcf_sample_names <- colnames(geno(vcf)[["AF"]])
    ctrl_idx <- match(normal_name, vcf_sample_names)
    if (is.na(ctrl_idx)) {
      stop("Normal sample '", normal_name, "' not found in VCF columns: ",
           paste(vcf_sample_names, collapse = ", "))
    }
    geno_roles <- character(2)
    geno_roles[ctrl_idx] <- "ctrl"
    geno_roles[3 - ctrl_idx] <- "tumor"
  } else {
    geno_roles <- "tumor"
  }

  af.df <- as.data.frame(geno(vcf)[["AF"]]) %>%
    mutate_all( ~ unlist(.)) %>%
    setNames(paste0("AF_", geno_roles))

  gt.df <- as.data.frame(geno(vcf)[["GT"]]) %>%
    mutate_all( ~ unlist(.)) %>%
    setNames(paste0("GT_", geno_roles))

  dp.df <- as.data.frame(geno(vcf)[["DP"]]) %>%
    mutate_all( ~ unlist(.)) %>%
    setNames(paste0("DP_", geno_roles))

  ad.df <- as.data.frame(geno(vcf)[["AD"]])

  # Extract REF counts (first element)
  ad_ref.df <- ad.df %>%
    mutate(across(everything(), ~ sapply(., function(x) x[1]))) %>%
    setNames(paste0("AD_REF_", geno_roles))

  # Extract ALT counts (second element)
  ad_alt.df <- ad.df %>%
    mutate(across(everything(), ~ sapply(., function(x) x[2]))) %>%
    setNames(paste0("AD_ALT_", geno_roles))

  csq <- extractCSQ(vcf)

  result <- cbind(filter.df, gt.df, ad_ref.df, ad_alt.df, af.df, dp.df, csq) %>%
    rownames_to_column("Variant") %>%
    mutate(Position = word(Variant, 1, sep = "_"),
           .before = "Variant") %>%
    mutate(Variant = word(Variant, 2, sep = "_")) %>%
    dplyr::select(
      Position, Variant, FILTER,
      matches("^AF_"), matches("^GT_"), matches("^DP_"), matches("^AD_"),
      everything()
    )

  result_snv <- as.data.frame(lapply(result, function(col) {
    if(is.character(col)) {
      return(sapply(col, decode_url))
    } else {
      return(col)
    }
  }))

  # --- Parse CNVkit cns for SNV-CNV overlap ----
  CNVs.df <- read_tsv(input_cns_cnv, show_col_types = FALSE)

  # Check if allele-specific CN columns exist (not available in tumor-only mode)
  has_allele_cn <- all(c("cn1", "cn2") %in% colnames(CNVs.df))

  # Convert to GRanges to overlap with the variants
  cnv_gr <- GRanges(
    seqnames = CNVs.df$chromosome,
    ranges = IRanges(start = CNVs.df$start, end = CNVs.df$end),
    cn = CNVs.df$cn,
    cn1 = if (has_allele_cn) CNVs.df$cn1 else NA_integer_,
    cn2 = if (has_allele_cn) CNVs.df$cn2 else NA_integer_
  )

  var_chr <- word(result_snv$Position, 1, sep = ":")
  var_pos <- as.integer(word(result_snv$Position, 2, sep = ":"))

  var_gr <- GRanges(
    seqnames = var_chr,
    ranges = IRanges(start = var_pos, end = var_pos)
  )

  # Find overlapping CNV segments for each variant
  overlaps <- findOverlaps(var_gr, cnv_gr, select = "first")
  tumor_cn <- cnv_gr$cn[overlaps]

  # Handle missing overlaps - default based on chromosome and sex
  tumor_cn <- ifelse(
    is.na(tumor_cn),
    ifelse(
      var_chr %in% c("chrX", "X") & sample_sex == "male",
      1,  # X in males
      ifelse(
        var_chr %in% c("chrY", "Y"),
        ifelse(sample_sex == "male", 1, 0),  # Y in males=1, females=0
        2  # Autosomes default to diploid
      )
    ),
    tumor_cn
  )

  # Calculate CCF using the generic tumor column names
  tumor_af <- result_snv[["AF_tumor"]]
  tumor_gt <- result_snv[["GT_tumor"]]

  # Normalize GT format - replace phasing delimiter with /
  tumor_gt_normalized <- gsub("\\|", "/", tumor_gt)

  expected_mutant_copies <- ifelse(
    tumor_gt_normalized == "1/1",
    tumor_cn,  # Homozygous - all copies mutant
    1  # Heterozygous (0/1, 1/0) - 1 copy mutant
  )

  # Calculate CCF
  # Formula: CCF = (observed_AF × tumor_cn) / (purity × expected_mutant_copies)
  #
  # Logic:
  # - observed_AF = (mutant_reads) / (total_reads)
  # - In pure tumor: AF = expected_mutant_copies / tumor_cn
  # - With purity p: AF = (p × expected_mutant_copies/tumor_cn + (1-p) × 0)
  #                     = p × expected_mutant_copies / tumor_cn
  # - But if mutation is subclonal (present in fraction f of cancer cells):
  #   AF = p × f × expected_mutant_copies / tumor_cn
  # - Solving for f (which is CCF):
  #   CCF = AF × tumor_cn / (p × expected_mutant_copies)
  # Calculate CCF
  # Formula: CCF = (observed_AF × purity × tumor_cn) / expected_mutant_copies

  ccf <- (tumor_af * tumor_cn) / (purity * expected_mutant_copies)
  ccf <- pmin(ccf, 1.0)  # Cap at 1.0 (100%)

  # Add CCF columns to result
  result_snv$expected_mutant_copies <- expected_mutant_copies
  result_snv$CCF <- round(ccf, 3)

  # Add normal CN
  normal_cn <- ifelse(
    var_chr %in% c("chrX", "X") & sample_sex == "male", 1,  # X in males
    ifelse(
      var_chr %in% c("chrY", "Y"),
      ifelse(sample_sex == "male", 1, 0),  # Y in males=1, females=0
      2  # Autosomes default to diploid
    ))
  result_snv$normal_cn <- normal_cn

  # Add tumor CN
  result_snv$tumor_cn <- tumor_cn

  # Add major and minor allele CNs
  tumor_cn1 <- cnv_gr$cn1[overlaps]
  result_snv$tumor_cn1 <- tumor_cn1

  tumor_cn2 <- cnv_gr$cn2[overlaps]
  result_snv$tumor_cn2 <- tumor_cn2

  # Reorder columns
  result_snv <- result_snv %>%
    dplyr::select(
      # Core variant info
      Position, Variant, FILTER,
      # Sample-level genotype data
      matches("^AF_"), matches("^GT_"), matches("^DP_"), matches("^AD_"),
      # Copy number & clonality
      normal_cn, tumor_cn, tumor_cn1, tumor_cn2, expected_mutant_copies, CCF,
      # Gene / consequence (VEP CSQ, --pick'd transcript)
      any_of(c(
        "SYMBOL", "Gene", "Consequence", "IMPACT", "BIOTYPE", "Feature",
        "HGVSc", "HGVSp", "Protein_position", "Amino_acids",
        "CANONICAL", "MANE", "MANE_SELECT"
      )),
      # Population frequencies
      any_of(c("AF", "MAX_AF", "MAX_AF_POPS")),
      starts_with("gnomADe_"), starts_with("gnomADg_"),
      # Clinical & somatic databases
      any_of(c("Existing_variation", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED")),
      # Functional predictions
      any_of(c("SIFT", "PolyPhen")),
      # Catch-all for any future additions
      everything()
    )
} # end of non-empty VCF branch

# --- Parse CNVkit cns (always needed for output) ----
if (!exists("CNVs.df")) {
  CNVs.df <- read_tsv(input_cns_cnv, show_col_types = FALSE)
}

# Parse the consensus (Manta + DELLY) SV annotations from AnnotSV. Report one
# row per SV (AnnotSV "full" mode) with a somatic-relevant column subset; the
# per-gene "split" rows stay in the on-disk TSV but would flood the sheet (a
# single large SV overlaps hundreds of genes). Decode SURVIVOR's SUPP_VEC
# (bit 1 = Manta, bit 2 = DELLY) into a readable SV_callers column.
#
# SV_chrom2 (INFO CHR2) is the only place the partner chromosome survives for
# TRA/BND calls — SV_chrom/SV_start/SV_end from AnnotSV collapse a
# translocation to its first breakend, indistinguishable from an intra-
# chromosomal event without it. Equal to SV_chrom for every other SV type.
#
# SV_DR_ref/SV_DR_alt are the tumor sample's own FORMAT DR ("# supporting
# reference,variant reads in that order", SURVIVOR spec), passed through
# unmodified from whichever caller(s) support the consensus call — not a
# native AF field (SURVIVOR/AnnotSV never compute one); derive one downstream
# as SV_DR_alt / (SV_DR_ref + SV_DR_alt) where SV_DR_ref + SV_DR_alt > 0.
sv_report_cols <- c(
  "SV_callers", "AnnotSV_ID", "SV_chrom", "SV_chrom2", "SV_start", "SV_end",
  "SV_length", "SV_type", "Gene_count", "Gene_name", "FILTER",
  "AnnotSV_ranking_score", "ACMG_class", "AnnotSV_ranking_criteria",
  "SV_DR_ref", "SV_DR_alt"
)

SVs.raw <- tryCatch(
  read_tsv(input_tsv_sv, show_col_types = FALSE),
  error = function(e) tibble()
)

if (nrow(SVs.raw) > 0 && "Annotation_mode" %in% names(SVs.raw)) {
  SVs.df <- SVs.raw %>% dplyr::filter(Annotation_mode == "full")

  caller_names <- c("Manta", "DELLY")
  supp_vec <- ifelse(
    grepl("SUPP_VEC=", SVs.df$INFO),
    sub(".*SUPP_VEC=([01]+).*", "\\1", SVs.df$INFO),
    NA_character_
  )
  SVs.df$SV_callers <- vapply(supp_vec, function(v) {
    if (is.na(v)) return(NA_character_)
    bits <- strsplit(v, "")[[1]] == "1"
    paste(caller_names[seq_along(bits)][bits], collapse = ";")
  }, character(1))

  SVs.df$SV_chrom2 <- ifelse(
    grepl("CHR2=", SVs.df$INFO),
    sub(".*CHR2=([^;]+).*", "\\1", SVs.df$INFO),
    NA_character_
  )

  # TRA/BND: AnnotSV can emit two "full" rows per event, a primary record plus
  # a synthesized "BNDrescue" mate at the other breakend (tagged in INFO). The
  # rescue row copies the primary's INFO verbatim, so its own CHR2 is
  # self-referential (points at itself) rather than at the true partner.
  # Collapse each primary/rescue pair (sharing the raw VCF ID column) into the
  # primary row, whose CHR2 is correct, and drop the now-redundant rescue row.
  # A rescue row without a surviving primary mate is kept — it is the only
  # record of that event — but its SV_chrom2 is unresolvable from its own
  # INFO, so it is cleared to NA rather than left silently wrong.
  is_rescue <- grepl("BNDrescue", SVs.df$INFO)
  if (any(is_rescue)) {
    primary_chrom_by_id <- setNames(SVs.df$SV_chrom[!is_rescue], SVs.df$ID[!is_rescue])
    rescue_has_primary <- is_rescue & SVs.df$ID %in% names(primary_chrom_by_id)
    SVs.df$SV_chrom2[is_rescue & !rescue_has_primary] <- NA_character_
    SVs.df <- SVs.df[!rescue_has_primary, ]
  }

  # Tumor genotype column is named by sample ID (SURVIVOR/AnnotSV carry the
  # VCF sample columns through verbatim), matching this rule's own {sample}
  # wildcard. FORMAT key order is not guaranteed fixed, so look up DR's
  # position per row rather than assuming it.
  tumor_col <- snakemake@wildcards[["sample"]]
  if (tumor_col %in% names(SVs.df) && "FORMAT" %in% names(SVs.df)) {
    format_keys <- strsplit(SVs.df$FORMAT, ":")
    geno_vals   <- strsplit(SVs.df[[tumor_col]], ":")
    dr_idx <- vapply(format_keys, function(k) {
      i <- match("DR", k)
      if (is.na(i)) NA_integer_ else i
    }, integer(1))
    dr_str <- mapply(function(v, i) {
      if (is.na(i) || i > length(v)) NA_character_ else v[i]
    }, geno_vals, dr_idx)
    dr_split <- strsplit(dr_str, ",")
    SVs.df$SV_DR_ref <- suppressWarnings(as.integer(vapply(
      dr_split, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)
    )))
    SVs.df$SV_DR_alt <- suppressWarnings(as.integer(vapply(
      dr_split, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)
    )))
  } else {
    SVs.df$SV_DR_ref <- NA_integer_
    SVs.df$SV_DR_alt <- NA_integer_
  }

  # Fixed report schema (stable across samples for downstream row-binding).
  SVs.df <- SVs.df %>% dplyr::select(dplyr::any_of(sv_report_cols))
  for (m in setdiff(sv_report_cols, names(SVs.df))) SVs.df[[m]] <- NA_character_
  SVs.df <- SVs.df %>% dplyr::select(dplyr::all_of(sv_report_cols))
} else {
  SVs.df <- as_tibble(setNames(
    rep(list(character()), length(sv_report_cols)), sv_report_cols
  ))
}

# Per-sample QC row: which purity source was used, plus PureCN's raw
# estimate for comparison even on samples where it wasn't used (see
# resolve_purity_source in purecn.smk).
qc_row <- purity_info %>%
  mutate(sample_sex = sample_sex, normal = normal_name, .before = 1)

# Combine and write the output
result_list <- list(
  "SNVs" = result_snv,
  "CNVs" = CNVs.df,
  "SVs" = SVs.df,
  "Sample_QC" = qc_row
)

output_file_path <- snakemake@output[["xlsx"]]
write.xlsx(result_list, output_file_path, overwrite = TRUE)

# Also write per-type CSVs for downstream row-binding (avoids xlsx type guessing)
write_csv(result_snv, snakemake@output[["snv_csv"]])
write_csv(CNVs.df,    snakemake@output[["cnv_csv"]])
write_csv(SVs.df,     snakemake@output[["sv_csv"]])
write_csv(qc_row,     snakemake@output[["qc_csv"]])
