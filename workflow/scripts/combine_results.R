#!/usr/bin/env Rscript

# Parse SNV vcf from Mutect2+Funcotator (GATK 4.5.0.0-4.6.1.0)
suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(VariantAnnotation)
  library(GenomicRanges)
})

extractFUNCOTATION <-
  function(vcf,
           fields = c()) {
    annotation_colnames <-
      info(vcf@metadata$header)["FUNCOTATION", "Description"] %>%
      stringr::str_remove("Functional annotation from the Funcotator tool.  Funcotation fields are: ") %>%
      stringr::str_split(pattern = "\\|") %>%
      unlist()
    
    if (length(fields) == 0) {
      fields <- annotation_colnames
    }
    
    as.data.frame(info(vcf)[["FUNCOTATION"]]) %>%
      dplyr::select("FUNCOTATION" = value) %>%
      dplyr::mutate(
        FUNCOTATION = word(FUNCOTATION, 1, sep = "\\]"),
        FUNCOTATION = str_remove_all(FUNCOTATION, "[\\[\\]]")
      ) %>%
      tidyr::separate(FUNCOTATION,
                      into = annotation_colnames,
                      sep = "\\|",
                      fill = "right") %>%
      dplyr::select(all_of(fields))
  }

decode_url <- function(x) {
  if(is.character(x)) {
    x <- gsub("_%20_", " ", x)
    x <- gsub("_%2C_", ",", x)
    x <- gsub("_%7C_", "|", x)
    return(x)
  } else {
    return(x)
  }
}

# Parse input files paths
input_tsv_sv <- snakemake@input[["sv_tsv"]]
input_file_vcf <- snakemake@input[["snv_vcf"]]
input_cns_cnv <- snakemake@input[["cnv_cns"]]
purity <- snakemake@params[["purity"]]
sample_sex <- snakemake@params[["sample_sex"]]

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
  
filter.df <- rowRanges(vcf) %>% 
  as.data.frame() %>% 
  dplyr::select("FILTER")
  
af.df <- as.data.frame(geno(vcf)[["AF"]]) %>%
  mutate_all( ~ unlist(.)) %>%
  setNames(., paste0("AF_", colnames(.)))

gt.df <- as.data.frame(geno(vcf)[["GT"]]) %>%
  mutate_all( ~ unlist(.)) %>%
  setNames(., paste0("GT_", colnames(.)))

dp.df <- as.data.frame(geno(vcf)[["DP"]]) %>%
  mutate_all( ~ unlist(.)) %>%
  setNames(., paste0("DP_", colnames(.)))

ad.df <- as.data.frame(geno(vcf)[["AD"]])

# Extract REF counts (first element)
ad_ref.df <- ad.df %>%
  mutate(across(everything(), ~ sapply(., function(x) x[1]))) %>%
  setNames(paste0("AD_REF_", colnames(.)))

# Extract ALT counts (second element)
ad_alt.df <- ad.df %>%
  mutate(across(everything(), ~ sapply(., function(x) x[2]))) %>%
  setNames(paste0("AD_ALT_", colnames(.)))

funcotation <- extractFUNCOTATION(vcf) %>% 
  dplyr::select(
    contains("Gencode_"),
    contains("HGNC_")
  )

result <- cbind(filter.df, gt.df, ad_ref.df, ad_alt.df, af.df, dp.df, funcotation) %>% 
  rownames_to_column("Variant") %>%
  mutate(Position = word(Variant, 1, sep = "_"),
         .before = "Variant") %>%
  mutate(Variant = word(Variant, 2, sep = "_")) %>%
  dplyr::select(
    Position,
    Variant,
    FILTER,
    contains("AF_"),
    contains("GT_"),
    contains("DP_"),
    contains("AD_"),
    contains("Gencode_"),
    "HGNC_HGNC_ID",
    "HGNC_Approved_name",
    "HGNC_Locus_type",
    "HGNC_Locus_group",
    "HGNC_Alias_symbols",
    "HGNC_Alias_names",
    "HGNC_Chromosome",
    "HGNC_Accession_numbers",
    "HGNC_Enzyme_IDs",
    "HGNC_NCBI_Gene_ID",
    "HGNC_Ensembl_gene_ID",
    "HGNC_Pubmed_IDs",
    "HGNC_RefSeq_IDs",
    "HGNC_Gene_group_ID",
    "HGNC_Gene_group_name",
    "HGNC_UniProt_ID(supplied_by_UniProt)",
    "HGNC_Ensembl_ID(supplied_by_Ensembl)"
  )

result_snv <- as.data.frame(lapply(result, function(col) {
  if(is.character(col)) {
    return(sapply(col, decode_url))
  } else {
    return(col)
  }
}))

# --- Parse CNVkit cns ----
CNVs.df <- read_tsv(input_cns_cnv, show_col_types = FALSE)

cnv_gr <- GRanges(
  seqnames = CNVs.df$chromosome,
  ranges = IRanges(start = CNVs.df$start, end = CNVs.df$end),
  cn = CNVs.df$cn
)

var_chr <- word(result_snv$Position, 1, sep = ":")
var_pos <- as.integer(word(result_snv$Position, 2, sep = ":"))

var_gr <- GRanges(
  seqnames = var_chr,
  ranges = IRanges(start = var_pos, end = var_pos)
)

# Find overlapping CNV segments for each variant
overlaps <- findOverlaps(var_gr, cnv_gr, select = "first")
local_cn <- cnv_gr$cn[overlaps]
# Handle missing overlaps - default based on chromosome and sex
local_cn <- ifelse(
  is.na(local_cn),
  ifelse(
    var_chr %in% c("chrX", "X") & sample_sex == "male",
    1,  # X in males
    ifelse(
      var_chr %in% c("chrY", "Y"),
      ifelse(sample_sex == "male", 1, 0),  # Y in males=1, females=0
      2  # Autosomes default to diploid
    )
  ),
  local_cn
)

# Get tumor AF (last column in AF columns, which is the tumor sample)
af_cols <- grep("^AF_", colnames(result_snv))
tumor_af <- result_snv[, af_cols[length(af_cols)]]  # Last AF column is tumor

# Get tumor GT (last column in GT columns)
gt_cols <- grep("^GT_", colnames(result_snv))
tumor_gt <- result_snv[, gt_cols[length(gt_cols)]]

# Normalize GT format - replace phasing delimiter with /
tumor_gt_normalized <- gsub("\\|", "/", tumor_gt)

expected_mutant_copies <- ifelse(
  tumor_gt_normalized == "1/1",
  local_cn,  # Homozygous - all copies mutant
  1  # Heterozygous (0/1, 1/0) - 1 copy mutant
)

# Calculate CCF
# Formula: CCF = (observed_AF × local_CN) / (purity × expected_mutant_copies)
# 
# Logic: 
# - observed_AF = (mutant_reads) / (total_reads)
# - In pure tumor: AF = expected_mutant_copies / local_CN
# - With purity p: AF = (p × expected_mutant_copies/local_CN + (1-p) × 0)
#                     = p × expected_mutant_copies / local_CN
# - But if mutation is subclonal (present in fraction f of cancer cells):
#   AF = p × f × expected_mutant_copies / local_CN
# - Solving for f (which is CCF):
#   CCF = AF × local_CN / (p × expected_mutant_copies)
# Calculate CCF
# Formula: CCF = (observed_AF × purity × local_CN) / expected_mutant_copies

ccf <- (tumor_af * local_cn) / (purity * expected_mutant_copies)
ccf <- pmin(ccf, 1.0)  # Cap at 1.0 (100%)

# Add CCF columns to result
result_snv$local_CN <- local_cn
result_snv$expected_mutant_copies <- expected_mutant_copies
result_snv$CCF <- round(ccf, 3)

# Reorder columns
result_snv <- result_snv %>%
  dplyr::select(
    Position,
    Variant,
    FILTER,
    contains("AF_"),
    contains("GT_"),
    local_CN,
    expected_mutant_copies,
    CCF,
    everything()
  )

# Parse the .tsv with the SVs from Manta
SVs.df <- read_tsv(input_tsv_sv, show_col_types = F)

# Combine and write the output
result_list <- list(
  "SNVs" = result_snv,
  "CNVs" = CNVs.df,
  "SVs" = SVs.df
)

output_file_path <- snakemake@output[["xlsx"]]
write.xlsx(result_list, output_file_path, overwrite = TRUE)
