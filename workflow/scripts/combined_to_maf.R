#!/usr/bin/env Rscript
# combined_to_maf.R -- derive a cohort MAF (GDC/maftools-compatible) from the
# already-built results/combined/combined_snvs.tsv. Isolated downstream of the
# merge step: a failure here leaves every other combined_* output intact.
#
# The SO-term -> Variant_Classification map and the effect-priority ordering are
# transcribed from mskcc/vcf2maf (vcf2maf.pl). Alleles are re-derived to MAF
# convention (leading/trailing-base trimming, "-" for simple indels, recomputed
# Start/End) since the VCF carries VCF-style anchored alleles.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

in_tsv  <- snakemake@input[["snv_tsv"]]
out_maf <- snakemake@output[["maf"]]

# --- vcf2maf effect priority (lower = more severe) -------------------------
# Used to collapse VEP's "&"-joined consequence list to a single SO term.
# Transcribed from mskcc/vcf2maf >= 1.6.21 (%effect_priority in vcf2maf.pl),
# restricted to terms VEP emits: the SnpEff-only synonyms vcf2maf also carries
# (conservative_inframe_*, 5_prime_UTR_premature_start_codon_gain_variant) are
# deliberately absent. Refresh both this table and so_to_maf() together -- any
# term missing from either is collected and warned about once per run.
effect_priority <- c(
  transcript_ablation = 1, exon_loss_variant = 1,
  splice_donor_variant = 2, splice_acceptor_variant = 2,
  stop_gained = 3, frameshift_variant = 3, stop_lost = 3,
  start_lost = 4, initiator_codon_variant = 4,
  disruptive_inframe_insertion = 5, disruptive_inframe_deletion = 5,
  inframe_insertion = 5, inframe_deletion = 5, protein_altering_variant = 5,
  missense_variant = 6, conservative_missense_variant = 6,
  rare_amino_acid_variant = 6,
  transcript_amplification = 7,
  splice_region_variant = 8, splice_donor_5th_base_variant = 8,
  splice_donor_region_variant = 8, splice_polypyrimidine_tract_variant = 8,
  stop_retained_variant = 9, start_retained_variant = 9, synonymous_variant = 9,
  incomplete_terminal_codon_variant = 10,
  coding_sequence_variant = 11, mature_miRNA_variant = 11, exon_variant = 11,
  `5_prime_UTR_variant` = 12, `3_prime_UTR_variant` = 12,
  non_coding_exon_variant = 13, non_coding_transcript_exon_variant = 13,
  non_coding_transcript_variant = 14, nc_transcript_variant = 14,
  intron_variant = 14, intragenic_variant = 14, INTRAGENIC = 14,
  NMD_transcript_variant = 15,
  upstream_gene_variant = 16, downstream_gene_variant = 16,
  TF_binding_site_variant = 17, TFBS_ablation = 17, TFBS_amplification = 17,
  regulatory_region_variant = 17, regulatory_region_ablation = 17,
  regulatory_region_amplification = 17, regulatory_region = 17,
  feature_elongation = 18, feature_truncation = 18,
  intergenic_variant = 19, intergenic_region = 19
)

# Unmapped SO terms are collected here and reported once by report_unmapped(),
# not warned about inline: pick_effect and so_to_maf run once per variant, so a
# single unmapped term in a real cohort would emit one warning per row and R
# collapses anything past 50 into "There were 50 or more warnings", losing the
# term names. The two buckets are distinct failures -- a term absent from
# effect_priority, versus a ranked term with no so_to_maf branch.
unmapped <- new.env(parent = emptyenv())
unmapped$priority <- character()
unmapped$class    <- character()

note_unmapped <- function(bucket, terms) {
  assign(bucket, union(get(bucket, envir = unmapped), terms), envir = unmapped)
  invisible(NULL)
}

report_unmapped <- function() {
  if (length(unmapped$priority))
    warning("[combined_to_maf] SO term(s) absent from effect_priority, ranked ",
            "last within their consequence list: ",
            paste(sort(unmapped$priority), collapse = ", "),
            "; refresh effect_priority and so_to_maf from vcf2maf.pl",
            call. = FALSE)
  if (length(unmapped$class))
    warning("[combined_to_maf] SO term(s) with no Variant_Classification, ",
            "written as Targeted_Region: ",
            paste(sort(unmapped$class), collapse = ", "),
            "; add them to so_to_maf", call. = FALSE)
  invisible(NULL)
}

pick_effect <- function(consequence) {
  if (is.na(consequence) || consequence == "") return("")
  terms <- str_split(consequence, "&", simplify = FALSE)[[1]]
  terms <- terms[terms != ""]
  if (length(terms) == 0) return("")
  pr <- effect_priority[terms]
  # An unmapped term still ranks last, as in vcf2maf, but must be reported here:
  # once the fill below runs it loses to any known term in the same list and
  # never reaches so_to_maf, so this is the only place it is visible.
  if (anyNA(pr)) note_unmapped("priority", terms[is.na(pr)])
  pr[is.na(pr)] <- 99L
  terms[which.min(pr)]
}

# --- SO term -> MAF Variant_Classification (vcf2maf GetVariantClassification)
so_to_maf <- function(effect, var_type, inframe) {
  if (effect %in% c("splice_acceptor_variant", "splice_donor_variant",
                    "transcript_ablation", "exon_loss_variant"))
    return("Splice_Site")
  if (effect == "stop_gained") return("Nonsense_Mutation")
  fs <- effect == "frameshift_variant" ||
    (effect == "protein_altering_variant" && !inframe)
  if (fs && var_type == "DEL") return("Frame_Shift_Del")
  if (fs && var_type == "INS") return("Frame_Shift_Ins")
  if (effect == "stop_lost") return("Nonstop_Mutation")
  if (effect %in% c("initiator_codon_variant", "start_lost"))
    return("Translation_Start_Site")
  if (effect %in% c("inframe_insertion", "disruptive_inframe_insertion") ||
      (effect == "protein_altering_variant" && inframe && var_type == "INS"))
    return("In_Frame_Ins")
  if (effect %in% c("inframe_deletion", "disruptive_inframe_deletion") ||
      (effect == "protein_altering_variant" && inframe && var_type == "DEL"))
    return("In_Frame_Del")
  if (effect %in% c("missense_variant", "coding_sequence_variant",
                    "conservative_missense_variant", "rare_amino_acid_variant"))
    return("Missense_Mutation")
  if (effect %in% c("transcript_amplification", "intron_variant",
                    "intragenic_variant", "INTRAGENIC"))
    return("Intron")
  if (effect %in% c("splice_region_variant", "splice_donor_5th_base_variant",
                    "splice_donor_region_variant",
                    "splice_polypyrimidine_tract_variant"))
    return("Splice_Region")
  if (effect %in% c("incomplete_terminal_codon_variant", "synonymous_variant",
                    "stop_retained_variant", "start_retained_variant",
                    "NMD_transcript_variant"))
    return("Silent")
  if (effect %in% c("mature_miRNA_variant", "exon_variant",
                    "non_coding_exon_variant", "non_coding_transcript_exon_variant",
                    "non_coding_transcript_variant", "nc_transcript_variant"))
    return("RNA")
  if (effect %in% c("5_prime_UTR_variant")) return("5'UTR")
  if (effect == "3_prime_UTR_variant") return("3'UTR")
  if (effect %in% c("TF_binding_site_variant", "TFBS_ablation",
                    "TFBS_amplification", "regulatory_region_variant",
                    "regulatory_region_ablation",
                    "regulatory_region_amplification", "regulatory_region",
                    "feature_elongation", "feature_truncation",
                    "intergenic_variant", "intergenic_region"))
    return("IGR")
  if (effect == "upstream_gene_variant") return("5'Flank")
  if (effect == "downstream_gene_variant") return("3'Flank")
  if (effect != "") note_unmapped("class", effect)
  "Targeted_Region"
}

# --- VCF-style allele -> MAF allele + Start/End + Variant_Type --------------
# Trims a shared suffix then a shared prefix (keeps left alignment), emits "-"
# for the empty side of a simple indel, and recomputes coordinates.
# This coordinate/allele transformation is intentional and deviates from the
# VCF-style Position/Variant fields in combined_snvs.tsv, so MAF coordinates are
# NOT a valid join key back to that table. Join on
# Tumor_Sample_Barcode + vcf_position + vcf_variant instead.
normalize_allele <- function(pos, ref, alt) {
  pos <- as.integer(pos)
  ref <- toupper(ref)
  alt <- toupper(strsplit(alt, ",", fixed = TRUE)[[1]][1])  # first ALT only

  if (nchar(ref) == 1L && nchar(alt) == 1L)
    return(list(start = pos, end = pos, ref = ref, alt = alt, type = "SNP"))

  # trim common suffix (never below 1 base on either side)
  while (nchar(ref) > 0L && nchar(alt) > 0L &&
         substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt)) &&
         (nchar(ref) > 1L || nchar(alt) > 1L)) {
    ref <- substr(ref, 1L, nchar(ref) - 1L)
    alt <- substr(alt, 1L, nchar(alt) - 1L)
  }
  # trim common prefix (advance start per trimmed base)
  while (nchar(ref) > 0L && nchar(alt) > 0L &&
         substr(ref, 1L, 1L) == substr(alt, 1L, 1L) &&
         (nchar(ref) > 1L || nchar(alt) > 1L)) {
    ref <- substr(ref, 2L, nchar(ref))
    alt <- substr(alt, 2L, nchar(alt))
    pos <- pos + 1L
  }

  rl <- nchar(ref); al <- nchar(alt)
  if (rl == 0L)                        # insertion: flanking bases pos-1, pos
    return(list(start = pos - 1L, end = pos, ref = "-", alt = alt, type = "INS"))
  if (al == 0L)                        # deletion
    return(list(start = pos, end = pos + rl - 1L, ref = ref, alt = "-", type = "DEL"))
  if (rl == al) {
    type <- switch(as.character(rl), "1" = "SNP", "2" = "DNP", "3" = "TNP", "ONP")
    return(list(start = pos, end = pos + rl - 1L, ref = ref, alt = alt, type = type))
  }
  # complex indel: classify by net length change
  if (al > rl)
    return(list(start = pos - 1L, end = pos, ref = "-", alt = alt, type = "INS"))
  list(start = pos, end = pos + rl - 1L, ref = ref, alt = alt, type = "DEL")
}

# --- HGVSp 3-letter -> HGVSp_Short (1-letter) for lollipop plots ------------
aa3to1 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q",
            Glu = "E", Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K",
            Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
            Tyr = "Y", Val = "V", Ter = "*", Sec = "U", Pyl = "O")
shorten_hgvsp <- function(x) {
  if (is.na(x) || x == "") return("")
  x <- sub("^[^:]*:", "", x)  # drop ENSP..: transcript prefix if present
  m <- gregexpr("[A-Z][a-z]{2}", x, perl = TRUE)
  regmatches(x, m) <- lapply(regmatches(x, m), function(z) {
    hit <- z %in% names(aa3to1)
    z[hit] <- aa3to1[z[hit]]
    z
  })
  x
}

maf_cols <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
  "Start_Position", "End_Position", "Strand", "Variant_Classification",
  "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
  "Tumor_Sample_Barcode", "HGVSc", "HGVSp", "HGVSp_Short",
  "Gene", "Transcript_ID", "t_depth", "t_ref_count", "t_alt_count", "tumor_vaf",
  "vcf_position", "vcf_variant"
)

snv <- suppressWarnings(read_tsv(in_tsv, show_col_types = FALSE,
                                 col_types = cols(.default = col_character())))

write_maf <- function(df) {
  cat("#version 2.4\n", file = out_maf)
  suppressWarnings(write_tsv(df, out_maf, na = "", append = TRUE, col_names = TRUE))
}

if (nrow(snv) == 0) {
  message("[combined_to_maf] no SNV rows in cohort; writing header-only MAF")
  write_maf(as_tibble(setNames(rep(list(character()), length(maf_cols)), maf_cols)))
  quit(save = "no", status = 0)
}

g <- function(name) if (name %in% names(snv)) snv[[name]] else rep(NA_character_, nrow(snv))

chrom <- word(snv$Position, 1, sep = fixed(":"))
pos   <- as.integer(word(snv$Position, 2, sep = fixed(":")))
ref_in <- word(snv$Variant, 1, sep = fixed("/"))
alt_in <- word(snv$Variant, 2, sep = fixed("/"))

norm <- Map(normalize_allele, pos, ref_in, alt_in)
inframe <- (abs(nchar(ref_in) - nchar(alt_in)) %% 3L) == 0L

effect <- vapply(g("Consequence"), pick_effect, character(1), USE.NAMES = FALSE)
var_type <- vapply(norm, `[[`, character(1), "type")
vclass <- mapply(so_to_maf, effect, var_type, inframe, USE.NAMES = FALSE)
report_unmapped()

# Hugo_Symbol falls back to the Ensembl gene ID, not to a shared "Unknown".
# maftools keys every per-gene summary on Hugo_Symbol, so one literal for all
# unnamed features collapses unrelated loci into a single row in oncoplots and
# gene summaries. The ENSG keeps them distinct; "Unknown" is left for the rows
# VEP gave no feature at all, mostly intergenic calls, which are genuinely not
# attributable to a gene. The Gene column below carries the ENSG either way.
symbol <- g("SYMBOL")
gene   <- g("Gene")
hugo <- ifelse(!is.na(symbol) & symbol != "", symbol,
               ifelse(!is.na(gene) & gene != "", gene, "Unknown"))

maf <- tibble(
  Hugo_Symbol           = hugo,
  Entrez_Gene_Id        = 0L,
  Center                = ".",
  NCBI_Build            = "GRCh38",
  Chromosome            = chrom,
  Start_Position        = vapply(norm, `[[`, integer(1), "start"),
  End_Position          = vapply(norm, `[[`, integer(1), "end"),
  Strand                = "+",
  Variant_Classification = vclass,
  Variant_Type          = var_type,
  Reference_Allele      = vapply(norm, `[[`, character(1), "ref"),
  Tumor_Seq_Allele1     = vapply(norm, `[[`, character(1), "ref"),
  Tumor_Seq_Allele2     = vapply(norm, `[[`, character(1), "alt"),
  Tumor_Sample_Barcode  = g("sID"),
  HGVSc                 = g("HGVSc"),
  HGVSp                 = g("HGVSp"),
  HGVSp_Short           = vapply(g("HGVSp"), shorten_hgvsp, character(1), USE.NAMES = FALSE),
  Gene                  = gene,  # Ensembl gene ID, also the Hugo_Symbol fallback
  Transcript_ID         = g("Feature"),
  t_depth               = g("DP_tumor"),
  t_ref_count           = g("AD_REF_tumor"),
  t_alt_count           = g("AD_ALT_tumor"),
  tumor_vaf             = g("AF_tumor"),
  vcf_position          = g("Position"),
  vcf_variant           = g("Variant")
)

write_maf(maf)
message(sprintf("[combined_to_maf] wrote %d variants across %d samples -> %s",
                nrow(maf), length(unique(maf$Tumor_Sample_Barcode)), out_maf))
