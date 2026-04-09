#!/usr/bin/env Rscript
# Bundle combined SNV/CNV/SV RDS files into a single list RDS for Wesseract

snv <- readRDS(snakemake@input[["snv"]])
cnv <- readRDS(snakemake@input[["cnv"]])
sv  <- readRDS(snakemake@input[["sv"]])

saveRDS(list(SNV = snv, CNV = cnv, SV = sv), snakemake@output[["rds"]])
