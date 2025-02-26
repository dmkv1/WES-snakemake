suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(VariantAnnotation)
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

input_file_vcf <- snakemake@input[["vcf"]]

vcf <- suppressWarnings(VariantAnnotation::readVcf(input_file_vcf, "hg38"))
vcf <- vcf[rowRanges(vcf)$FILTER == "PASS"]
vcf <- vcf[!(grepl(names(vcf), pattern = "chrUn_")),]
  
funcotation <- extractFUNCOTATION(vcf) %>% 
    dplyr::select(
      contains("Gencode_"),
      contains("HGNC_")
    )
  
vaf.df <- as.data.frame(geno(vcf)[["VAF"]]) %>%
  mutate_all( ~ unlist(.)) %>%
  setNames(., paste0("VAF_", colnames(.)))

gt.df <- as.data.frame(geno(vcf)[["GT"]]) %>%
  mutate_all( ~ unlist(.)) %>%
  setNames(., paste0("GT_", colnames(.))) 
  

result <- cbind(vaf.df, gt.df) %>%
  cbind(., funcotation) %>%
  rownames_to_column("Variant") %>%
  mutate(Position_hg38 = word(Variant, 1, sep = "_"),
         .before = "Variant") %>%
  mutate(Variant = word(Variant, 2, sep = "_")) %>%
  mutate(FILTER = NA) %>% 
  dplyr::select(
    Position_hg38,
    Variant,
    "FILTER",
    VAF_NORMAL, VAF_TUMOR,
    GT_NORMAL, GT_TUMOR,
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

result <- as.data.frame(lapply(result, function(col) {
  if(is.character(col)) {
    return(sapply(col, decode_url))
  } else {
    return(col)
  }
}))

output.table <- snakemake@output[["xlsx"]]

write.xlsx(result, output.table, overwrite = TRUE)