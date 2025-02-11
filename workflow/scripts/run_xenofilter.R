suppressMessages(library(tidyverse))
suppressMessages(library(XenofilteR))
suppressMessages(library(BiocParallel))
bpp <- SerialParam()

wd <- getwd()

snakemake_output <- file.path(wd, snakemake@output[[1]])
run_dir <- dirname(snakemake_output)
sample_name <- snakemake@wildcards$sample

sample_data <- data.frame(
  graft = file.path(wd, snakemake@input[["human_bam"]]),
  host = file.path(wd, snakemake@input[["mouse_bam"]])
)
XenofilteR(
  sample.list = sample_data,
  destination.folder = run_dir,
  bp.param = bpp,
  output.names = sample_name
)

filtered_dir <- file.path(run_dir, "Filtered_bams")
Xenofilter_output_bam <- paste0(sample_name, "_Filtered.bam")
file.rename(
  from = file.path(filtered_dir, Xenofilter_output_bam),
  to = snakemake_output
)
file.rename(
  from = file.path(filtered_dir, paste0(Xenofilter_output_bam, ".bai")),
  to = paste0(snakemake_output, ".bai")
)
file.rename(
  from = file.path(filtered_dir, "XenofilteR.log"),
  to = file.path(run_dir, "XenofilteR.log")
)
if (length(list.files(filtered_dir)) == 0) {
  unlink(filtered_dir, recursive = TRUE)
}