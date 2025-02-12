suppressMessages(library(tidyverse))
suppressMessages(library(XenofilteR))
suppressMessages(library(BiocParallel))
bpp <- SerialParam()

# Get absolute paths
wd <- getwd()
# Expected outputs
sample_name <- snakemake@wildcards$sample
snakemake_output_bam <- file.path(wd, snakemake@output$bam[[1]])
snakemake_output_bai <- file.path(wd, snakemake@output$bai[[1]])
snakemake_output_log <- file.path(wd, snakemake@output$log[[1]])
run_dir <- dirname(snakemake_output_bam)
# Xenofilter default output paths
filtered_dir <- file.path(run_dir, "Filtered_bams")
Xenofilter_output_bam <- paste0(sample_name, "_Filtered.bam")
Xenofilter_output_bai <- paste0(sample_name, "_Filtered.bam.bai")

# Pair graft and host
sample_data <- data.frame(
  graft = file.path(wd, snakemake@input[["graft_bam"]]),
  host = file.path(wd, snakemake@input[["host_bam"]])
)
# Run XenofilteR
XenofilteR(
  sample.list = sample_data,
  destination.folder = run_dir,
  bp.param = bpp,
  output.names = sample_name
)

# Move the outputs to conform to Snakemake schema
file.rename(
  from = file.path(filtered_dir, Xenofilter_output_bam),
  to = snakemake_output_bam
)
file.rename(
  from = file.path(filtered_dir, paste0(Xenofilter_output_bam, ".bai")),
  to = snakemake_output_bai
)
file.rename(
  from = file.path(filtered_dir, "XenofilteR.log"),
  to = snakemake_output_log
)
if (length(list.files(filtered_dir)) == 0) {
  unlink(filtered_dir, recursive = TRUE)
}