# Caps base qualities at MAX_BASE_QUALITY before Manta calling. BALSAMIC applies
# the same cap; uncapped binned/odd NextSeq qualities (see TODO.md #1) can skew
# Manta's scoring. Only the Manta-facing BAM is affected.
import pysam

MAX_BASE_QUALITY = 70

with pysam.AlignmentFile(snakemake.input.bam, "rb") as infile:
    with pysam.AlignmentFile(snakemake.output.bam, "wb", template=infile) as outfile:
        for read in infile:
            quals = read.query_qualities
            if quals is not None:
                read.query_qualities = [min(q, MAX_BASE_QUALITY) for q in quals]
            outfile.write(read)

pysam.index(snakemake.output.bam)
