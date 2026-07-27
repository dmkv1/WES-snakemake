# Reference resources

Download locations and preparation commands for the external reference files the
pipeline expects. Paths match the `refs:` / `params:` entries in `config.yaml`.

## DELLY exclude regions

Telomere/centromere/unplaced-contig regions excluded during DELLY short-read SV
calling (`params.delly.exclude`, passed as `delly sr -x`). Per the DELLY author's
short-read recommendation. `chr`-prefixed, matching `Homo_sapiens_assembly38.fasta`.

- Target path: `/mnt/data/NGS/refs/delly/human.hg38.excl.tsv`

```bash
mkdir -p /mnt/data/NGS/refs/delly
wget -O /mnt/data/NGS/refs/delly/human.hg38.excl.tsv \
  https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv
```
