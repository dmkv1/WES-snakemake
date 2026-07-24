#!/usr/bin/env bash
# Experiment: does letting PureCN do its OWN segmentation (native CBS, noise-
# calibrated) instead of inheriting CNVkit's 615-segment .seg fix the
# over-segmentation / bad-purity problem?
#
# Identical to the production purecn_run EXCEPT: no --seg-file, and
# --fun-segmentation CBS (production uses --seg-file + Hclust). Everything else
# (tumor .cnr, BAF vcf, normalDB, mapping-bias, blacklist, seed) is unchanged.
#
# Run on PDX samples: these are ~pure human tumor (mouse stroma removed by
# xengsort), so ground-truth purity should be ~0.85-1.0. Production gave 0.32 /
# 0.18 -> if native segmentation recovers high purity + <300 segments, the
# CNVkit->PureCN seg handoff is the culprit.
set -euo pipefail
cd "$(dirname "$0")"

ENV=.snakemake/conda/f3f12da6c28cf2d85124c11fb35a0cef_
RS="$ENV/bin/Rscript"
PURECN=$("$RS" -e 'cat(system.file("extdata","PureCN.R",package="PureCN"))')
PON=/mnt/data/NGS/WES/PON/WES-PON-smk/results/PON/purecn/SureSelectV6UTR
BLACKLIST=/mnt/data/NGS/refs/UCSC/hg38_simpleRepeats.bed
SEG_FUN="${1:-CBS}"          # override: ./exp_...sh PSCBS  (if installed)
OUTROOT=work/purecn_native

# sample:run:sex   (both PDX, V6+UTR)
for spec in P004_PDX_VFND1:P004:F P013_PDX_VFND4:P013:M; do
  IFS=: read -r sample run sex <<< "$spec"
  cnr="work/cnvkit/$run/$sample/$sample.filtered.cnr"
  vcf="work/cnvkit/$run/$sample/$sample.hetsnp.vcf"
  sfx=$([ "$sex" = "F" ] && echo f || echo m)
  ndb="$PON/normalDB_SureSelectV6UTR_${sfx}_hg38.rds"
  mbias="$PON/mapping_bias_SureSelectV6UTR_${sfx}_hg38.rds"
  out="$OUTROOT/$sample"; mkdir -p "$out"

  # PureCN native segmentation needs coverage in ITS format (GATK3-style
  # Target/total_coverage/on_target), NOT the CNVkit .cnr (whose `depth` column
  # PureCN misreads as read counts -> everything fails the <100-read filter).
  # Reproduce the PON's purecn_coverage transform (total_coverage = depth*width)
  # so the tumor matches the NormalDB's read-count scale.
  cov="$out/$sample.purecn_cov.txt"
  awk -F'\t' 'NR==1{print "Target\ttotal_coverage\ton_target"; next}
    {ot=($4=="Antitarget")?"False":"True";
     printf "%s:%d-%d\t%s\t%s\n", $1, $2+1, $3, $6*($3-$2), ot}' "$cnr" > "$cov"

  echo "=== $sample (sex $sex) native seg=$SEG_FUN -> $out ==="
  "$RS" "$PURECN" \
      --out "$out" \
      --sampleid "$sample" \
      --tumor "$cov" \
      --sex "$sex" \
      --vcf "$vcf" \
      --genome hg38 \
      --normaldb "$ndb" \
      --mapping-bias-file "$mbias" \
      --snp-blacklist "$BLACKLIST" \
      --fun-segmentation "$SEG_FUN" \
      --min-base-quality 20 \
      --post-optimize \
      --cores 4 \
      --force --seed 42 \
      > "$out/$sample.native.log" 2>&1 || echo "  (PureCN exited nonzero - see log)"
done

echo
echo "===================== COMPARISON ====================="
echo "(segments = 'Found N segments' PureCN actually fit on; PDX ground-truth purity ~0.85-1.0)"
printf "%-18s %-12s %-9s %-8s %-7s %s\n" sample variant segments purity ploidy flagged
# grep the count PureCN reports after its own short-seg removal / merge
seg_count() { grep -oE "Found [0-9]+ segments" "$1" 2>/dev/null | tail -1 | grep -oE "[0-9]+" || echo NA; }
for spec in P004_PDX_VFND1:P004 P013_PDX_VFND4:P013; do
  IFS=: read -r sample run <<< "$spec"
  # production (Hclust on CNVkit seg)
  pc="work/purecn/$run/$sample/$sample.csv"
  pl="work/logs/purecn_${run}_${sample}.log"
  [ -f "$pc" ] && awk -F',' -v S="$sample" -v N="$(seg_count "$pl")" \
     'NR==2{gsub(/"/,"");printf "%-18s %-12s %-9s %-8s %-7.4f %s\n",S,"prod-Hclust",N,$2,$3,$6}' "$pc"
  # native
  nc="$OUTROOT/$sample/$sample.csv"
  if [ -f "$nc" ]; then
    awk -F',' -v S="$sample" -v V="native-$SEG_FUN" -v N="$(seg_count "$OUTROOT/$sample/$sample.native.log")" \
      'NR==2{gsub(/"/,"");printf "%-18s %-12s %-9s %-8s %-7.4f %s\n",S,V,N,$2,$3,$6}' "$nc"
  else
    printf "%-18s %-12s %-9s %s\n" "$sample" "native-$SEG_FUN" NA "(no csv - failed, see native.log)"
  fi
done
