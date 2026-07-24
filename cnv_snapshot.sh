#!/usr/bin/env bash
# Snapshot the state of the CNVkit + PureCN arm (PON artifacts and caller
# tumor outputs) to a labeled markdown file, for before/after comparison
# across the off-target (antitarget) change. Usage:
#   ./cnv_snapshot.sh <label>      e.g. pre-rebuild | post-rebuild
set -euo pipefail
cd "$(dirname "$0")"
LABEL="${1:?usage: cnv_snapshot.sh <label>}"
PON=/mnt/data/NGS/WES/PON/WES-PON-smk/results
DEV=work
DST=/mnt/data/NGS/WES/PON/cnv_snapshots
mkdir -p "$DST"
OUT="$DST/snapshot_${LABEL}.md"

{
echo "# CNV/PureCN snapshot: $LABEL"
echo "_generated $(date '+%Y-%m-%d %H:%M')_"
echo
echo "## PON artifacts (WES-PON-smk)"
echo '```'
for k in SureSelectV6UTR SureSelectV8UTR; do
  for f in targets antitargets; do
    b="$PON/PON/cnvkit/$k/$f.bed"
    printf "%-22s %-12s lines=%-8s %s\n" "$k" "$f.bed" "$(wc -l < "$b" 2>/dev/null)" "$(stat -c %y "$b" 2>/dev/null | cut -d. -f1)"
  done
  for s in m f; do
    r="$PON/PON/cnvkit/$k/reference_${s}.cnn"
    [ -f "$r" ] && awk -F'\t' -v K="$k" -v S="$s" 'NR>1{n++;if($4=="Antitarget")a++}END{printf "%-22s reference_%s   bins=%-7d antitarget=%d\n",K,S,n,a+0}' "$r"
  done
done
echo
for ok in "$PON"/purecn/*/interval_check_*.ok; do [ -f "$ok" ] && cat "$ok"; done
echo
echo "NormalDB / mapping_bias:"
find "$PON"/PON/purecn -name "*.rds" -printf '  %f  %s bytes  %TY-%Tm-%Td %TH:%TM\n' | sort
echo '```'
echo
echo "## Caller tumor outputs (WES-snakemake-dev)"
echo '```'
printf "%-26s %-8s %-10s %-7s %-7s %-9s %s\n" sample cnr_bins antitarget purity ploidy flagged comment
for d in $(find "$DEV"/purecn -maxdepth 2 -mindepth 2 -type d | sort); do
  run=$(basename "$(dirname "$d")"); s=$(basename "$d")
  csv="$d/$s.csv"; cnr="$DEV/cnvkit/$run/$s/$s.filtered.cnr"
  bins=NA; anti=NA
  [ -f "$cnr" ] && read bins anti < <(awk -F'\t' 'NR>1{n++;if($4=="Antitarget")a++}END{print n, a+0}' "$cnr")
  if [ -f "$csv" ]; then
    awk -F',' -v S="$s" -v B="$bins" -v A="$anti" 'NR==2{
      gsub(/"/,"");
      printf "%-26s %-8s %-10s %-7s %-7.4f %-9s %s\n", S, B, A, $2, $3, $6, $9}' "$csv"
  else
    printf "%-26s %-8s %-10s %s\n" "$s" "$bins" "$anti" "(no PureCN csv - failed/absent)"
  fi
done
echo
echo "cnr mtimes:"; find "$DEV"/cnvkit -name "*.filtered.cnr" -printf '  %TY-%Tm-%Td %TH:%TM  %p\n' | sort
echo '```'
} > "$OUT"
echo "wrote $OUT"
