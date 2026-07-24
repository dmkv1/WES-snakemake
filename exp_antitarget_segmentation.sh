#!/usr/bin/env bash
# Does keeping antitargets improve CNVkit segmentation for these WES+UTR kits?
# Re-segments existing .cnr files WITH vs WITHOUT antitarget bins (BAF-informed,
# as in production) and compares segment count + on-target residual noise.
# Uses existing artifacts only (no re-alignment). Run on host.
set -euo pipefail
cd "$(dirname "$0")"

OUT="tmp/exp_antitarget"; mkdir -p "$OUT"

# cnvkit image (auto-detect among cached SIFs of this repo or the PON repo)
SIF=""
for s in .snakemake/singularity/*.simg ../WES-PON-smk/.snakemake/singularity/*.simg; do
    [ -f "$s" ] || continue
    if singularity exec "$s" cnvkit.py version >/dev/null 2>&1; then SIF="$s"; break; fi
done
[ -n "$SIF" ] || { echo "FAIL: no cnvkit image found"; exit 1; }
run() { singularity exec -B /mnt/data/NGS "$SIF" "$@"; }

seg_and_report() {
    local tag="$1" cnr="$2" vcf="$3" id="$4"
    run cnvkit.py segment "$cnr" -m cbs -v "$vcf" -i "$id" -o "$OUT/$tag.cns" >/dev/null 2>&1
    local nseg; nseg=$(tail -n +2 "$OUT/$tag.cns" | wc -l)
    # on-target residual MAD: |bin.log2 - seg.log2| over target bins only
    local mad
    mad=$(run cnvkit.py segmetrics "$cnr" -s "$OUT/$tag.cns" --mad -o /dev/stdout 2>/dev/null \
          | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{w+=$(h["weight"])*$(h["mad"]);sw+=$(h["weight"])}END{if(sw>0)printf "%.4f",w/sw; else print "NA"}')
    printf "  %-16s segments=%-5d weighted_seg_MAD=%s\n" "$tag" "$nseg" "$mad"
}

for s in "P004:P004_PT_DG" "P020:P020_PT_DLBCL"; do
    run="${s%%:*}"; id="${s##*:}"
    cnr="work/cnvkit/$run/$id/$id.filtered.cnr"
    vcf="work/cnvkit/$run/$id/$id.hetsnp.vcf"
    [ -f "$cnr" ] || { echo "skip $id (no cnr)"; continue; }
    # target-only variant: drop Antitarget rows, keep header
    awk -F'\t' 'NR==1 || $4!="Antitarget"' "$cnr" > "$OUT/${id}.targetonly.cnr"
    echo "== $id (antitarget mean depth: $(awk -F'\t' 'NR>1 && $4=="Antitarget"{d+=$6;n++}END{printf "%.2f",d/n}' "$cnr")x) =="
    seg_and_report "${id}.with"    "$cnr"                         "$vcf" "$id"
    seg_and_report "${id}.without" "$OUT/${id}.targetonly.cnr"    "$vcf" "$id"
done
echo
echo "Interpretation: fewer segments + similar/lower MAD without antitargets"
echo "=> antitargets add breakpoints/noise without improving fit -> drop them."
