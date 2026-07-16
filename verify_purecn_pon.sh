#!/usr/bin/env bash
# Sanity-checks the PureCN NormalDB test run at /mnt/scratch/DM/WES-PON-smk
# before promoting it from scratch to the production PON path. Run on host
# (needs bcftools + the purecn conda env for the R sanity check).
set -uo pipefail

PON=/mnt/scratch/DM/WES-PON-smk/results/PON/purecn

for kit in SureSelectV6UTR SureSelectV8UTR; do
  echo "=== $kit ==="

  echo "-- coverage_files.list samples --"
  cov_samples=$(sed 's/.*\///; s/\.txt$//' "$PON/$kit/coverage_files.list" | sort)
  echo "$cov_samples"
  n_cov=$(echo "$cov_samples" | wc -l)

  echo "-- joint VCF samples --"
  vcf_samples=$(bcftools query -l "$PON/$kit/normals_${kit}.joint.vcf.gz" | sort)
  echo "$vcf_samples"
  n_vcf=$(echo "$vcf_samples" | wc -l)

  echo "-- diff (should be empty) --"
  diff <(echo "$cov_samples") <(echo "$vcf_samples")

  echo "coverage list: $n_cov samples | joint VCF: $n_vcf samples"

  echo "-- mapping_bias_hq_sites line count vs low_coverage_targets line count --"
  wc -l "$PON/$kit/mapping_bias_hq_sites_${kit}_hg38.bed" "$PON/$kit/low_coverage_targets_${kit}_hg38.bed" 2>/dev/null
  echo
done

echo "=== cross-kit overlap check (should normally be empty unless a normal was run on both kits) ==="
comm -12 \
  <(sed 's/.*\///; s/\.txt$//' "$PON/SureSelectV6UTR/coverage_files.list" | sort) \
  <(sed 's/.*\///; s/\.txt$//' "$PON/SureSelectV8UTR/coverage_files.list" | sort)
echo

echo "=== RDS sanity (run in the purecn conda env) ==="
for kit in SureSelectV6UTR SureSelectV8UTR; do
  Rscript -e "
    x <- readRDS('$PON/$kit/normalDB_${kit}_hg38.rds')
    cat('$kit normalDB:', class(x), '\n')
    m <- readRDS('$PON/$kit/mapping_bias_${kit}_hg38.rds')
    cat('$kit mapping_bias:', class(m), '\n')
  "
done
