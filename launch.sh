#!/usr/bin/env bash
set -euo pipefail

# Run settings (cores, conda, singularity, bind mounts, etc.) live in
# profiles/default/config.yaml. config.yaml is loaded via the Snakefile.
# Pass -n (or any extra args) through to snakemake, e.g. ./launch.sh -n
snakemake \
    --profile profiles/default \
    "$@" 2>&1 | tee snakemake.log
