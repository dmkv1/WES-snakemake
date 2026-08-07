#!/usr/bin/env bash
set -euo pipefail

# Run settings (cores, conda, singularity, bind mounts, etc.) live in
# profiles/default/config.yaml. config.yaml is loaded via the Snakefile.
# Pass -n (or any extra args) through to snakemake, e.g. ./launch.sh -n
trap 'rm -f snakemake.pid' EXIT

{
    snakemake \
        --profile profiles/default \
        "$@" &
    echo $! > snakemake.pid
    wait $!
} 2>&1 | tee snakemake.log
