#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f snakemake.pid ]]; then
    echo "No snakemake.pid found. The pipeline does not run from launch.sh." >&2
    exit 1
fi

kill -TERM "$(cat snakemake.pid)"
