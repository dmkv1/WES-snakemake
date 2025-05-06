#!/bin/bash
set -e

nohup snakemake --use-conda --cores 32 --rerun-incomplete > snakemake.log 2>&1 & disown
echo $! > snakemake.pid
