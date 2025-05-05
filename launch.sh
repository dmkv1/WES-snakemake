#!/bin/bash
set -e

nohup snakemake --use-conda --cores 32 > snakemake.log 2>&1 & disown
echo $! > snakemake.pid
