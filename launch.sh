#!/bin/bash
set -e

nohup snakemake --profile profiles/default > snakemake.log 2>&1 & disown
echo $! > snakemake.pid
