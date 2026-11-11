#!/usr/bin/env bash
set -euo pipefail

# Group-writable outputs, so any group member can rerun jobs over them.
umask 002

# Run settings (cores, memory, conda, singularity, etc.) live in a workflow
# profile under profiles/. Without one, Snakemake uses profiles/default, the
# safe floor. Select another with --workflow-profile; a --profile would be
# overridden key by key by profiles/default.
#   ./launch.sh -n                               dry run, default profile
#   ./launch.sh --workflow-profile <name>        run with profiles/<name>/
# config.yaml is loaded via the Snakefile. All args pass through to snakemake.
trap 'rm -f snakemake.pid' EXIT

{
    snakemake "$@" &
    echo $! > snakemake.pid
    wait $!
} 2>&1 | tee snakemake.log
