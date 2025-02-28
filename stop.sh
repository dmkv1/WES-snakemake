#!/bin/bash
set -e

kill -TERM $(cat snakemake.pid)
