#!/bin/bash
# Wrapper to run Rscript with the conda compiler in PATH.
# Use instead of calling Rscript directly:
#   ./run_r.sh pipeline/02_cnv/run_infercnv.R
#
# Or source this file to set up your shell:
#   source run_r.sh --setup

CONDA_BIN="/Users/jonathanzhao/miniconda3/envs/paralog-dep/bin"
export PATH="$CONDA_BIN:/usr/bin:/bin:/usr/sbin:/sbin"
export R_MAX_VSIZE=64Gb    # raise R vector heap limit (default is ~16GB)

if [[ "$1" == "--setup" ]]; then
  echo "PATH set. You can now call Rscript directly."
else
  "$CONDA_BIN/Rscript" "$@"
fi
