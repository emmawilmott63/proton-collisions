#!/bin/bash

# Smearing study: for each relative pT smearing level, rebuild the candidate
# tree from the SAME particle file, train a fresh BDT, apply it, and print
# data/tmva_summary.csv at the end.
#
# This is a thin wrapper around `make study`, which only redoes the steps
# whose inputs changed (an interrupted run resumes where it stopped).
#
# Usage (from the project root):
#   ./scripts/run_smear_study.sh [particleFile] [level level ...]
# Example:
#   ./scripts/run_smear_study.sh data/d0_inc_5_100.root 0 0.01 0.03 0.05 0.10
#
# Extra options for each step come from the environment, e.g.
#   BUILD_OPTS="--smear-eta=0.002 --smear-phi=0.002" \
#   TRAIN_OPTS="--depth=4" APPLY_OPTS="--lumi-pb=10" \
#   ./scripts/run_smear_study.sh

INPUT=${1:-data/d0_inc_5_100.root}
[ $# -gt 0 ] && shift
LEVELS=${*:-0 0.01 0.03 0.05 0.10}

exec make study INPUT="$INPUT" LEVELS="$LEVELS"
