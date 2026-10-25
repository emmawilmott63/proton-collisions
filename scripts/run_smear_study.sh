#!/bin/bash

# Smearing study: for each relative pT smearing level, rebuild the candidate
# tree from the SAME particle file, train a fresh BDT, apply it, and collect
# one summary line per level in data/tmva_summary.csv.
#
# Usage (from the project root):
#   ./scripts/run_smear_study.sh [particleFile] [level level ...]
# Example:
#   ./scripts/run_smear_study.sh data/d0_inc_5_100.root 0 0.01 0.03 0.05 0.10
#
# Extra options for every step come from the environment, e.g.
#   BUILD_OPTS="--smear-eta=0.002 --smear-phi=0.002" \
#   TRAIN_OPTS="--depth=4"  APPLY_OPTS="--lumi-pb=10" \
#   ./scripts/run_smear_study.sh

INPUT=${1:-data/d0_inc_5_100.root}
[ $# -gt 0 ] && shift
LEVELS=${@:-0 0.01 0.03 0.05 0.10}

for EXE in ./buildtree ./trainTMVA ./applyTMVA; do
    if [ ! -x "$EXE" ]; then
        echo "Error: $EXE not found. Run 'make' in the project root first."
        exit 1
    fi
done
if [ ! -f "$INPUT" ]; then
    echo "Error: $INPUT not found."
    exit 1
fi

mkdir -p data logs plots

for S in $LEVELS; do
    TAG="smear$(echo "$S" | tr '.' 'p')"
    CAND="data/candidates_${TAG}.root"
    LOG="logs/${TAG}"

    echo "=== relative pT smearing ${S} (tag ${TAG}) ==="

    ./buildtree "$INPUT" --smear-pt="$S" --tag="$TAG" $BUILD_OPTS \
        > "${LOG}_build.log" 2>&1 \
        || { echo "buildtree failed, see ${LOG}_build.log"; exit 1; }

    ./trainTMVA "$CAND" --dataset="ds_${TAG}" $TRAIN_OPTS \
        > "${LOG}_train.log" 2>&1 \
        || { echo "trainTMVA failed, see ${LOG}_train.log"; exit 1; }

    ./applyTMVA "$CAND" --dataset="ds_${TAG}" $APPLY_OPTS \
        > "${LOG}_apply.log" 2>&1 \
        || { echo "applyTMVA failed, see ${LOG}_apply.log"; exit 1; }

    grep -E "Signal mass|ROC integral|Best cut" "${LOG}_apply.log"
done

echo
echo "Summary (data/tmva_summary.csv):"
column -s, -t data/tmva_summary.csv
