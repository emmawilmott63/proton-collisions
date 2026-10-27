#!/bin/bash

# Ablation study: retrain the BDT with one feature group removed at a time
# and compare against the full feature set. This is the most reliable
# importance measure (the BDT cannot fall back on a correlated partner),
# at the cost of one training per group.
#
# Usage (from the project root):
#   ./scripts/run_ablation.sh [candidatesFile]
# Example:
#   ./scripts/run_ablation.sh data/candidates_inc_5_100.root
#
# Groups (see tmvaCommon.h): daughters, pair, jet, isolation, event
# Extra options for every training come from the environment:
#   TRAIN_OPTS="--depth=5 --min-node=0.5"  APPLY_OPTS="--lumi-pb=10"
#   FEATURE_GROUPS="jet isolation"   only these groups
#
# Results: one row per run in data/tmva_summary.csv, tags <base>_all and
# <base>_no_<group>. Compare the auc and signif_best columns.

CAND=${1:-data/candidates_inc_5_100.root}
FEATURE_GROUPS=${FEATURE_GROUPS:-"daughters pair jet isolation event"}

for EXE in ./trainTMVA ./applyTMVA; do
    if [ ! -x "$EXE" ]; then
        echo "Error: $EXE not found. Run 'make' in the project root first."
        exit 1
    fi
done
if [ ! -f "$CAND" ]; then
    echo "Error: $CAND not found."
    exit 1
fi

BASE=$(basename "$CAND" .root)
BASE=${BASE#candidates_}

mkdir -p data logs plots

run_one() {
    TAG=$1
    shift
    echo "=== $TAG ==="
    ./trainTMVA "$CAND" --dataset="ds_${TAG}" --tag="$TAG" "$@" $TRAIN_OPTS \
        > "logs/${TAG}_train.log" 2>&1 \
        || { echo "trainTMVA failed, see logs/${TAG}_train.log"; return 1; }
    ./applyTMVA "$CAND" --dataset="ds_${TAG}" --tag="$TAG" "$@" $APPLY_OPTS \
        > "logs/${TAG}_apply.log" 2>&1 \
        || { echo "applyTMVA failed, see logs/${TAG}_apply.log"; return 1; }
    grep -E "ROC integral|Best cut" "logs/${TAG}_apply.log"
}

run_one "ablation_${BASE}_all" || exit 1
for G in $FEATURE_GROUPS; do
    run_one "ablation_${BASE}_no_${G}" --drop-group="$G" || exit 1
done

echo
echo "Ablation summary (higher auc / signif_best = better):"
( head -1 data/tmva_summary.csv; grep "^ablation_${BASE}_" data/tmva_summary.csv ) | column -s, -t
