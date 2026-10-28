#!/bin/bash

# K-fold cross-validation of one configuration: for each fold, train on the
# other folds and score the held-out fold. Gives mean +- spread of every
# metric (needed to tell real differences from noise). One summary row per
# fold, tagged <NAME>_f<fold>, in TMVA-output/results/tmva_summary.csv.
#
# Usage (from the project root):
#   ./scripts/run_cv.sh NAME [candidatesFile]
# Example:
#   ./scripts/run_cv.sh cv_baseline data/candidates_inc_5_100.root
#   TRAIN_OPTS="--depth=5 --min-node=0.5" ./scripts/run_cv.sh cv_depth5
#
# Environment (all optional):
#   FOLDS=5               number of folds
#   HOLDOUT=1             1 = also set aside the holdout events (event%7==6),
#                         which tuning must never touch. 0 = use all events.
#   METHOD=BDT            method title (passed to train and apply)
#   TRAIN_OPTS="..."      extra options for trainTMVA, e.g. --depth=5,
#                         --bkg-subsample=0.2 (faster), --also=Fisher,MLP
#   APPLY_OPTS="..."      extra options for applyTMVA, e.g. --lumi-pb=10
#   FORCE=1               redo folds that already have results
#
# The summary figures and tables come from: python3 scripts/make_report.py

NAME=${1:?"Usage: $0 NAME [candidatesFile]"}
CAND=${2:-data/candidates_inc_5_100.root}
FOLDS=${FOLDS:-5}
HOLDOUT=${HOLDOUT:-1}
METHOD=${METHOD:-BDT}

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

mkdir -p TMVA-output/logs

for ((i = 0; i < FOLDS; i++)); do
    TAG="${NAME}_f${i}"
    STAMP="TMVA-output/plots/tmva_${METHOD}_${TAG}_significance.png"
    if [ -z "$FORCE" ] && [ -f "$STAMP" ]; then
        echo "=== $TAG: done already (FORCE=1 to redo) ==="
        continue
    fi
    echo "=== $TAG (fold $i of $FOLDS) ==="
    SPLIT="--nfolds=$FOLDS --fold=$i --holdout=$HOLDOUT"

    ./trainTMVA "$CAND" --dataset="ds_${TAG}" --tag="$TAG" --method="$METHOD" \
        $SPLIT $TRAIN_OPTS > "TMVA-output/logs/${TAG}_train.log" 2>&1 \
        || { echo "trainTMVA failed, see TMVA-output/logs/${TAG}_train.log"; exit 1; }

    # every method trained in this run (the main one plus --also methods)
    METHODS="$METHOD"
    for OPT in $TRAIN_OPTS; do
        case "$OPT" in --also=*) METHODS="$METHOD ${OPT#--also=}";; esac
    done
    METHODS=$(echo "$METHODS" | tr ',' ' ')

    for M in $METHODS; do
        ./applyTMVA "$CAND" --dataset="ds_${TAG}" --tag="$TAG" --method="$M" \
            $SPLIT $APPLY_OPTS > "TMVA-output/logs/${TAG}_${M}_apply.log" 2>&1 \
            || { echo "applyTMVA ($M) failed, see TMVA-output/logs/${TAG}_${M}_apply.log"; exit 1; }
        grep -E "ROC integral|Best cut" "TMVA-output/logs/${TAG}_${M}_apply.log" | sed "s/^/  [$M] /"
    done
done

echo
echo "Done: $NAME. Next: python3 scripts/make_report.py"
