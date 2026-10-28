#!/bin/bash

# One-parameter scan with cross-validation: for every value, run K-fold CV
# (scripts/run_cv.sh) with --PARAM=VALUE added to trainTMVA. Tags are
# scan_<param>_<value>_f<fold>, which make_report.py turns into one figure per
# parameter with error bars.
#
# Usage (from the project root):
#   ./scripts/run_scan.sh PARAM "VALUE VALUE ..." [candidatesFile]
# Examples:
#   ./scripts/run_scan.sh depth "2 3 4 5 6"
#   ./scripts/run_scan.sh ntrees "100 200 400 800"
#   ./scripts/run_scan.sh min-node "0.25 0.5 1 2 4"
#   ./scripts/run_scan.sh boost-type "AdaBoost RealAdaBoost Grad"
#
# PARAM is any numeric or text option of trainTMVA (depth, ntrees, min-node,
# ncuts, ada-beta, shrinkage, bag-fraction, boost-type, separation, ...).
# Values must not contain '_' or spaces.
#
# Environment: the same as run_cv.sh (FOLDS, HOLDOUT, TRAIN_OPTS, APPLY_OPTS,
# FORCE). Keep the other settings fixed across a scan by putting them in
# TRAIN_OPTS. A full scan is (values x folds) trainings: start with
#   FOLDS=3 TRAIN_OPTS="--bkg-subsample=0.2"
# and repeat the interesting region at full statistics.

PARAM=${1:?"Usage: $0 PARAM \"VALUE VALUE ...\" [candidatesFile]"}
VALUES=${2:?"Usage: $0 PARAM \"VALUE VALUE ...\" [candidatesFile]"}
CAND=${3:-data/candidates_inc_5_100.root}

BASE_TRAIN_OPTS=$TRAIN_OPTS
for V in $VALUES; do
    case "$V" in *_*) echo "Error: value '$V' must not contain '_'"; exit 1;; esac
    echo "##### $PARAM = $V #####"
    TRAIN_OPTS="--${PARAM}=${V} ${BASE_TRAIN_OPTS}" \
        ./scripts/run_cv.sh "scan_${PARAM}_${V}" "$CAND" || exit 1
done

echo
echo "Scan finished. Next: python3 scripts/make_report.py"
