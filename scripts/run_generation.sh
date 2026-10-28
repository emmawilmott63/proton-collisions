#!/bin/bash

# Generate pThat-binned samples in parallel and merge them with hadd.
#
# Usage (from the project root):
#   ./scripts/run_generation.sh [inc|charm] [events_per_bin] [parallel_jobs]
# Example:
#   ./scripts/run_generation.sh inc 1000000 2
#
# Optional environment variables:
#   EDGES="5 10 20 40 60 100"     pThat bin edges in GeV/c (default shown)
#   EVENTS_PER_BIN="5000000 2000000 1000000 1000000 1000000"
#                                  one entry per bin; overrides events_per_bin
#                                  (the 5-10 bin carries ~97% of the cross
#                                  section, so it limits the yield precision)
#   SEED_BASE=51000               bin number i (1, 2, ...) uses seed SEED_BASE+i
#   FORCE=1                       regenerate bins that already finished
#
# A bin that finished successfully leaves a <file>.ok marker and is skipped
# on the next run (so an interrupted job can be resumed). Use FORCE=1 or
# change EDGES / EVENTS / SEED_BASE together with FORCE=1 to redo it.

SAMPLE=${1:-inc}
NEVENTS=${2:-10000}
MAX_JOBS=${3:-2}

EDGES=${EDGES:-"5 10 20 40 60 100"}
SEED_BASE=${SEED_BASE:-51000}

if [ "$SAMPLE" = "inc" ]; then
    EXE="./jets-inc"
    PREFIX="d0_inc"
elif [ "$SAMPLE" = "charm" ]; then
    EXE="./jets-charm"
    PREFIX="d0_charm"
else
    echo "Usage: $0 [inc|charm] [events_per_bin] [parallel_jobs]"
    exit 1
fi

if [ ! -x "$EXE" ]; then
    echo "Error: $EXE not found or not executable."
    echo "Run this script from your project root directory."
    exit 1
fi

if ! command -v hadd > /dev/null 2>&1; then
    echo "Error: hadd not found (is ROOT set up in this shell?)."
    exit 1
fi

read -r -a EDGE_ARR <<< "$EDGES"
NBINS=$(( ${#EDGE_ARR[@]} - 1 ))
if [ "$NBINS" -lt 1 ]; then
    echo "Error: EDGES needs at least two values."
    exit 1
fi

if [ -n "$EVENTS_PER_BIN" ]; then
    read -r -a EV_ARR <<< "$EVENTS_PER_BIN"
    if [ "${#EV_ARR[@]}" -ne "$NBINS" ]; then
        echo "Error: EVENTS_PER_BIN has ${#EV_ARR[@]} entries but there are $NBINS bins."
        exit 1
    fi
else
    EV_ARR=()
    for ((b = 0; b < NBINS; b++)); do
        EV_ARR+=("$NEVENTS")
    done
fi

mkdir -p data logs

run_bin() {
    MIN=$1
    MAX=$2
    SEED=$3
    N=$4

    OUT="data/${PREFIX}_${MIN}_${MAX}.root"
    LOG="logs/${PREFIX}_${MIN}_${MAX}.log"

    if [ -z "$FORCE" ] && [ -s "$OUT" ] && [ -f "${OUT}.ok" ]; then
        echo "Skipping pThat ${MIN}-${MAX} GeV (already done; FORCE=1 to redo)."
        return 0
    fi

    rm -f "${OUT}.ok"
    echo "Starting pThat ${MIN}-${MAX} GeV: ${N} events, seed ${SEED}..."

    "$EXE" "$N" "$SEED" "$OUT" \
        "--pthat-min=$MIN" \
        "--pthat-max=$MAX" \
        > "$LOG" 2>&1

    STATUS=$?

    if [ "$STATUS" -ne 0 ] || [ ! -s "$OUT" ]; then
        echo "FAILED: pThat ${MIN}-${MAX}. Check $LOG"
        return 1
    fi

    touch "${OUT}.ok"
    echo "Finished pThat ${MIN}-${MAX} GeV."
}

# Run each bin in the background, limiting concurrency.
FILES=()
PIDS=()

for ((b = 0; b < NBINS; b++)); do
    MIN=${EDGE_ARR[$b]}
    MAX=${EDGE_ARR[$((b + 1))]}
    SEED=$((SEED_BASE + b + 1))

    FILES+=("data/${PREFIX}_${MIN}_${MAX}.root")

    run_bin "$MIN" "$MAX" "$SEED" "${EV_ARR[$b]}" &
    PIDS+=($!)

    # Wait when we've reached the parallel-job limit.
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
done

# Wait for all remaining jobs.
FAILED=0
for PID in "${PIDS[@]}"; do
    if ! wait "$PID"; then
        FAILED=1
    fi
done

if [ "$FAILED" -ne 0 ]; then
    echo "At least one generation job failed. Check logs/."
    exit 1
fi

MERGED="data/${PREFIX}_${EDGE_ARR[0]}_${EDGE_ARR[$NBINS]}.root"

echo "Merging ROOT files..."
hadd -f "$MERGED" "${FILES[@]}" || { echo "hadd failed."; exit 1; }

echo "Done: $MERGED"
