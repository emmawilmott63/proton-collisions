#!/bin/bash

# Usage:
# ./scripts/run_generation.sh inc 10000 2
# Arguments: sample type, events per bin, parallel jobs

SAMPLE=${1:-inc}
NEVENTS=${2:-10000}
MAX_JOBS=${3:-2}

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

mkdir -p data logs

run_bin() {
    MIN=$1
    MAX=$2
    SEED=$3

    OUT="data/${PREFIX}_${MIN}_${MAX}.root"
    LOG="logs/${PREFIX}_${MIN}_${MAX}.log"

    echo "Starting pThat ${MIN}-${MAX} GeV..."

    "$EXE" "$NEVENTS" "$SEED" "$OUT" \
        "--pthat-min=$MIN" \
        "--pthat-max=$MAX" \
        > "$LOG" 2>&1

    STATUS=$?

    if [ "$STATUS" -ne 0 ] || [ ! -s "$OUT" ]; then
        echo "FAILED: pThat ${MIN}-${MAX}. Check $LOG"
        return 1
    fi

    echo "Finished pThat ${MIN}-${MAX} GeV."
}

# Run each bin in the background, limiting concurrency.
PIDS=""

for BIN in "5 10 51001" "10 20 51002" "20 40 51003" \
           "40 60 51004" "60 100 51005"; do

    set -- $BIN
    MIN=$1
    MAX=$2
    SEED=$3

    run_bin "$MIN" "$MAX" "$SEED" &
    PIDS="$PIDS $!"
    
    # Wait when we've reached the parallel-job limit.
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
done

# Wait for all remaining jobs.
FAILED=0
for PID in $PIDS; do
    if ! wait "$PID"; then
        FAILED=1
    fi
done

if [ "$FAILED" -ne 0 ]; then
    echo "At least one generation job failed. Check logs/."
    exit 1
fi

echo "Merging ROOT files..."

hadd -f "data/${PREFIX}_5_100.root" \
    "data/${PREFIX}_5_10.root" \
    "data/${PREFIX}_10_20.root" \
    "data/${PREFIX}_20_40.root" \
    "data/${PREFIX}_40_60.root" \
    "data/${PREFIX}_60_100.root"

echo "Done: data/${PREFIX}_5_100.root"
