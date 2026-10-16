#!/usr/bin/env bash
# Generate a large sample in parallel, then merge it into one file.
#
# Usage: ./run_generate.sh <charm|inc> <jobs> <events-per-job> [generator options]
# Example: ./run_generate.sh inc 8 250000 --pthat-min=5 --pthat-max=10
#
# Run it from the project folder (where ./jets-charm and ./jets-inc live).
# Each job uses a different seed. Pieces go to data/parts/, and the merged
# result is data/d0_charm.root or data/d0_inclusive.root (same names the
# other programs expect). Needs ROOT's `hadd` on your PATH.

set -eo pipefail

mode="$1"
jobs="$2"
events="$3"

if [ -z "$mode" ] || [ -z "$jobs" ] || [ -z "$events" ]; then
    echo "Usage: $0 <charm|inc> <jobs> <events-per-job> [generator options]"
    exit 1
fi
shift 3

case "$mode" in
    charm) prog=./jets-charm; name=charm;     base=1000 ;;
    inc)   prog=./jets-inc;   name=inclusive; base=2000 ;;
    *) echo "Mode must be 'charm' or 'inc'"; exit 1 ;;
esac

mkdir -p data/parts
rm -f data/parts/${name}_*.root data/parts/${name}_*.log

pids=()
for ((j = 1; j <= jobs; j++)); do
    seed=$((base + j))
    "$prog" "$events" "$seed" "data/parts/${name}_${j}.root" "$@" \
        > "data/parts/${name}_${j}.log" 2>&1 &
    pids+=($!)
done

failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
done

if [ "$failed" -ne 0 ]; then
    echo "At least one job failed; check data/parts/${name}_*.log"
    exit 1
fi

hadd -f "data/d0_${name}.root" data/parts/${name}_*.root

echo
echo "Merged ${jobs} x ${events} events into data/d0_${name}.root"
grep -h "Events with a D0" data/parts/${name}_*.log
