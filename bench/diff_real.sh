#!/usr/bin/env bash
set -euo pipefail

# Compare chainc against Kent chainCleaner on real data and report wall time
# plus peak RSS for both implementations.
#
# Usage: diff_real.sh <label> <chain> <ref.2bit> <query.2bit> <tSizes> <qSizes> <netfile|->
#   Last arg "-" means both implementations net the input internally.
#
# Env:
#   CHAIN_CLEANER  Kent chainCleaner binary (default: chainCleaner)
#   CHAIN_SORT     Kent chainSort binary (default: chainSort)
#   Nonet mode additionally needs chainNet and NetFilterNonNested.perl on PATH.

label="$1"
chain="$2"
ref="$3"
query="$4"
t_sizes="$5"
q_sizes="$6"
netfile="$7"

CHAIN_CLEANER="${CHAIN_CLEANER:-chainCleaner}"
CHAIN_SORT="${CHAIN_SORT:-chainSort}"
export CHAIN_CLEANER CHAIN_SORT

mkdir -p "out/$label"

if [[ "$netfile" != "-" ]]; then
    /usr/bin/time -v -o "out/$label/kent.time" \
        "$CHAIN_CLEANER" \
        "$chain" \
        "$ref" \
        "$query" \
        "out/$label/kent.chain" \
        "out/$label/kent.bed" \
        -net="$netfile" \
        -tSizes="$t_sizes" \
        -qSizes="$q_sizes" \
        -linearGap=loose \
        -newChainIDDict="out/$label/kent.ids"

    /usr/bin/time -v -o "out/$label/chainc.time" \
        target/release/chainc \
        -c "$chain" \
        -r "$ref" \
        -q "$query" \
        -o "out/$label/chainc.chain" \
        -b "out/$label/chainc.bed" \
        -R "$t_sizes" \
        -Q "$q_sizes" \
        -g loose \
        -i "out/$label/chainc.ids" \
        -n "$netfile"
else
    /usr/bin/time -v -o "out/$label/kent.time" \
        "$CHAIN_CLEANER" \
        "$chain" \
        "$ref" \
        "$query" \
        "out/$label/kent.chain" \
        "out/$label/kent.bed" \
        -tSizes="$t_sizes" \
        -qSizes="$q_sizes" \
        -linearGap=loose \
        -newChainIDDict="out/$label/kent.ids"

    /usr/bin/time -v -o "out/$label/chainc.time" \
        target/release/chainc \
        -c "$chain" \
        -r "$ref" \
        -q "$query" \
        -o "out/$label/chainc.chain" \
        -b "out/$label/chainc.bed" \
        -R "$t_sizes" \
        -Q "$q_sizes" \
        -g loose \
        -i "out/$label/chainc.ids"
fi

echo "=== $label ==="
grep -H 'Maximum resident set size (kbytes)' "out/$label/kent.time" "out/$label/chainc.time"
grep -H 'Elapsed (wall clock) time' "out/$label/kent.time" "out/$label/chainc.time"

cmp "out/$label/kent.bed" "out/$label/chainc.bed"
cmp "out/$label/kent.ids" "out/$label/chainc.ids"
cmp <(grep '^#' "out/$label/kent.chain") <(grep '^#' "out/$label/chainc.chain")
cmp <(grep -v '^#' "out/$label/kent.chain" | sort) <(grep -v '^#' "out/$label/chainc.chain" | sort)
