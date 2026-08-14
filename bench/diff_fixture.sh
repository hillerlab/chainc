#!/usr/bin/env bash
set -euo pipefail

# Compare chainc against Kent chainCleaner on the unit fixtures, where output
# is byte-identical for a supplied NET.
#
# Usage: diff_fixture.sh <stem> <ref.2bit> <query.2bit>
#   stem ∈ {single, pair}
#
# Env:
#   CHAIN_CLEANER  Kent chainCleaner binary (default: chainCleaner)
#   CHAIN_SORT     Kent chainSort binary (default: chainSort)

stem="$1"
ref="$2"
query="$3"

CHAIN_CLEANER="${CHAIN_CLEANER:-chainCleaner}"
CHAIN_SORT="${CHAIN_SORT:-chainSort}"
export CHAIN_CLEANER CHAIN_SORT

mkdir -p "out/fixture/$stem"

"$CHAIN_CLEANER" \
    "tests/fixtures/$stem.chain" \
    "$ref" \
    "$query" \
    "out/fixture/$stem/kent.chain" \
    "out/fixture/$stem/kent.bed" \
    -net="tests/fixtures/$stem.net" \
    -linearGap=loose

target/release/chainc \
    -c "tests/fixtures/$stem.chain" \
    -r "$ref" \
    -q "$query" \
    -o "out/fixture/$stem/chainc.chain" \
    -b "out/fixture/$stem/chainc.bed" \
    -n "tests/fixtures/$stem.net" \
    -g loose

cmp "out/fixture/$stem/kent.bed" "out/fixture/$stem/chainc.bed"
cmp "out/fixture/$stem/kent.chain" "out/fixture/$stem/chainc.chain"
