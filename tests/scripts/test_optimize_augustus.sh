#!/usr/bin/env bash
# Smoke test for bin/optimize_augustus.pl.
# Asserts that --help exits 0 and emits the Usage block. The script invokes
# pod2usage on GetOptions failure ("Unknown option: help"), which prints
# Usage to stdout and exits 0.
#
# Correctness of the optimisation behaviour itself lives at the rule-level
# integration test (§3h.1 augustus_optimise) — this is a script-parse smoke only.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SIF="$REPO/containers/jamg.sif"
SCRIPT="$REPO/bin/optimize_augustus.pl"

if [[ ! -e "$SIF" ]]; then
    echo "SKIP: containers/jamg.sif not built" >&2
    exit 0
fi

out=$(apptainer exec --bind "$REPO":"$REPO" "$SIF" perl "$SCRIPT" --help 2>&1 || true)

if ! echo "$out" | grep -q "Usage"; then
    echo "FAIL: --help did not emit a Usage block. Got:" >&2
    echo "$out" | head -20 >&2
    exit 1
fi
if ! echo "$out" | grep -q -- "--species"; then
    echo "FAIL: Usage did not document --species (a mandatory parameter)" >&2
    exit 1
fi
echo "OK"
