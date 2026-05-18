#!/usr/bin/env bash
# Phase 2 smoke test for bin/jamg + workflow/. Verifies the full upstream
# wiring works:
#   1. bin/jamg loads and parses the example config.
#   2. The config validates against workflow/config/config.schema.yaml.
#   3. collect_binds walks the config without raising.
#   4. The snakemake invocation reaches DAG construction.
#
# Phase 2 has no producing rule for OGS.gff3 (Phase 3 deliverable), so DAG
# construction will surface a MissingInputException for that file. That
# specific error proves every step above succeeded; the test passes on
# that signature. Once Phase 3 lands real rules, the assertion flips to
# "dry-run exits 0 and prints the DAG".
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: $SIF not found. Run 'make sif' first." >&2
    exit 1
fi

# Capture stdout+stderr and the snakemake exit code.
out=$(apptainer exec --bind "$REPO":"$REPO" "$SIF" \
    "$REPO/bin/jamg" run \
        --config "$REPO/workflow/config/config.example.yaml" \
        --dry-run 2>&1) && rc=0 || rc=$?

# Phase-2 expected signature: snakemake reached DAG construction and
# correctly identified OGS.gff3 as the unsatisfied root input.
if grep -qE "MissingInputException in rule all" <<< "$out" \
   && grep -qE "OGS\.gff3" <<< "$out"; then
    echo "OK: bin/jamg dry-run reached DAG construction (Phase 2 signature)"
    exit 0
fi

echo "FAIL: dry-run did not reach the expected DAG-construction stage" >&2
echo "----- captured output (rc=$rc) -----" >&2
echo "$out" >&2
exit 1
