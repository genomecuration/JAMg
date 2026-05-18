#!/usr/bin/env bash
# Phase 3g integration test for workflow/rules/proteins.smk.
# blastx softmasked genome -> mini SwissProt subset, convert to GFF via
# blast2gff.py, emit augustus hints via gff2hints.pl.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/proteins"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: $SIF not found. Run 'make sif RM_LIB_HOST=...' first." >&2
    exit 1
fi

rm -rf "$OUT"

apptainer exec --bind "$REPO":"$REPO" "$SIF" \
    "$REPO/bin/jamg" run --config "$CFG" --cores 4 --until proteins_hints

test -s "$OUT/swissprot.blastx.gff3"       || { echo "blastx GFF missing/empty";   exit 1; }
test -e "$OUT/swissprot.blastx.gff3.hints" || { echo "augustus hints missing";     exit 1; }

n=$(grep -cvE '^(#|$)' "$OUT/swissprot.blastx.gff3"; true)
if (( n < 1 )); then
    echo "swissprot.blastx.gff3 has zero feature rows"
    exit 1
fi
echo "OK: proteins rule produced $n blast features"
