#!/usr/bin/env bash
# Phase 3a integration test for workflow/rules/repeats.smk.
# Asserts that `bin/jamg run --until repeats_merge` materialises the merged
# GFF, the soft/hardmasked FASTAs, the BLAST-DB sentinel, and the augustus
# hint file, and that the merged GFF carries at least one feature.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/repeats"
GENOME_BASENAME=mini-genome.fasta

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: $SIF not found. Run 'make sif' first." >&2
    exit 1
fi

rm -rf "$OUT"

apptainer exec --bind "$REPO":"$REPO" "$SIF" \
    "$REPO/bin/jamg" run --config "$CFG" --cores 4 --until repeats_merge

test -s "$OUT/all_repeat_masks.gff3"                  || { echo "merged GFF missing";          exit 1; }
test -s "$OUT/${GENOME_BASENAME}.softmasked"          || { echo "softmasked FASTA missing";    exit 1; }
test -s "$OUT/${GENOME_BASENAME}.hardmasked"          || { echo "hardmasked FASTA missing";    exit 1; }
test -s "$OUT/${GENOME_BASENAME}.softmasked.fai"      || { echo "softmasked .fai missing";     exit 1; }
test -s "$OUT/${GENOME_BASENAME}.hardmasked.fai"      || { echo "hardmasked .fai missing";     exit 1; }
test -e "$OUT/${GENOME_BASENAME}.softmasked.blastdb.done" || { echo "BLAST-DB sentinel missing"; exit 1; }
test -s "$OUT/all_repeat_masks.gff3.hints"            || { echo "augustus hints missing";      exit 1; }

# At least one feature in the merged GFF (non-comment, non-empty line).
feats=$(grep -cP '^[^#]\S' "$OUT/all_repeat_masks.gff3" || true)
if (( feats < 1 )); then
    echo "merged GFF has no features ($feats)"
    exit 1
fi
echo "OK: repeats rule produced $feats features in all_repeat_masks.gff3"
