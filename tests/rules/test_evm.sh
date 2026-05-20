#!/usr/bin/env bash
# Phase 3i integration test: workflow/rules/evm.smk.
# Runs the full upstream chain (repeats + rnaseq + tgg + pasa + proteins +
# genemark + golden + augustus) then evm.smk end-to-end through evm_run.
# Verifies EVM.gff3 is produced and contains at least one EVidenceModeler
# feature row.

set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/evm"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing"; exit 1; }
done

rm -rf "$OUT"
cd "$REPO"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until evm_run

test -s "$OUT/EVM.gff3" || { echo "EVM.gff3 missing or empty"; exit 1; }

# At least one EVidenceModeler-emitted feature row (column 2 source 'EVM').
n=$(awk -F'\t' '$2 == "EVM" && $3 == "gene"' "$OUT/EVM.gff3" | wc -l)
if (( n < 1 )); then
    echo "FATAL: EVM.gff3 contains zero gene rows with source=EVM"
    head -20 "$OUT/EVM.gff3" >&2
    exit 1
fi
echo "OK: EVM.gff3 produced $n EVM gene rows"
