#!/usr/bin/env bash
# Phase 3h integration test for workflow/rules/augustus.smk.
# Runs the full upstream chain (repeats + rnaseq + tgg + pasa + proteins +
# golden) then augustus_preflight -> augustus.
#
# Validates:
#   (1) augustus_results.gff3 is produced and contains at least one feature.
#   (2) The augustus_preflight sentinel exists (species + extrinsic.cfg ok).
#
# Uses `pixi run` because the DAG spans jamg.sif (most rules), pasa.sif
# (pasa_setup_db, pasa_align, pasa_compare_transdecoder), and trinity.sif
# (tgg_trinity). `pixi run` routes each rule to its declared container via
# --software-deployment-method apptainer.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/augustus"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    if [[ ! -f "$sif" ]]; then
        echo "ERROR: $sif not found. Run 'make sifs RM_LIB_HOST=...' first." >&2
        exit 1
    fi
done

rm -rf "$OUT"

cd "$REPO"
pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until augustus

test -e "$OUT/.preflight.ok"          || { echo "preflight sentinel missing"; exit 1; }
test -s "$OUT/augustus_results.gff3"  || { echo "augustus_results.gff3 missing/empty"; exit 1; }

n=$(grep -cvE '^(#|$)' "$OUT/augustus_results.gff3"; true)
if (( n < 1 )); then echo "augustus GFF has zero feature rows"; exit 1; fi
echo "OK: augustus rule produced $n features in augustus_results.gff3"
