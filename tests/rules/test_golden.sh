#!/usr/bin/env bash
# Phase 3e integration test for workflow/rules/golden.smk.
# Runs the full upstream chain (repeats + pasa) then golden_genes ->
# golden_hints.
#
# CURRENT STATUS (2026-05-19): expected to FAIL end-to-end until Phase 4
# rewrite of the driver lands. The v1 bridge
# bin/prepare_golden_genes_for_predictors.pl requires AAT.pl + the
# 'filter' binary (the latter is an ELF compiled tool from
# 3rd_party/bin/filter); neither is bundled in jamg.sif. Three other
# integration fixes for the v1 bridge are already in golden.smk's
# current shell (--gmap_dir for the read-only SIF path, --augustus
# pointing at the actual install location, and all 5 PASA inputs
# satisfying check_for_options).
#
# Per systematic-debugging skill: 3 sequential fixes revealed new
# missing dependencies each time = architectural problem with the v1
# bridge. Plan §4 explicitly rewrites the driver to remove these
# v1-era deps; bundling AAT into jamg.sif as a workaround would be
# fighting that plan.
#
# Uses `pixi run` because the DAG spans both jamg.sif (golden + repeats)
# and pasa.sif (pasa_align / compare_transdecoder upstream).
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/golden"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    if [[ ! -f "$sif" ]]; then
        echo "ERROR: $sif not found. Run 'make sifs RM_LIB_HOST=...' first." >&2
        exit 1
    fi
done

rm -rf "$OUT"

cd "$REPO"
pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until golden_hints

test -s "$OUT/final_golden_genes.gff3.nr.golden.gff3"       || { echo "golden GFF missing/empty"; exit 1; }
test -e "$OUT/final_golden_genes.gff3.nr.golden.gff3.hints" || { echo "golden hints missing"; exit 1; }

n=$(grep -cvE '^(#|$)' "$OUT/final_golden_genes.gff3.nr.golden.gff3"; true)
if (( n < 1 )); then echo "golden GFF has zero feature rows"; exit 1; fi
echo "OK: golden rule produced $n features in final_golden_genes.gff3.nr.golden.gff3"
