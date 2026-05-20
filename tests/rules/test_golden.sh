#!/usr/bin/env bash
# Phase 3e integration test for workflow/rules/golden.smk.
# Runs the full upstream chain (repeats + pasa) then golden_genes ->
# golden_hints.
#
# CURRENT STATUS: PASSING. JID 9994 (2026-05-20, --norefine mode) ran
# the full upstream chain through golden_hints in 5:08 and produced
# 712 features. Phase 4 (`bb3e82a5`) replaced the v1 bridge with
# `bin/prepare_golden_genes.pl` + PerlLib/Golden/*.pm (pure-Perl
# validator, no AAT, no fathom). The Phase 3e+3j integration fixes
# in golden.smk pass --norefine to exonerate; --refine remains
# multi-minute per query on the mini-fixture and is not exercised
# here. The snakemake-wrapped chain downstream of golden (augustus
# + EVM + ogs) is verified separately by tests/rules/test_evm.sh.
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
