#!/usr/bin/env bash
# Phase 3i integration test: workflow/rules/evm.smk.
# Pre-stages every upstream rule's output from test_suite/fixtures/snapshot/
# then runs ONLY the EVM rules (preflight + tag_* + abinitio_cat +
# stage_repeats + run). Skips the slow upstream (repeats, rnaseq, tgg,
# proteins, pasa, golden, augustus, genemark) via snakemake mtime check.
# Target runtime: under 2 minutes.
#
# Prerequisites:
#   - tools/build_fixtures.sh has been run (mini-* fixtures present)
#   - tools/regen_snapshot.sh has been run (snapshot/ tree exists)
#     OR a prior successful full-DAG run + tools/snapshot_dag_outputs.sh
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/evm"
SNAP="$REPO/test_suite/fixtures/snapshot"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing"; exit 1; }
done
[[ -d "$SNAP" ]] || {
    echo "FATAL: no snapshot at $SNAP" >&2
    echo "Run tools/regen_snapshot.sh first (regenerates from a seed full-DAG run)." >&2
    exit 1
}

cd "$REPO"
rm -rf test_suite/output .snakemake

# Pre-stage every upstream rule's outputs. EVM consumes: repeats
# (softmasked + all_repeat_masks.gff3), proteins, pasa, golden, augustus,
# rnaseq (transcript alignments via tgg), and genemark predictions.
bash tools/stage_fixtures.sh repeats rnaseq tgg proteins pasa golden augustus

# The genemark fixture is a separate committed file (subset of a 1 Mb run);
# see tests/rules/test_genemark.sh. Stage it explicitly so the snapshot
# regeneration path (which uses the committed fixture too) and the per-rule
# test path are consistent.
GM_OUT="$REPO/test_suite/output/genemark"
mkdir -p "$GM_OUT"
cp "$REPO/test_suite/fixtures/genemark.gff3" "$GM_OUT/genemark.gff3"
cp "$REPO/test_suite/fixtures/genemark.gtf"  "$GM_OUT/genemark.gtf"
: > "$GM_OUT/.preflight.ok"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" \
    --until evm_run

test -s "$OUT/EVM.gff3" || { echo "EVM.gff3 missing or empty"; exit 1; }

# At least one EVidenceModeler-emitted gene row (source 'EVM').
n=$(awk -F'\t' '$2 == "EVM" && $3 == "gene"' "$OUT/EVM.gff3" | wc -l)
if (( n < 1 )); then
    echo "FATAL: EVM.gff3 contains zero gene rows with source=EVM"
    head -20 "$OUT/EVM.gff3" >&2
    exit 1
fi
echo "OK: EVM.gff3 produced $n EVM gene rows"
