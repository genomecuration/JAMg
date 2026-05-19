#!/usr/bin/env bash
# Phase 3c integration test for workflow/rules/tgg.smk.
# Trinity genome-guided assembly consumes the sorted RNA-seq BAM from
# rnaseq.smk and produces Trinity-GG.fasta. Trinity lives in trinity.sif;
# this test depends on both jamg.sif and trinity.sif.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
TRINITY_SIF="$REPO/containers/trinity.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/tgg"

for sif in "$JAMG_SIF" "$TRINITY_SIF"; do
    if [[ ! -f "$sif" ]]; then
        echo "ERROR: $sif not found. Run 'make sifs RM_LIB_HOST=...' first." >&2
        exit 1
    fi
done

rm -rf "$OUT"

# Run from the host (via pixi default env, which has snakemake 9.21.0) so
# snakemake's `--software-deployment-method apptainer` can dispatch per-rule
# containers. Phase 3c uses trinity.sif which is a SEPARATE SIF from jamg.sif;
# nesting apptainer-in-apptainer is not supported, so inside-SIF invocation
# does not work for this rule.
cd "$REPO"
pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until tgg_trinity

test -s "$OUT/Trinity-GG.fasta" || { echo "Trinity-GG.fasta missing or empty"; exit 1; }
# `; true` normalises grep's exit status so $() captures only the count line.
n=$(grep -c '^>' "$OUT/Trinity-GG.fasta"; true)
if (( n < 1 )); then
    echo "Trinity-GG.fasta has zero transcripts"
    exit 1
fi
echo "OK: tgg rule produced $n Trinity-GG transcripts"
