#!/usr/bin/env bash
# Phase 3f integration test for workflow/rules/genemark.smk.
# Runs the preflight (extracts gmes tarball, ungzips key, ABI-checks
# gmes_petap.pl --help) and the main genemark rule (gmes_petap.pl --ET
# on the softmasked genome + RNA-seq intron evidence). Produces a
# canonical genemark.gff3 consumed downstream by EVM.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"
# 1 Mb genome (gmes_petap --ET needs enough transcript coverage to converge).
# Every other rule test uses the 100 kb default mini-config.yaml.
# Distinct outdir (test_suite/output-1mb) so this test can run in parallel
# with test_evm / test_ogs / test_full_dag (which use test_suite/output)
# without colliding on shared output paths.
CFG="$REPO/test_suite/mini-config-1mb.yaml"
OUT="$REPO/test_suite/output-1mb/genemark"
FIXTURE_DIR="$REPO/test_suite/fixtures"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: $SIF not found. Run 'make sif RM_LIB_HOST=...' first." >&2
    exit 1
fi

for f in ~/software/gmes_linux_64_4.tar.gz ~/software/gm_key_64.gz; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: $f not found. The GeneMark archives are required for Phase 3f." >&2
        exit 1
    fi
done

rm -rf "$OUT"

apptainer exec --bind "$REPO":"$REPO" --bind "$HOME":"$HOME" "$SIF" \
    "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until genemark

test -e "$OUT/.preflight.ok"     || { echo "preflight sentinel missing";          exit 1; }
test -s "$OUT/genemark.gff3"     || { echo "genemark.gff3 missing/empty";         exit 1; }
test -s "$OUT/genemark.gtf"      || { echo "genemark.gtf missing/empty";          exit 1; }
test -x "$OUT/extracted/gmes_linux_64_4/gmes_petap.pl" \
    || { echo "extracted gmes_petap.pl not executable"; exit 1; }
test -s "$OUT/home/.gm_key"      || { echo "ungzipped .gm_key missing/empty";     exit 1; }

# At least one gene prediction
n=$(grep -cE '^[^#]\S' "$OUT/genemark.gff3"; true)
if (( n < 1 )); then
    echo "genemark.gff3 has zero non-comment features"
    exit 1
fi
echo "OK: genemark rule produced $n features in genemark.gff3"

# Emit the 100 kb subset as a fixture for tests/rules/test_evm.sh,
# tests/rules/test_ogs.sh, and tests/e2e/test_full_dag.sh. Those tests
# run on the 100 kb default fixture (X_mini coords 1..100000); they
# pre-stage this file into test_suite/output/genemark/ so snakemake
# skips the slow genemark rule via mtime comparison.
#
# Filter is $5<=100000 (end-coord must fit the 100 kb window). A gene
# spanning 99000-110000 is dropped entirely; its individual child rows
# (exon/CDS) with $5<=100000 survive as orphans, but sort_gff3.pl drops
# orphans silently downstream, so the fixture is structurally clean by
# the time EVM consumes it. Filtering by start coord ($4<=100000) would
# admit invalid rows that overshoot the 100 kb genome boundary.
mkdir -p "$FIXTURE_DIR"
awk -F'\t' '/^#/ {print; next} $1=="X_mini" && $4>=1 && $5<=100000 {print}' \
    "$OUT/genemark.gff3" > "$FIXTURE_DIR/genemark.gff3"
awk -F'\t' '/^#/ {print; next} $1=="X_mini" && $4>=1 && $5<=100000 {print}' \
    "$OUT/genemark.gtf"  > "$FIXTURE_DIR/genemark.gtf"
subset_n=$(grep -cE '^[^#]\S' "$FIXTURE_DIR/genemark.gff3" || true)
echo "OK: emitted 100 kb subset fixture ($subset_n features) -> $FIXTURE_DIR/genemark.gff3"
