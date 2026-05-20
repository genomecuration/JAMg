#!/usr/bin/env bash
# Phase 6 end-to-end smoke. Runs the full DAG via Layout B (snakemake's
# Slurm executor fans each rule out as its own sbatch via the profile at
# workflow/profiles/slurm/). Invoke FROM THE HEAD NODE (lazebnik), NOT
# wrapped in sbatch: this script's snakemake process submits each rule
# and waits for completion. Total wall clock = DAG critical path, not
# the sum of all rule times.
#
# Intentionally does NOT pre-stage upstream outputs (except the
# committed genemark fixture, which is required because gmes_petap
# cannot converge on the 100 kb extract). All other rules run from
# scratch. This is the END-TO-END BASELINE, not a deployability smoke.
# For fast per-rule smokes that pre-stage upstream from the snapshot,
# see tests/rules/test_*.sh.
#
# Asserts the OGS gene count is within +/-1 of the baseline at
# test_suite/expected/OGS.gff3 (committed by `make smoke` once the
# pipeline is stable). Also asserts no legacy artefacts (zff, geneid,
# fathom) in the output.

set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
TRINITY_SIF="$REPO/containers/trinity.sif"
CFG="$REPO/test_suite/mini-config.yaml"
EXPECTED="$REPO/test_suite/expected/OGS.gff3"

for sif in "$JAMG_SIF" "$PASA_SIF" "$TRINITY_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing — run 'make sifs' first"; exit 1; }
done

GM_FIXTURE="$REPO/test_suite/fixtures/genemark.gff3"
if [[ ! -s "$GM_FIXTURE" ]]; then
    echo "FATAL: genemark fixture missing at $GM_FIXTURE" >&2
    echo "Run 'bash tests/rules/test_genemark.sh' first to produce it." >&2
    exit 1
fi

cd "$REPO"
rm -rf test_suite/output

# Pre-stage genemark fixture (see tests/rules/test_evm.sh for rationale).
# The full DAG on the 100 kb default fixture cannot run gmes_petap --ES
# (needs >=1 Mb); the fixture is a 100 kb subset of a 1 Mb prediction.
GM_FIXTURE_GTF="$REPO/test_suite/fixtures/genemark.gtf"
[[ -s "$GM_FIXTURE_GTF" ]] || {
    echo "FATAL: genemark .gtf fixture missing at $GM_FIXTURE_GTF" >&2
    echo "Run 'bash tests/rules/test_genemark.sh' first to (re)produce both fixtures." >&2
    exit 1
}
mkdir -p "$REPO/test_suite/output/genemark"
cp "$GM_FIXTURE" "$REPO/test_suite/output/genemark/genemark.gff3"
cp "$GM_FIXTURE_GTF" "$REPO/test_suite/output/genemark/genemark.gtf"
: > "$REPO/test_suite/output/genemark/.preflight.ok"

# Layout B: snakemake-slurm submits each rule as its own sbatch via the
# profile at workflow/profiles/slurm/. Per-rule threads/mem_mb map to
# --cpus-per-task / --mem on the sbatch line. This script runs on the
# head node and waits for all per-rule jobs to complete.
pixi run "$REPO/bin/jamg" run --config "$CFG" --executor slurm --jobs 32

test -s "$REPO/test_suite/output/OGS.gff3" || { echo "no OGS.gff3"; exit 1; }

# Gene-count drift check (±1) against the committed baseline.
if [[ -f "$EXPECTED" ]]; then
    expected_genes=$(grep -cP '\tgene\t' "$EXPECTED")
    actual_genes=$(grep -cP '\tgene\t' "$REPO/test_suite/output/OGS.gff3")
    if (( actual_genes < expected_genes - 1 || actual_genes > expected_genes + 1 )); then
        echo "FATAL: OGS gene count drift: got $actual_genes, expected ~$expected_genes"
        exit 1
    fi
    echo "OK: gene count $actual_genes within ±1 of baseline $expected_genes"
else
    echo "NOTE: no baseline at $EXPECTED; recording current gene count for visibility"
    actual_genes=$(grep -cP '\tgene\t' "$REPO/test_suite/output/OGS.gff3")
    echo "  $actual_genes genes in $REPO/test_suite/output/OGS.gff3"
fi

# No legacy validator / predictor artefacts in the workflow output.
if find "$REPO/test_suite/output" \
        -name '*.zff' -o -name '*.geneid*' -o -name '*.fathom*' 2>/dev/null \
        | grep -q .; then
    echo "FATAL: legacy validator/predictor artefacts found in workflow output"
    find "$REPO/test_suite/output" \
         -name '*.zff' -o -name '*.geneid*' -o -name '*.fathom*' >&2
    exit 1
fi

echo "OK: full DAG produced OGS.gff3 with no legacy artefacts"
