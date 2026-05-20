#!/usr/bin/env bash
# Phase 6 end-to-end smoke. Runs the full DAG on the mini-fixtures and
# asserts the OGS gene count is within +/-1 of the baseline at
# test_suite/expected/OGS.gff3 (committed by `make smoke` once the
# pipeline is stable).
#
# Also asserts no legacy artefacts (zff, geneid, fathom) in the output.
#
# This is the END-TO-END BASELINE: it intentionally does NOT pre-stage
# upstream outputs (only genemark, which needs the 1 Mb fixture). All
# other rules run from scratch on the 100 kb default fixture. Runtime
# is therefore long (15-30 min on Layout A; ~5-10 min on Layout B once
# Task 10 lands). This test is the gold standard, not a deployability
# smoke. For fast per-rule smokes that pre-stage upstream from the
# snapshot, see tests/rules/test_*.sh.

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
touch -d '2038-01-15' "$REPO/test_suite/output/genemark/genemark.gff3" \
                     "$REPO/test_suite/output/genemark/genemark.gtf" \
                     "$REPO/test_suite/output/genemark/.preflight.ok"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --mtime-only

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
