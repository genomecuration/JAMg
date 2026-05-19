#!/usr/bin/env bash
# Phase 3h integration test for workflow/rules/augustus.smk.
# Asserts that `bin/jamg run --until augustus_preflight` materialises
# the preflight sentinel, which means:
#   (1) `species: drosophila` resolves in augustus's species library,
#   (2) `workflow/config/extrinsic.cfg` [SOURCES] line is the v2 token set,
#   (3) extrinsic.cfg parses cleanly under Augustus 3.5 against /dev/null.
#
# The full augustus rule's DAG inputs include pasa hints + golden hints
# (per plan §3h), so an end-to-end `--until augustus` run is blocked
# until pasa_align lands. This test isolates the preflight portion,
# which has no pasa/golden dependency.
#
# TODO(post-pasa-fix): extend this test to also exercise `--until augustus`
# (end-to-end). The augustus rule body (workflow/rules/augustus.smk lines
# ~80-120) is currently not covered by any integration test — the
# code-review on commit-X flagged a `cd` + relative-output redirect bug
# that the preflight-only test could not catch, and similar pattern bugs
# in the rule body will go undetected until end-to-end coverage lands.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/augustus"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: $SIF not found. Run 'make sif' first." >&2
    exit 1
fi

rm -rf "$OUT"

apptainer exec --bind "$REPO":"$REPO" "$SIF" \
    "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until augustus_preflight

test -e "$OUT/.preflight.ok" || { echo "preflight sentinel missing"; exit 1; }

echo "OK: augustus_preflight succeeded (species + extrinsic.cfg checks all passed)"
