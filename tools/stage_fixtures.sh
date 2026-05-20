#!/usr/bin/env bash
# Pre-stage rule outputs from test_suite/fixtures/snapshot/ into
# test_suite/output/ so snakemake's mtime check skips those rules.
# Usage: bash tools/stage_fixtures.sh <subdir>...
#   e.g. bash tools/stage_fixtures.sh repeats rnaseq tgg proteins pasa golden augustus
set -euo pipefail
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
cd "$REPO"

SRC=test_suite/fixtures/snapshot
DST=test_suite/output

[[ -d "$SRC" ]] || {
    echo "FATAL: no snapshot at $SRC" >&2
    echo "       run tools/regen_snapshot.sh from a fresh clone, or" >&2
    echo "       tools/snapshot_dag_outputs.sh after a successful DAG run" >&2
    exit 1
}

mkdir -p "$DST"
for sub in "$@"; do
    [[ -d "$SRC/$sub" ]] || { echo "FATAL: snapshot has no $sub/ subdir" >&2; exit 1; }
    rm -rf "$DST/$sub"
    cp -a "$SRC/$sub" "$DST/$sub"
    # Snapshot files are produced AFTER source fixtures (which build_fixtures.sh
    # pins to 2020-01-01). cp -a preserves snapshot mtimes; no touch needed here.
    echo "staged $sub/"
done
