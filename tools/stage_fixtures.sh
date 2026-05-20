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
    # Set the staged outputs' mtime to "now". This is later than any
    # source-fixture mtime snakemake will encounter (test_suite/mini-*),
    # so snakemake's mtime check treats staged outputs as
    # up-to-date and skips re-running the staged rules. Downstream rules
    # whose outputs come from these staged inputs produce files with
    # mtime "now+ε", which is strictly later than the staged inputs --
    # snakemake's clock-skew detector does NOT fire (it requires input
    # mtime to be in the future relative to system time).
    #
    # The earlier `-d '2038-01-15'` trick clamped to the 32-bit time_t
    # cap to make staged outputs unconditionally newest, but it broke
    # any downstream rule that needed to produce new output: snakemake
    # would refuse to accept the new file because its mtime (now) was
    # older than the staged input (2038). See HANDOVER for the bug
    # history.
    #
    # `-h` is mandatory: snapshot/repeats/*/mini-genome.fasta is a
    # symlink to test_suite/mini-genome.fasta; without -h, the touch
    # would silently overwrite the source fixture's mtime.
    find "$DST/$sub" -exec touch -h {} +
    echo "staged $sub/"
done
