#!/usr/bin/env bash
# Smoke test for bin/optimize_extrinsic_augustus.pl.
# Asserts that --help exits and emits a Usage block.
#
# (NOTE — script bug, not test scope: bin/optimize_extrinsic_augustus_merge_cfg.pl
# `die unless $files[1] && -s $files[1]` at line 8 with NO --help / Pod::Usage,
# so it cannot satisfy a "--help runs cleanly" assertion. Verifying only the
# main optimize_extrinsic_augustus.pl per spec §3h.1 — that script DOES emit
# Usage via pod2usage on bad opts.)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SIF="$REPO/containers/jamg.sif"
SCRIPT="$REPO/bin/optimize_extrinsic_augustus.pl"

if [[ ! -e "$SIF" ]]; then
    echo "SKIP: containers/jamg.sif not built" >&2
    exit 0
fi

out=$(apptainer exec --bind "$REPO":"$REPO" "$SIF" perl "$SCRIPT" --help 2>&1 || true)
if ! echo "$out" | grep -q "Usage"; then
    echo "FAIL: --help did not emit Usage. Got:" >&2
    echo "$out" | head -20 >&2
    exit 1
fi
if ! echo "$out" | grep -q -- "--species"; then
    echo "FAIL: Usage did not document --species" >&2
    exit 1
fi
echo "OK"
