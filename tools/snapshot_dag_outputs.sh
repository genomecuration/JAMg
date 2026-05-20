#!/usr/bin/env bash
# Copy a successful test_suite/output/ tree into test_suite/fixtures/snapshot/
# so per-rule tests can pre-stage upstream outputs and skip those rules via
# snakemake's --rerun-triggers mtime. Run after a verified successful
# full-DAG run (e.g. test_evm passing or Task 2 seed completing).
# The snapshot is gitignored; regenerate via tools/regen_snapshot.sh.
set -euo pipefail
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
cd "$REPO"

SRC=test_suite/output
DST=test_suite/fixtures/snapshot

[[ -d "$SRC" ]] || { echo "FATAL: $SRC missing; run the full DAG first" >&2; exit 1; }
mkdir -p "$DST"

# Required per-rule directories. Tests/rules/test_*.sh depend on these.
REQUIRED=(repeats rnaseq tgg proteins pasa golden augustus)
# Optional: evm and ogs depend on EVidenceModeler running cleanly, which
# fails on the 100 kb fixture due to overlapping PASA predictions; see
# docs/limitations.md (or HANDOVER.md). genemark is pre-staged from a
# separate committed fixture and not part of this snapshot.
OPTIONAL=(evm genemark)

for sub in "${REQUIRED[@]}"; do
    [[ -d "$SRC/$sub" ]] || { echo "FATAL: $SRC/$sub missing" >&2; exit 1; }
    rm -rf "$DST/$sub"
    cp -a "$SRC/$sub" "$DST/$sub"
done

for sub in "${OPTIONAL[@]}"; do
    if [[ -d "$SRC/$sub" ]]; then
        rm -rf "$DST/$sub"
        cp -a "$SRC/$sub" "$DST/$sub"
    else
        echo "NOTE: $SRC/$sub missing; snapshot will be partial (test_ogs will be blocked)"
    fi
done

# Manifest for future staleness detection.
(cd "$DST" && find . -type f -not -name 'MANIFEST.sha256' -exec sha256sum {} +) \
    | sort -k2 > "$DST/MANIFEST.sha256"

# Explicit non-empty check on EVM output that test_ogs.sh consumes (only
# enforced when the seed run actually included EVM; partial snapshots without
# EVM are valid for test_evm but block test_ogs).
if [[ -d "$DST/evm" ]]; then
    [[ -s "$DST/evm/EVM.gff3" ]] \
        || { echo "FATAL: snapshot's evm/EVM.gff3 is empty" >&2; exit 1; }
fi

# PASA sqlite portability check. The pass1 sqlite is what ogs_pasa_compare_*
# consumes; a corrupt/non-portable snapshot here breaks test_ogs.sh silently.
if [[ -s "$DST/pasa/pasa.sqlite.pass1.bz2" ]]; then
    TMPDB=$(mktemp --suffix .sqlite)
    bunzip2 -c "$DST/pasa/pasa.sqlite.pass1.bz2" > "$TMPDB"
    if ! apptainer exec containers/pasa.sif sqlite3 "$TMPDB" \
            "SELECT COUNT(*) FROM status;" >/dev/null 2>&1; then
        rm -f "$TMPDB"
        echo "FATAL: pasa.sqlite.pass1.bz2 is not portable (no 'status' table or corrupt)" >&2
        echo "       this snapshot WILL break test_ogs.sh's ogs_pasa_compare_* rules" >&2
        exit 1
    fi
    rm -f "$TMPDB"
    echo "OK: pasa sqlite portability check passed"
fi

echo "OK: snapshot at $DST"
du -sh "$DST"/* 2>/dev/null
echo "Total: $(du -sh "$DST" | awk '{print $1}')"
