#!/usr/bin/env bash
# Phase 3f integration test for workflow/rules/genemark.smk.
# Runs the preflight (extracts gmes tarball, ungzips key, ABI-checks
# gmes_petap.pl --help) and the main genemark rule (gmes_petap.pl --ET
# on the softmasked genome + RNA-seq intron evidence). Produces a
# canonical genemark.gff3 consumed downstream by EVM.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/genemark"

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
    "$REPO/bin/jamg" run --config "$CFG" --cores 4 --until genemark

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
