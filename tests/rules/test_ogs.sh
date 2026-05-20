#!/usr/bin/env bash
# Phase 3j integration test: workflow/rules/ogs.smk.
# Pre-stages every upstream rule's output from test_suite/fixtures/snapshot/
# (including evm/EVM.gff3) then runs ONLY the OGS rules. Target runtime:
# under 2 minutes.
#
# Prerequisites:
#   - tools/build_fixtures.sh has been run
#   - tools/regen_snapshot.sh has been run (snapshot/ tree exists with evm/)
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output"
SNAP="$REPO/test_suite/fixtures/snapshot"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing"; exit 1; }
done
[[ -d "$SNAP" ]] || {
    echo "FATAL: no snapshot at $SNAP" >&2
    echo "Run tools/regen_snapshot.sh first." >&2
    exit 1
}
[[ -s "$SNAP/evm/EVM.gff3" ]] || {
    echo "FATAL: snapshot's evm/EVM.gff3 is empty; cannot test OGS without it" >&2
    exit 1
}

cd "$REPO"
rm -rf test_suite/output .snakemake

# Stage every upstream rule's outputs (including EVM, which is the immediate
# upstream of OGS).
bash tools/stage_fixtures.sh repeats rnaseq tgg proteins pasa golden augustus evm

# Genemark fixture (separate from the snapshot; subset of a 1 Mb run).
GM_OUT="$REPO/test_suite/output/genemark"
mkdir -p "$GM_OUT"
cp "$REPO/test_suite/fixtures/genemark.gff3" "$GM_OUT/genemark.gff3"
cp "$REPO/test_suite/fixtures/genemark.gtf"  "$GM_OUT/genemark.gtf"
: > "$GM_OUT/.preflight.ok"
touch -d '2038-01-15' "$GM_OUT/genemark.gff3" "$GM_OUT/genemark.gtf" "$GM_OUT/.preflight.ok"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" \
    --until ogs_emit

fail=0
for product in OGS.gff3 OGS.mRNA.fasta OGS.CDS.fasta OGS.pep.fasta OGS.gtf OGS.bed; do
    if [[ ! -s "$OUT/$product" ]]; then
        echo "FATAL: $product missing or empty"
        fail=1
    fi
done
(( fail == 0 )) || exit 1

n_genes=$(awk -F'\t' '$3 == "gene"' "$OUT/OGS.gff3" | wc -l)
n_pep=$(grep -c '^>' "$OUT/OGS.pep.fasta")
(( n_genes >= 1 )) || { echo "FATAL: OGS.gff3 has 0 gene rows"; exit 1; }
(( n_pep   >= 1 )) || { echo "FATAL: OGS.pep.fasta has 0 sequences"; exit 1; }
echo "OK: OGS produced $n_genes gene rows + $n_pep peptide sequences"
