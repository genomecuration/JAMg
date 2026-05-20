#!/usr/bin/env bash
# Phase 3j integration test: workflow/rules/ogs.smk.
# Runs the full upstream chain through evm_run then ogs.smk end-to-end.
# Verifies the six canonical OGS products land under $OUTDIR/.

set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing"; exit 1; }
done

# Don't `rm -rf $OUT` — upstream stages are heavy. Snakemake re-runs only
# downstream-of-changes rules. Tests/runme.sh handles full-clean runs.
cd "$REPO"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until ogs_emit

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
