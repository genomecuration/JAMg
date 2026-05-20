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
GM_OUT="$REPO/test_suite/output/genemark"
GM_FIXTURE="$REPO/test_suite/fixtures/genemark.gff3"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing"; exit 1; }
done

# Pre-stage genemark fixture (see comment in tests/rules/test_evm.sh for
# rationale). gmes_petap --ES needs >=1 Mb; the 100 kb default cannot
# run it, so we subset a 1 Mb genemark prediction into the 100 kb output
# dir and let snakemake skip the rule via mtime.
if [[ ! -s "$GM_FIXTURE" ]]; then
    echo "FATAL: genemark fixture missing at $GM_FIXTURE" >&2
    echo "Run 'bash tests/rules/test_genemark.sh' first to produce it." >&2
    exit 1
fi
cd "$REPO"
# test_ogs chains from test_evm; assert the precondition explicitly rather
# than silently consuming a corrupted output dir from a cancelled test_evm.
[[ -s "$OUT/evm/EVM.gff3" ]] || {
    echo "FATAL: $OUT/evm/EVM.gff3 missing or empty" >&2
    echo "test_ogs.sh requires test_evm.sh to have completed cleanly first." >&2
    exit 1
}
GM_FIXTURE_GTF="$REPO/test_suite/fixtures/genemark.gtf"
[[ -s "$GM_FIXTURE_GTF" ]] || {
    echo "FATAL: genemark .gtf fixture missing at $GM_FIXTURE_GTF" >&2
    echo "Run 'bash tests/rules/test_genemark.sh' first to (re)produce both fixtures." >&2
    exit 1
}
mkdir -p "$GM_OUT"
cp "$GM_FIXTURE" "$GM_OUT/genemark.gff3"
cp "$GM_FIXTURE_GTF" "$GM_OUT/genemark.gtf"
: > "$GM_OUT/.preflight.ok"
touch -d '2038-01-15' "$GM_OUT/genemark.gff3" "$GM_OUT/genemark.gtf" "$GM_OUT/.preflight.ok"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until ogs_emit --mtime-only

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
