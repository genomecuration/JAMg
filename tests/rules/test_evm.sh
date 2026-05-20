#!/usr/bin/env bash
# Phase 3i integration test: workflow/rules/evm.smk.
# Runs the full upstream chain (repeats + rnaseq + tgg + pasa + proteins +
# genemark + golden + augustus) then evm.smk end-to-end through evm_run.
# Verifies EVM.gff3 is produced and contains at least one EVidenceModeler
# feature row.

set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/evm"
GM_OUT="$REPO/test_suite/output/genemark"
GM_FIXTURE="$REPO/test_suite/fixtures/genemark.gff3"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    [[ -f "$sif" ]] || { echo "ERROR: $sif missing"; exit 1; }
done

# Pre-stage the genemark output (subset of a 1 Mb run) into the workflow
# output dir so snakemake skips the slow genemark rule on the 100 kb
# default fixture. gmes_petap --ES self-training needs >=1 Mb to
# converge; subsetting a 1 Mb prediction to the 100 kb window keeps
# this test fast (~1-2 min vs ~15+ min) without exercising genemark.
# Run tests/rules/test_genemark.sh first to (re)generate the fixture.
if [[ ! -s "$GM_FIXTURE" ]]; then
    echo "FATAL: genemark fixture missing at $GM_FIXTURE" >&2
    echo "Run 'bash tests/rules/test_genemark.sh' first to produce it." >&2
    exit 1
fi
cd "$REPO"
# Wipe the whole output tree (not just $OUT) so stale intermediate files
# from a prior cancelled run or a different config (e.g. test_genemark on
# the 1 Mb fixture) cannot leak through snakemake's mtime-only rerun
# check. RepeatMasker caches alignments to <input>.cat.gz and reuses
# them silently — this killed JID 10020.
rm -rf test_suite/output

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
# 2038-01-15 is the 32-bit time_t boundary on this filesystem; later
# dates silently clamp. 2038 is still far-future enough for snakemake.
touch -d '2038-01-15' "$GM_OUT/genemark.gff3" "$GM_OUT/genemark.gtf" "$GM_OUT/.preflight.ok"

pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until evm_run --mtime-only

test -s "$OUT/EVM.gff3" || { echo "EVM.gff3 missing or empty"; exit 1; }

# At least one EVidenceModeler-emitted feature row (column 2 source 'EVM').
n=$(awk -F'\t' '$2 == "EVM" && $3 == "gene"' "$OUT/EVM.gff3" | wc -l)
if (( n < 1 )); then
    echo "FATAL: EVM.gff3 contains zero gene rows with source=EVM"
    head -20 "$OUT/EVM.gff3" >&2
    exit 1
fi
echo "OK: EVM.gff3 produced $n EVM gene rows"
