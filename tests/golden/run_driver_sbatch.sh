#!/usr/bin/env bash
# Phase 4 driver integration helper. Runs bin/prepare_golden_genes.pl inside
# jamg.sif against the mini fixtures. Used by tests/golden/test_prepare_golden_genes.t
# when JAMG_GOLDEN_SBATCH is set, OR submitted directly via sbatch.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="$REPO/containers/jamg.sif"
GENOME="$REPO/test_suite/mini-genome.fasta"
MRNA="$REPO/test_suite/mini-transcripts.fasta"
SOFT="$REPO/test_suite/output/repeats/mini-genome.fasta.softmasked"
OUT="${1:-$TMP/golden-driver-test-$$}"

mkdir -p "$OUT"
cd "$REPO"

[[ -s "$SIF"   ]] || { echo "ERROR: $SIF missing"   >&2; exit 1; }
[[ -s "$GENOME" ]] || { echo "ERROR: $GENOME missing" >&2; exit 1; }
[[ -s "$MRNA"  ]] || { echo "ERROR: $MRNA missing"  >&2; exit 1; }
[[ -s "$SOFT"  ]] || { echo "ERROR: $SOFT missing — run tests/rules/test_repeats.sh first" >&2; exit 1; }

THREADS="${SLURM_CPUS_PER_TASK:-4}"
echo "Running prepare_golden_genes.pl with -threads $THREADS, outdir=$OUT"

exec apptainer exec --bind "$REPO":"$REPO" --bind "$OUT":"$OUT" "$SIF" \
    perl "$REPO/bin/prepare_golden_genes.pl" \
        -genome "$GENOME" -mrna "$MRNA" -softmasked "$SOFT" \
        -augustus /opt/jamg/share/Augustus \
        -outdir "$OUT" -threads "$THREADS"
