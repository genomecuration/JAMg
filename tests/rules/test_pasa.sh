#!/usr/bin/env bash
# Phase 3d integration test for workflow/rules/pasa.smk.
# Runs the full PASA pre-EVM chain: setup_db -> align -> compare_transdecoder
# -> compare_load_1 -> compare_load_2 -> hints. Asserts pass3.bz2 exists
# (consumed by 3j) and pasa_assemblies.gff3 is non-empty.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/pasa"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    if [[ ! -f "$sif" ]]; then
        echo "ERROR: $sif not found. Run 'make sifs RM_LIB_HOST=...' first." >&2
        exit 1
    fi
done

rm -rf "$OUT"

cd "$REPO"
pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until pasa_hints

test -s "$OUT/pasa_assemblies.gff3"                          || { echo "pasa_assemblies.gff3 missing/empty";  exit 1; }
test -e "$OUT/pasa.transdecoder.genome.gff3"                 || { echo "transdecoder GFF missing";           exit 1; }
test -e "$OUT/polyAsites.fasta"                              || { echo "polyAsites.fasta missing";           exit 1; }
test -e "$OUT/pasa.sqlite.align.bz2"                         || { echo "align snapshot missing";             exit 1; }
test -e "$OUT/pasa.sqlite.pass1.bz2"                         || { echo "pass1.bz2 snapshot missing";        exit 1; }
test -e "$OUT/pasa.sqlite.pass2.bz2"                         || { echo "pass2.bz2 snapshot missing";        exit 1; }
test -e "$OUT/pasa.sqlite.pass3.bz2"                         || { echo "pass3.bz2 snapshot missing";        exit 1; }
test -e "$OUT/gene_structures_post_PASA_updates.round1.gff3" || { echo "round1 gene structures missing";    exit 1; }
test -e "$OUT/gene_structures_post_PASA_updates.round2.gff3" || { echo "round2 gene structures missing";    exit 1; }
test -e "$OUT/polyAsites.hints"                              || { echo "polyA hints missing";               exit 1; }
test -e "$OUT/pasa_assemblies.gff3.hints"                    || { echo "assembly hints missing";            exit 1; }

assemblies=$(grep -cvE '^(#|$)' "$OUT/pasa_assemblies.gff3"; true)
echo "OK: pasa rule chain produced $assemblies assembly features and the three snapshot files"
