#!/usr/bin/env bash
# Phase 3d integration test for workflow/rules/pasa.smk.
# Runs the §3d PASA chain: setup_db -> align -> compare_transdecoder -> hints.
# Asserts the post-PASA artefacts that EVM (§3i) consumes are present:
# transdecoder.genome.gff3, pasa_assemblies.gff3, polyAsites.fasta + hints.
#
# pasa_compare_load_1/_2 are NOT tested here — they were moved to plan §3j
# (post-EVM) because PASA's `-A -L --annots` mode requires a pre-existing
# DIFFERENT annotation to update, which only exists post-EVM (per PASA
# upstream wiki PASA_genome_annotation.md). See pasa.smk's NOTE block where
# pasa_compare_load_1 used to be defined.
set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
JAMG_SIF="$REPO/containers/jamg.sif"
PASA_SIF="$REPO/containers/pasa.sif"
CFG="$REPO/test_suite/mini-config.yaml"
OUT="$REPO/test_suite/output/pasa"

for sif in "$JAMG_SIF" "$PASA_SIF"; do
    if [[ ! -f "$sif" ]]; then
        echo "ERROR: $sif not found. Run 'make sifs' first." >&2
        exit 1
    fi
done

rm -rf "$OUT"

cd "$REPO"
pixi run "$REPO/bin/jamg" run --config "$CFG" --cores "${SLURM_CPUS_PER_TASK:-20}" --until pasa_hints

test -s "$OUT/pasa_assemblies.gff3"             || { echo "pasa_assemblies.gff3 missing/empty"; exit 1; }
test -s "$OUT/pasa.transdecoder.genome.gff3"    || { echo "transdecoder.genome.gff3 missing/empty"; exit 1; }
test -e "$OUT/polyAsites.fasta"                 || { echo "polyAsites.fasta missing"; exit 1; }
test -e "$OUT/pasa.sqlite.align.bz2"            || { echo "align snapshot missing"; exit 1; }
test -e "$OUT/pasa.sqlite.pass1.bz2"            || { echo "pass1.bz2 snapshot missing"; exit 1; }
test -e "$OUT/polyAsites.hints"                 || { echo "polyA hints missing"; exit 1; }
test -e "$OUT/pasa_assemblies.gff3.hints"       || { echo "assembly hints missing"; exit 1; }

assemblies=$(grep -cvE '^(#|$)' "$OUT/pasa_assemblies.gff3"; true)
echo "OK: pasa chain produced $assemblies assembly features + transdecoder.genome.gff3 + hints"
