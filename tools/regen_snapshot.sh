#!/usr/bin/env bash
# Regenerate test_suite/fixtures/snapshot/ from scratch.
# Required before any per-rule test (tests/rules/test_*.sh, except
# test_genemark) can be run from a fresh clone; the snapshot itself is
# gitignored to keep the repo lightweight, and the per-rule tests rely on
# its presence to pre-stage upstream rule outputs.
set -euo pipefail
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
cd "$REPO"

# 1. Build the polyester-driven RNA-seq fixtures + genome extracts if missing.
[[ -s test_suite/mini-genome.fasta && -s test_suite/mini-rnaseq.bam ]] \
    || bash tools/build_fixtures.sh

# 2. Clean state for the seed full-DAG run.
rm -rf test_suite/output test_suite/output-1mb .snakemake

# 3. Pre-stage the committed genemark fixture (the 1 Mb genemark prediction
# subset to the 100 kb window) so the full DAG on the 100 kb fixture can
# skip genemark. test_genemark.sh regenerates this fixture from the 1 Mb
# config when needed.
#
# Stage to genemark.staged.{gtf,gff3} (NOT the declared outputs). Snakemake
# removes its declared output files before invoking the rule, so staging
# under the canonical names is destroyed. The genemark rule's shell guard
# (workflow/rules/genemark.smk) cps the .staged. files into the declared
# positions and exits 0 when staging is detected.
mkdir -p test_suite/output/genemark
cp test_suite/fixtures/genemark.gff3 test_suite/output/genemark/genemark.staged.gff3
cp test_suite/fixtures/genemark.gtf  test_suite/output/genemark/genemark.staged.gtf

# 4. Run the seed full DAG (everything through ogs_emit), producing the
# outputs Task 3 snapshots.
pixi run bin/jamg run --config test_suite/mini-config.yaml --cores "${SLURM_CPUS_PER_TASK:-20}" \
    --until ogs_emit

# 5. Snapshot.
bash tools/snapshot_dag_outputs.sh

echo "OK: snapshot regenerated at test_suite/fixtures/snapshot/"
echo "    per-rule tests can now run via tools/stage_fixtures.sh (Task 4)"
