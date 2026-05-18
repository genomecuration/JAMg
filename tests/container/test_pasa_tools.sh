#!/usr/bin/env bash
# tests/container/test_pasa_tools.sh
# Verify the pasa.sif (upstream Docker image pasapipeline/pasapipeline:2.5.3)
# has the runtime entrypoints the workflow's rules/pasa.smk will call.
#
# Usage: bash tests/container/test_pasa_tools.sh [path/to/pasa.sif]

set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="${1:-"$REPO_ROOT/containers/pasa.sif"}"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: pasa SIF does not exist (expected at $SIF)" >&2
    echo "Run 'make sifs' (note plural) to build it." >&2
    exit 1
fi

TOOLS=(Launch_PASA_pipeline.pl seqclean accession_extractor.pl)
FAILURES=0
for tool in "${TOOLS[@]}"; do
    if apptainer exec "$SIF" which "$tool" >/dev/null 2>&1; then
        echo "OK: $tool"
    else
        echo "MISSING: $tool"
        FAILURES=$(( FAILURES + 1 ))
    fi
done

if [[ "$FAILURES" -gt 0 ]]; then
    echo "RESULT: $FAILURES tool(s) MISSING in $SIF" >&2
    exit 1
fi

echo "OK: ${#TOOLS[@]} tools present in $SIF"
