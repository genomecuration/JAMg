#!/usr/bin/env bash
# tests/container/test_trinity_tools.sh
# Verify the trinity.sif (upstream Docker image trinityrnaseq/trinityrnaseq:2.15.2)
# has the runtime entrypoints the workflow's rules/tgg.smk will call.
#
# Usage: bash tests/container/test_trinity_tools.sh [path/to/trinity.sif]

set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="${1:-"$REPO_ROOT/containers/trinity.sif"}"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: trinity SIF does not exist (expected at $SIF)" >&2
    echo "Run 'make sifs' (note plural) to build it." >&2
    exit 1
fi

TOOLS=(Trinity jellyfish salmon bowtie2 samtools)
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
