#!/usr/bin/env bash
# tests/container/test_tools_present.sh
# Verify that every tool in the JAMg SIF is (a) locatable via which and
# (b) executable with a version flag that exits 0 (or 1 -- some tools exit 1
# for --version but still print a version string, e.g. augustus).
#
# Usage: bash tests/container/test_tools_present.sh [path/to/jamg.sif]
#
# Review-#1 N31: also runs a version check per tool, with allowances for
# tools that use -V or -v instead of --version.

set -euo pipefail

# Resolve SIF path relative to repo root regardless of CWD.
REPO_ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
SIF="${1:-"$REPO_ROOT/containers/jamg.sif"}"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: containers/jamg.sif does not exist (expected at $SIF)" >&2
    echo "Run 'make sif' first." >&2
    exit 1
fi

# Tools that accept --version but may exit non-zero (still print version).
# We allow exit codes 0 and 1 for these.
ALLOW_NONZERO_VERSION="augustus"

# Tools that use -V instead of --version (GMAP suite).
USE_CAPITAL_V="gsnap gmap_build"

# Tools that use -v instead of --version.
USE_LOWERCASE_V=""

# Tools where version checking is skipped -- scripts that accept no flags
# and print usage on any unknown flag.  Presence (which) is enough.
# trf prints its own help banner on no args and exits non-zero on --version,
# so presence is enough (code-review N3 addition).
SKIP_VERSION="EVidenceModeler Launch_PASA_pipeline.pl TransDecoder.LongOrfs rmOutToGFF3.pl cdbfasta cdbyank trf"

TOOLS=(
    augustus
    EVidenceModeler
    Launch_PASA_pipeline.pl
    TransDecoder.LongOrfs
    RepeatMasker
    RepeatModeler
    trf                       # code-review N3: RepeatMasker configure passes
                              # -trf_prgm=$PREFIX/bin/trf so trf must exist
    Trinity
    STAR
    minimap2
    gsnap
    gmap_build
    blastx
    samtools
    bcftools
    bedtools
    tabix
    bgzip
    gffread
    cdbfasta
    cdbyank
    snakemake
    exonerate
    ParaFly
    rmOutToGFF3.pl
)

FAILURES=0

for tool in "${TOOLS[@]}"; do
    # --- presence check ---
    if ! apptainer exec "$SIF" which "$tool" >/dev/null 2>&1; then
        echo "MISSING (which): $tool"
        FAILURES=$(( FAILURES + 1 ))
        continue
    fi

    # --- version check (skip list) ---
    if echo " $SKIP_VERSION " | grep -q " $tool "; then
        echo "OK (which only): $tool"
        continue
    fi

    # Determine version flag.
    if echo " $USE_CAPITAL_V " | grep -q " $tool "; then
        VFLAG="-V"
    elif echo " $USE_LOWERCASE_V " | grep -q " $tool "; then
        VFLAG="-v"
    else
        VFLAG="--version"
    fi

    # Run version check; capture exit code.
    if echo " $ALLOW_NONZERO_VERSION " | grep -q " $tool "; then
        apptainer exec "$SIF" "$tool" "$VFLAG" >/dev/null 2>&1 || true
        echo "OK (version/$VFLAG, non-zero allowed): $tool"
    else
        if apptainer exec "$SIF" "$tool" "$VFLAG" >/dev/null 2>&1; then
            echo "OK (version/$VFLAG): $tool"
        else
            echo "FAILED (version/$VFLAG exit non-zero): $tool"
            FAILURES=$(( FAILURES + 1 ))
        fi
    fi
done

if [[ "$FAILURES" -gt 0 ]]; then
    echo ""
    echo "RESULT: $FAILURES tool(s) FAILED" >&2
    exit 1
fi

echo ""
echo "OK: ${#TOOLS[@]} tools present in $SIF"
