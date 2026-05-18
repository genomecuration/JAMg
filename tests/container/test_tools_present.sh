#!/usr/bin/env bash
# tests/container/test_tools_present.sh
# Verify that every tool in the JAMg SIF is (a) locatable via which and
# (b) executable with a version flag that exits 0 (or 1 -- some tools exit 1
# for --version but still print a version string, e.g. augustus).
#
# Usage: bash tests/container/test_tools_present.sh [path/to/jamg.sif]
#
# For each tool: (1) confirm `which` finds it, (2) run a version check
# tolerating the per-tool flag quirks captured in the lists below.

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
# - augustus: prints version, exits 1 with "no input given" hint.
# - RepeatMasker: prints "RepeatMasker version 4.2.3" then exits non-zero with
#   "Unknown option: version" warning and a usage banner.
# - exonerate: prints "exonerate from exonerate version 2.4.0" then exits 1.
ALLOW_NONZERO_VERSION="augustus RepeatMasker exonerate"

# Tools that use a single-dash -version (NCBI BLAST+ convention).
USE_SINGLE_DASH_VERSION="blastx"

# Tools that use -V instead of --version. (GSNAP accepts --version; only
# kept here if a tool truly needs -V.)
USE_CAPITAL_V=""

# Tools that use -v instead of --version.
USE_LOWERCASE_V=""

# Tools where version checking is skipped because the tool prints usage
# (not a version string) on any flag. Presence (which) is enough.
SKIP_VERSION="EVidenceModeler Launch_PASA_pipeline.pl TransDecoder.LongOrfs rmOutToGFF3.pl cdbfasta cdbyank trf ParaFly gmap_build"

TOOLS=(
    # Source-built into jamg.sif (containers/post-install.sh)
    augustus
    EVidenceModeler
    TransDecoder.LongOrfs
    cdbfasta
    cdbyank
    ParaFly
    # Slim bioconda set (containers/environment.yml)
    RepeatMasker
    RepeatModeler
    STAR
    minimap2
    gmap_build               # gmap suite (gsnap, gmap_build, etc.)
    gsnap
    samtools
    bcftools
    bedtools
    tabix
    bgzip
    gffread
    snakemake
    rmOutToGFF3.pl           # bioconda's RepeatMasker exposes util/ scripts
    # apt-installed in containers/jamg-base.def
    exonerate
    embossversion            # EMBOSS sentinel binary (the `emboss` package
                             # installs many; embossversion is the simplest)
    trf
    blastx                   # ncbi-blast+ from apt
    jellyfish
    salmon
    bowtie2
    # blat removed: not in debian trixie apt; v2 PASA uses gmap+minimap2 instead.
)

# Note: Trinity is in containers/trinity.sif (separate; tests/container/
# test_trinity_tools.sh). PASA is in containers/pasa.sif (separate;
# tests/container/test_pasa_tools.sh). They are NOT in this jamg.sif test
# because they're not in this SIF.

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
    elif echo " $USE_SINGLE_DASH_VERSION " | grep -q " $tool "; then
        VFLAG="-version"
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
