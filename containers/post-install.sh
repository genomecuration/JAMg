#!/usr/bin/env bash
# Source-build the analytic tools that are NOT in apt and NOT in the slim
# bioconda set: Augustus, EVidenceModeler, TransDecoder, cdbfasta, ParaFly.
#
# Tools handled elsewhere:
#   - exonerate, emboss, trf, ncbi-blast+: apt (jamg-base.sif).
#   - minimap2, gmap, samtools, bcftools, bedtools, gffread, htslib, star,
#     repeatmasker, repeatmodeler: bioconda (containers/environment.yml).
#   - Trinity, PASA: separate SIFs (trinity.def, pasa.def, upstream Docker).
#
# Build environment (set by the bash subshell wrapper in jamg.def):
#   - micromamba base env active.
#   - PATH=/opt/conda/bin:... (compiler + autotools chain).
#   - CC, CXX, AR, RANLIB exported by conda activate.d.
#
# Authoritative pin ledger: containers/source-pins.toml. The version
# variables below must agree with that file.

set -euo pipefail

PREFIX=/opt/conda
JAMG=/opt/jamg
SRC=/opt/src
mkdir -p "$SRC" "$JAMG/share"
cd "$SRC"

# -------------------------------------------------------------------------
# Augustus (gene prediction; consumed by workflow/rules/augustus.smk).
# Pinned to a master commit because the v3.5.0 release tag predates the
# GCC 15.2 fix (a Statement::bindInt64 signature mismatch in
# include/sqliteDB.hh) and fails to link under debian trixie's toolchain.
# CGP build is enabled (libgsl-dev / libsuitesparse-dev / liblpsolve55-dev
# are installed by jamg.def's %post). MYSQL=false omits the MySQL hint-DB
# code path; auxprogs (bam2hints, filterBam, bam2wig) require libbamtools
# which is intentionally not installed, so the `auxprogs` make target is
# tolerated as a soft failure.
# -------------------------------------------------------------------------
AUGUSTUS_REPO=https://github.com/Gaius-Augustus/Augustus.git
AUGUSTUS_SHA=cd36870303ec01e8995ebcca90e48bf194ea56e3
git clone "$AUGUSTUS_REPO" Augustus
cd Augustus
git checkout "$AUGUSTUS_SHA"
make -j"$(nproc)" augustus MYSQL=false ZIPINPUT=true
make -j"$(nproc)" auxprogs MYSQL=false ZIPINPUT=true || true   # tolerate auxprog failures
mkdir -p "$JAMG/share/Augustus"
cp -a bin scripts config "$JAMG/share/Augustus/"
# Expose binaries
for b in "$JAMG/share/Augustus/bin/"*; do
    [ -x "$b" ] && ln -sf "$b" "$JAMG/bin/$(basename "$b")"
done
# AUGUSTUS_CONFIG_PATH is exported from jamg.def's %environment so it is
# present under `apptainer exec` (not just login shells via /etc/profile.d).
cd "$SRC"

# -------------------------------------------------------------------------
# EVidenceModeler 2.1.0 (annotation evidence integrator; pure Perl, no compile).
# -------------------------------------------------------------------------
EVM_REPO=https://github.com/EVidenceModeler/EVidenceModeler.git
EVM_TAG=EVidenceModeler-v2.1.0
git clone --depth 1 --branch "$EVM_TAG" "$EVM_REPO" EVidenceModeler
mkdir -p "$JAMG/share/EVidenceModeler"
cp -a EVidenceModeler/* "$JAMG/share/EVidenceModeler/"
# EVM's top-level `EVidenceModeler` script is the user entry point
ln -sf "$JAMG/share/EVidenceModeler/EVidenceModeler" "$JAMG/bin/EVidenceModeler"
cd "$SRC"

# -------------------------------------------------------------------------
# TransDecoder 5.7.1 (ORF prediction from transcripts; mostly Perl + minor C).
# -------------------------------------------------------------------------
TRANSDECODER_REPO=https://github.com/TransDecoder/TransDecoder.git
TRANSDECODER_TAG=TransDecoder-v5.7.1
git clone --depth 1 --branch "$TRANSDECODER_TAG" "$TRANSDECODER_REPO" TransDecoder
cd TransDecoder
make
mkdir -p "$JAMG/share/TransDecoder"
# TransDecoder 5.7.1 ships: TransDecoder.LongOrfs, TransDecoder.Predict, util,
# PerlLib. Older versions shipped a pfam/ dir; 5.7.1 dropped it.
cp -a TransDecoder.LongOrfs TransDecoder.Predict util PerlLib \
    "$JAMG/share/TransDecoder/"
ln -sf "$JAMG/share/TransDecoder/TransDecoder.LongOrfs" "$JAMG/bin/TransDecoder.LongOrfs"
ln -sf "$JAMG/share/TransDecoder/TransDecoder.Predict" "$JAMG/bin/TransDecoder.Predict"
cd "$SRC"

# -------------------------------------------------------------------------
# cdbfasta / cdbyank (random-access indexed FASTA; trivial Makefile).
# Source-pinned to gpertea/cdbfasta@da8f5ba5 (containers/source-pins.toml).
# -------------------------------------------------------------------------
CDBFASTA_REPO=https://github.com/gpertea/cdbfasta.git
CDBFASTA_SHA=da8f5ba5
git clone "$CDBFASTA_REPO" cdbfasta
cd cdbfasta
git checkout "$CDBFASTA_SHA"
make
install -m 0755 cdbfasta cdbyank "$JAMG/bin/"
cd "$SRC"

# -------------------------------------------------------------------------
# ParaFly (shell-command parallel dispatcher used by legacy bin/*.pl).
# Pinned to ParaFly/ParaFly@44487e0f.
#
# ParaFly's autotools configure does NOT auto-detect or auto-link OpenMP;
# CXXFLAGS / CFLAGS / LDFLAGS must be exported with -fopenmp before
# `./configure` for the generated Makefile to pick them up. Without these,
# the final link fails with "undefined reference to omp_set_dynamic /
# omp_set_num_threads / omp_get_thread_num". This is a permanent property
# of the ParaFly build, not a temporary patch.
# -------------------------------------------------------------------------
PARAFLY_REPO=https://github.com/ParaFly/ParaFly.git
PARAFLY_SHA=44487e0f
git clone "$PARAFLY_REPO" ParaFly
cd ParaFly
git checkout "$PARAFLY_SHA"
CXXFLAGS="-fopenmp -O2" CFLAGS="-fopenmp -O2" LDFLAGS="-fopenmp" \
    ./configure --prefix="$JAMG"
make -j"$(nproc)"
make install
# EVidenceModeler hardcodes its ParaFly path as
# $PLUGINS_DIR/ParaFly/bin/ParaFly (EVidenceModeler line 327).
# Symlink so EVM's pipeline finds the binary without modification.
mkdir -p "$JAMG/share/EVidenceModeler/plugins/ParaFly/bin"
ln -sf "$JAMG/bin/ParaFly" \
       "$JAMG/share/EVidenceModeler/plugins/ParaFly/bin/ParaFly"
cd "$SRC"

# -------------------------------------------------------------------------
# RepeatMasker (bioconda) -- libraries staged from host, then non-interactive
# configure.
#
# bioconda's repeatmasker package does NOT ship the repeat libraries
# (RepeatMasker.lib, RepeatPeps.lib, taxonomy.dat.bz2, etc.). The host's
# pre-curated library set is bind-mounted at /rm_lib_host (by the Makefile's
# `apptainer build --bind <RM_LIB_HOST>:/rm_lib_host:ro` call) and copied into
# $PREFIX/share/RepeatMasker/Libraries/ before the configure step.
# -------------------------------------------------------------------------
if [ ! -d "/rm_lib_host" ]; then
    echo "FATAL: /rm_lib_host not bind-mounted; cannot stage RepeatMasker libraries." >&2
    echo "Build via 'make sif' (which sets --bind), not bare apptainer build." >&2
    exit 1
fi
mkdir -p "$PREFIX/share/RepeatMasker/Libraries"
# `cp -a` preserves symlinks verbatim. If the host library set contains
# symlinks pointing OUTSIDE /rm_lib_host (e.g. to a host filer path), those
# links would dangle inside the squashed SIF. Use `-aL` to dereference if
# that becomes an issue.
cp -a /rm_lib_host/. "$PREFIX/share/RepeatMasker/Libraries/"

"$PREFIX/share/RepeatMasker/configure" \
    -trf_prgm="/usr/bin/trf" \
    -default_search_engine=rmblast \
    -rmblast_dir="$PREFIX/bin" \
    -libdir="$PREFIX/share/RepeatMasker/Libraries" </dev/null

# Expose RepeatMasker util scripts on PATH (rmOutToGFF3.pl lives in util/).
# A bare symlink into /opt/conda/bin/ is broken because Perl resolves @INC from
# the SCRIPT'S DIRECTORY (not the symlink target), so CrossmatchSearchEngine.pm
# at $PREFIX/share/RepeatMasker/ is not found. Install thin wrapper scripts
# that set PERL5LIB and exec the real script via its absolute path.
for s in "$PREFIX/share/RepeatMasker/util/"*.pl; do
    [ -e "$s" ] || continue
    name=$(basename "$s")
    # Bioconda pre-creates symlinks at $PREFIX/bin/<name> pointing back into
    # share/RepeatMasker/util/. A `cat > <path>` that follows that symlink
    # would write the wrapper INTO the real util script (O_CREAT|O_TRUNC
    # follows symlinks). Unlink first so the heredoc creates a fresh file.
    rm -f "$PREFIX/bin/$name"
    cat > "$PREFIX/bin/$name" <<WRAPEOF
#!/bin/sh
export PERL5LIB="$PREFIX/share/RepeatMasker:\${PERL5LIB-}"
exec "$s" "\$@"
WRAPEOF
    chmod +x "$PREFIX/bin/$name"
done

# Cleanup: remove the source tree (jamg.def's %post does a broader sweep too,
# but doing it here saves disk during the build).
rm -rf "$SRC"

# Make everything world-readable + executable where appropriate.
chmod -R a+rX "$JAMG"
