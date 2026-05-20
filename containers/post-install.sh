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
# EVidenceModeler (annotation evidence integrator; pure Perl, no compile).
# Sourced from the alpapan/EVidenceModeler submodule, vendored into the SIF
# build context via jamg.def %files at /opt/jamg-src/EVidenceModeler. The
# `cp -a .../* …` glob skips the submodule's .git/ + dotfiles; the staging
# directory is removed at the end of this script so it does not bloat the
# final SIF.
# -------------------------------------------------------------------------
mkdir -p "$JAMG/share/EVidenceModeler"
cp -a /opt/jamg-src/EVidenceModeler/* "$JAMG/share/EVidenceModeler/"
ln -sf "$JAMG/share/EVidenceModeler/EVidenceModeler" "$JAMG/bin/EVidenceModeler"
rm -rf /opt/jamg-src/EVidenceModeler
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

# python3 entry-point shim. The SIF's %environment puts /usr/bin BEFORE
# /opt/conda/bin so that `env perl` resolves to the apt-installed system
# perl 5.40 with BioPerl + DB_File (a deliberate single-source-of-truth
# decision in jamg-base.def). For python3 the same ordering hits a wall:
# RepeatMasker's bioconda package ships scripts with `#!/usr/bin/env
# python3` (famdb.py, util/RM2Bed.py) that need h5py, and /usr/bin/python3
# does not have h5py while /opt/conda/bin/python3 does. Shimming a
# python3 link under /opt/jamg/bin (which precedes both) makes env find
# the conda python first without changing perl's resolution.
ln -sf /opt/conda/bin/python3 "$JAMG/bin/python3"

cd "$SRC"

# -------------------------------------------------------------------------
# RepeatMasker (bioconda) -- libraries staged from the SIF build context,
# then non-interactive configure.
#
# bioconda's repeatmasker package does NOT ship the repeat libraries
# (RepeatMasker.lib, RepeatPeps.lib, taxonomy.dat.bz2, etc.). The library
# set is vendored under containers/rm_libs/ (LFS-tracked) and staged into
# the SIF build at /rm_lib_host by jamg.def %files, then copied into
# $PREFIX/share/RepeatMasker/Libraries/ here before the configure step.
# The /rm_lib_host staging dir is removed at the end of jamg.def %post so
# the squashed SIF does not ship 667 MB of duplicates.
# -------------------------------------------------------------------------
if [ ! -d "/rm_lib_host" ]; then
    echo "FATAL: /rm_lib_host not populated; cannot stage RepeatMasker libraries." >&2
    echo "containers/rm_libs/ is likely empty. Run 'git lfs pull' to fetch the vendored libraries, then rebuild." >&2
    exit 1
fi
mkdir -p "$PREFIX/share/RepeatMasker/Libraries"
# `cp -a` preserves symlinks verbatim. Symlinks pointing OUTSIDE the
# vendored set would dangle inside the squashed SIF; if that becomes an
# issue, use `-aL` to dereference.
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
