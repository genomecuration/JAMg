#!/usr/bin/env bash
set -euo pipefail

PREFIX=/opt/conda
SRC=/opt/src
mkdir -p "$SRC" && cd "$SRC"

# Pinned upstream commits (verified 2026-05-18 by Explore agent):
EXONERATE_REPO=https://github.com/nathanweeks/exonerate.git
EXONERATE_SHA=22beeb02
CDBFASTA_REPO=https://github.com/gpertea/cdbfasta.git
CDBFASTA_SHA=da8f5ba5
PARAFLY_REPO=https://github.com/ParaFly/ParaFly.git
PARAFLY_SHA=44487e0f

# 1) exonerate (maintained fork; bioconda is 2017)
git clone "$EXONERATE_REPO" exonerate && cd exonerate
git checkout "$EXONERATE_SHA"
autoreconf -fi
./configure --prefix="$PREFIX"
make -j"$(nproc)"
make install
cd "$SRC"

# 2) cdbtools (cdbfasta/cdbyank; bioconda is 2016)
git clone "$CDBFASTA_REPO" cdbfasta && cd cdbfasta
git checkout "$CDBFASTA_SHA"
make
install -m 0755 cdbfasta cdbyank "$PREFIX/bin/"
cd "$SRC"

# 3) ParaFly (bioconda is 2013; legacy tools in bin/ still call it)
git clone "$PARAFLY_REPO" ParaFly && cd ParaFly
git checkout "$PARAFLY_SHA"
./configure --prefix="$PREFIX"
make -j"$(nproc)"
make install
cd "$SRC"

# Cleanup source trees to keep SIF small
rm -rf "$SRC"

# 4) Configure RepeatMasker non-interactively (review-#2 I16).
#    RepeatMasker 4.2.x accepts these args directly; no answer-stream parsing.
"$PREFIX/share/RepeatMasker/configure" \
  -trf_prgm="$PREFIX/bin/trf" \
  -default_search_engine=rmblast \
  -rmblast_dir="$PREFIX/bin" \
  -libdir="$PREFIX/share/RepeatMasker/Libraries" </dev/null

# 5) Expose RepeatMasker util scripts on PATH (rmOutToGFF3.pl lives in util/)
for s in "$PREFIX/share/RepeatMasker/util/"*.pl; do
  ln -sf "$s" "$PREFIX/bin/$(basename "$s")"
done

chmod -R a+rX /opt/jamg
