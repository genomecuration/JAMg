#!/usr/bin/env bash
# Reproducibly rebuild test_suite/mini-* fixtures from the bundled sources.
#
# Two fixture sets are produced:
#   - **100 kb default** (mini-genome.fasta + siblings): used by every
#     per-rule test EXCEPT genemark, and by the e2e smoke after the
#     genemark output is pre-staged from the 1 Mb fixture.
#   - **1 Mb genemark fixture** (mini-genome-1mb.fasta + siblings):
#     used ONLY by test_genemark.sh. GeneMark-ES self-training needs
#     >=1 Mb to converge, so genemark cannot run on the 100 kb default.
#     Tests downstream of genemark (test_evm, test_ogs, test_full_dag)
#     pre-stage a SUBSET of a 1 Mb genemark run (test_suite/fixtures/
#     genemark.gff3) and let snakemake skip the genemark rule via mtime.
#
# Sources used:
#   - test_suite/dmel-X-r5.53.fasta.bz2                                (mini-genome region)
#   - test_suite/Drosophila_official_annotations_cleaned.tar           (annotation -> mRNAs)
#   - databases/repeats/rnammer-SILVA.classified.nr95.renamed.fasta.bz2 (RNA repeats subset)
#   - https://ftp.uniprot.org/.../uniprot_sprot.fasta.gz               (cached at tools/cache/)
#
# Reproducibility notes:
#   - wgsim -S 42                                seed pinned
#   - python random.seed(42)                     seed pinned for species-repeats sampling
#   - SwissProt: fetched from UniProt's current-release URL, NOT a named-release
#     pin. On first run the SHA256 is recorded into tools/cache/.sha256; later
#     runs verify against that hash so subsequent builds are byte-identical to
#     the original build on this machine. Cross-machine reproducibility
#     requires copying both tools/cache/uniprot_sprot.fasta.gz AND .sha256.
#
#   - mini-repeats-species*.fasta is sampled FROM the same-size mini-genome
#     itself (no external Diptera RepeatModeler library is bundled).
#     RepeatMasker will match those windows at 100%, so up to ~5-10 kb of
#     the mini-genome will be masked when these are used. Phase-3 rule
#     tests assert on rule completion + output presence, NOT on prediction
#     quality metrics, so this masking aggressiveness is acceptable.
#
# Outputs:
#   test_suite/mini-genome.fasta      + .fai      (100 kb, X:100001-200000)
#   test_suite/mini-rnaseq.bam        + .bai      (~10k reads, ~20x coverage)
#   test_suite/mini-transcripts.fasta             (annotation in 100 kb window)
#   test_suite/mini-proteins.fasta    + BLAST DB  (SwissProt hits to 100 kb)
#   test_suite/mini-repeats-rna.fasta             (20 rnammer entries; size-independent)
#   test_suite/mini-repeats-species.fasta         (20 windows sampled from 100 kb)
#   test_suite/mini-genome-1mb.fasta   + .fai     (1 Mb, X:100001-1100000)
#   test_suite/mini-rnaseq-1mb.bam     + .bai     (~100k reads, ~20x coverage)
#   test_suite/mini-repeats-species-1mb.fasta     (20 windows sampled from 1 Mb)

set -euo pipefail

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
cd "$REPO"

SIF="$REPO/containers/jamg.sif"
if [[ ! -f "$SIF" ]]; then
    echo "FATAL: $SIF not found. Run 'make sif' first." >&2
    exit 1
fi
RUN=(apptainer exec --bind "$REPO":"$REPO" "$SIF")

W=test_suite/.work
mkdir -p "$W" tools/cache

# Decompress the source X-chromosome once; reused for both extracts.
SRC_FA=$(mktemp --suffix .fasta)
bzcat test_suite/dmel-X-r5.53.fasta.bz2 > "$SRC_FA"
"${RUN[@]}" samtools faidx "$SRC_FA" 2>/dev/null

# Untar + decompress the source annotation once; reused for both extracts.
tar -xf test_suite/Drosophila_official_annotations_cleaned.tar -C "$W" \
    melanogaster/dmel-all-no-analysis-r5.53.gff3.gff3.clean.bz2
bunzip2 -kf "$W/melanogaster/dmel-all-no-analysis-r5.53.gff3.gff3.clean.bz2"
SRC_GFF=$W/melanogaster/dmel-all-no-analysis-r5.53.gff3.gff3.clean

# ----------------------------------------------------------------------------
# build_genome_fixture <end_coord> <wgsim_N> <star_saindexnbases> \
#                      <genome_out> <bam_out> <species_repeats_out> <tag>
#
# Extracts X:100001-<end_coord> from $SRC_FA (renamed to >X_mini), simulates
# wgsim_N PE reads at 100 bp, aligns with STAR (--genomeSAindexNbases =
# star_saindexnbases), samples 20 species repeat windows from the extracted
# genome. wgsim_N should be chosen to give ~20x coverage on the extract.
# <tag> is a short label used in log filenames ("100kb", "1mb").
# ----------------------------------------------------------------------------
build_genome_fixture() {
    local end="$1"
    local wgsim_n="$2"
    local star_sai="$3"
    local genome_out="$4"
    local bam_out="$5"
    local species_out="$6"
    local tag="$7"          # short label for log lines

    # --- genome ---
    "${RUN[@]}" samtools faidx "$SRC_FA" "X:100001-$end" 2>/dev/null \
        | sed 's|^>X:.*|>X_mini|' > "$genome_out"
    "${RUN[@]}" samtools faidx "$genome_out" 2>/dev/null

    # --- wgsim PE reads + STAR alignment ---
    local r1="$W/${tag}_r1.fq" r2="$W/${tag}_r2.fq"
    rm -f "$r1" "$r2"
    "${RUN[@]}" wgsim -S 42 -N "$wgsim_n" -1 100 -2 100 -e 0.005 -r 0.001 -R 0 -X 0 \
        "$genome_out" "$r1" "$r2" \
        > "$W/${tag}_wgsim.mut" 2>"$W/${tag}_wgsim.err"

    local star_idx="$W/${tag}_star_idx"
    rm -rf "$star_idx"; mkdir -p "$star_idx"
    "${RUN[@]}" STAR --runMode genomeGenerate \
        --genomeDir "$star_idx" \
        --genomeFastaFiles "$genome_out" \
        --genomeSAindexNbases "$star_sai" --runThreadN 4 \
        >"$W/${tag}_star_index.log" 2>&1

    "${RUN[@]}" STAR --runMode alignReads --runThreadN 4 \
        --genomeDir "$star_idx" \
        --readFilesIn "$r1" "$r2" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "$W/${tag}_star_aln_" \
        --alignIntronMax 50000 \
        >"$W/${tag}_star_align.log" 2>&1

    cp "$W/${tag}_star_aln_Aligned.sortedByCoord.out.bam" "$bam_out"
    "${RUN[@]}" samtools index "$bam_out" 2>/dev/null

    # Sanity: at least 10% of reads mapped (lower bar than 1 Mb because
    # short fixtures show more edge-effects).
    local mapped
    mapped=$("${RUN[@]}" samtools flagstat "$bam_out" 2>/dev/null \
        | awk '/^[0-9]+ \+ [0-9]+ primary mapped/{print $1; exit}')
    local min=$(( wgsim_n * 2 / 10 ))   # 10% of 2*N (PE = 2 reads/pair)
    if ! [[ "$mapped" =~ ^[0-9]+$ ]] || (( mapped < min )); then
        echo "FATAL: only $mapped mapped reads in $bam_out (need >=$min for ~20x cov)" >&2
        exit 1
    fi

    # --- species-repeats: 20 sampled windows from this genome ---
    python3 - "$genome_out" "$species_out" <<'PY'
import random, sys
random.seed(42)
genome_path, out_path = sys.argv[1], sys.argv[2]
seq = ""
with open(genome_path) as fh:
    for line in fh:
        if line.startswith(">"): continue
        seq += line.strip()
N = len(seq)
with open(out_path, "w") as out:
    for i in range(20):
        L = random.randint(150, 500)
        start = random.randint(0, N - L)
        chunk = seq[start:start+L]
        out.write(f">species_repeat_{i+1}#LINE/L1\n")
        for j in range(0, len(chunk), 80):
            out.write(chunk[j:j+80] + "\n")
PY
}

# ----------------------------------------------------------------------------
# 100 kb default fixture (X:100001-200000)
# 1k pairs * 2 reads/pair * 100 bp = 200 kb sequencing on 100 kb ≈ 2x coverage.
# 2x is enough for STAR to align + Trinity-GG to assemble a handful of
# transcripts; we are NOT testing assembly quality, just rule wiring.
# Trinity-GG runtime scales with read count, not genome size, so the 10x
# read-count cut from the original 10k is the largest single test-time win.
# SAindexNbases 7 = min(14, log2(100_000)/2 - 1) ≈ 7.3 -> 7.
# ----------------------------------------------------------------------------
build_genome_fixture \
    200000 1000 7 \
    test_suite/mini-genome.fasta \
    test_suite/mini-rnaseq.bam \
    test_suite/mini-repeats-species.fasta \
    100kb

# ----------------------------------------------------------------------------
# 1 Mb genemark fixture (X:100001-1100000)
# 100k pairs * 2 reads/pair * 100 bp = 20 Mb sequencing on 1 Mb ≈ 20x coverage.
# SAindexNbases 9 = min(14, log2(1_000_000)/2 - 1) ≈ 8.96 -> 9 (10 overshoots
# the formula's cap and can segfault STAR's index builder on small genomes).
# Used only by tests/rules/test_genemark.sh; gmes_petap --ES needs >=1 Mb.
# ----------------------------------------------------------------------------
build_genome_fixture \
    1100000 100000 9 \
    test_suite/mini-genome-1mb.fasta \
    test_suite/mini-rnaseq-1mb.bam \
    test_suite/mini-repeats-species-1mb.fasta \
    1mb

# ----------------------------------------------------------------------------
# Annotation-derived fixtures (transcripts + proteins). Default is 100 kb;
# downstream tests on the 100 kb genome use these.
# ----------------------------------------------------------------------------
awk -F'\t' -v OFS='\t' -v off=100000 '
    /^#/ { print; next }
    $1=="X" && $4 >= 100001 && $5 <= 200000 && ($3=="gene"||$3=="mRNA"||$3=="exon"||$3=="CDS"||$3=="five_prime_UTR"||$3=="three_prime_UTR") {
        $1 = "X_mini"; $4 -= off; $5 -= off; print
    }
' "$SRC_GFF" > "$W/mini-annot.gff3"

"${RUN[@]}" gffread -w test_suite/mini-transcripts.fasta \
    -g test_suite/mini-genome.fasta "$W/mini-annot.gff3" 2>/dev/null

# --- mini-proteins.fasta (subset of SwissProt with hits to mini-genome) ----
SPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
SPROT_GZ=tools/cache/uniprot_sprot.fasta.gz
SPROT_FA=tools/cache/uniprot_sprot.fasta

if [[ ! -s "$SPROT_GZ" ]]; then
    curl -fsSL -o "$SPROT_GZ.partial" "$SPROT_URL"
    mv "$SPROT_GZ.partial" "$SPROT_GZ"
fi
if [[ -s tools/cache/.sha256 ]]; then
    sha256sum -c tools/cache/.sha256 \
        || { echo "FATAL: SwissProt SHA mismatch against tools/cache/.sha256" >&2; exit 1; }
else
    sha256sum "$SPROT_GZ" > tools/cache/.sha256
fi
[[ -s "$SPROT_FA" ]] || gunzip -k "$SPROT_GZ"

if [[ ! -s "$SPROT_FA.pdb" ]]; then
    "${RUN[@]}" makeblastdb -in "$SPROT_FA" -dbtype prot -parse_seqids -hash_index -title SwissProt >/dev/null
fi

"${RUN[@]}" blastx -query test_suite/mini-genome.fasta -db "$SPROT_FA" \
    -outfmt 6 -max_target_seqs 500 -num_threads "${SLURM_CPUS_PER_TASK:-20}" -evalue 1e-10 \
    -out "$W/mini_blastx.tsv" 2>/dev/null

cut -f2 "$W/mini_blastx.tsv" | sort -u > "$W/prot_ids.txt"
"${RUN[@]}" blastdbcmd -db "$SPROT_FA" -entry_batch "$W/prot_ids.txt" \
    -out test_suite/mini-proteins.fasta 2>/dev/null

if [[ ! -s test_suite/mini-proteins.fasta.pdb ]]; then
    "${RUN[@]}" makeblastdb -in test_suite/mini-proteins.fasta -dbtype prot \
        -parse_seqids -hash_index -title "mini SwissProt subset" >/dev/null
fi

# --- mini-repeats-rna.fasta (20 entries from rnammer-SILVA) ---------------
# Genome-size-independent; built once.
bzcat databases/repeats/rnammer-SILVA.classified.nr95.renamed.fasta.bz2 > "$W/rnammer.fa"
awk 'BEGIN{n=0} /^>/{n++; if(n>20) exit} {print}' "$W/rnammer.fa" > test_suite/mini-repeats-rna.fasta
rm -f "$W/rnammer.fa"

rm -f "$SRC_FA" "${SRC_FA}.fai"

echo "OK: fixtures rebuilt (100 kb default + 1 Mb genemark)"
ls -la test_suite/mini-*
