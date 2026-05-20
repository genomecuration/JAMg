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
#   - polyester seed=42                          PE RNA-seq simulation seed pinned
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
#   test_suite/mini-genome.fasta            + .fai      (100 kb, X:100001-200000)
#   test_suite/mini-transcripts.fasta                   (spliced transcripts in 100 kb window)
#   test_suite/mini-rnaseq.bam              + .bai      (polyester PE reads, 10x per-transcript)
#   test_suite/mini-proteins.fasta          + BLAST DB  (SwissProt hits to 100 kb)
#   test_suite/mini-repeats-rna.fasta                   (20 rnammer entries; size-independent)
#   test_suite/mini-repeats-species.fasta               (20 windows sampled from 100 kb)
#   test_suite/mini-genome-1mb.fasta        + .fai      (1 Mb, X:100001-1100000)
#   test_suite/mini-transcripts-1mb.fasta               (spliced transcripts in 1 Mb window)
#   test_suite/mini-rnaseq-1mb.bam          + .bai      (polyester PE reads, 10x per-transcript)
#   test_suite/mini-repeats-species-1mb.fasta           (20 windows sampled from 1 Mb)

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
# build_genome_fixture <start_coord> <end_coord> <target_cov> <star_saindexnbases> \
#                      <genome_out> <transcripts_out> <bam_out> <species_repeats_out> <tag>
#
# Extracts X:<start_coord>-<end_coord> from $SRC_FA (renamed to >X_mini),
# subsets the source GFF to that window, extracts spliced transcripts via
# gffread, then drives Bioconductor polyester to simulate paired-end RNA-seq
# reads FROM the transcripts FASTA at <target_cov> per-transcript coverage.
# STAR aligns those reads to the genome with splice-aware mode so the BAM
# carries real intron-spanning split reads (which augustus RNA-seq hints and
# genemark --ET both need). Also samples 20 species repeat windows from the
# extracted genome. <tag> is a short label used in log filenames.
# ----------------------------------------------------------------------------
build_genome_fixture() {
    local start="$1"
    local end="$2"
    local target_cov="$3"
    local star_sai="$4"
    local genome_out="$5"
    local transcripts_out="$6"
    local bam_out="$7"
    local species_out="$8"
    local tag="$9"          # short label for log lines

    local off=$(( start - 1 ))      # offset for GFF coord rebasing

    # --- genome ---
    "${RUN[@]}" samtools faidx "$SRC_FA" "X:$start-$end" 2>/dev/null \
        | sed 's|^>X:.*|>X_mini|' > "$genome_out"
    "${RUN[@]}" samtools faidx "$genome_out" 2>/dev/null

    # --- annotation window -> spliced transcripts FASTA ---
    awk -F'\t' -v OFS='\t' -v start="$start" -v end="$end" -v off="$off" '
        /^#/ { print; next }
        $1=="X" && $4 >= start && $5 <= end && ($3=="gene"||$3=="mRNA"||$3=="exon"||$3=="CDS"||$3=="five_prime_UTR"||$3=="three_prime_UTR") {
            $1 = "X_mini"; $4 -= off; $5 -= off; print
        }
    ' "$SRC_GFF" > "$W/${tag}_annot.gff3"

    "${RUN[@]}" gffread -w "$transcripts_out" \
        -g "$genome_out" "$W/${tag}_annot.gff3" 2>/dev/null

    # --- polyester paired-end RNA-seq simulation from transcripts ---
    # Per-transcript coverage = (reads_per_tx * 2 * 100 bp) / tx_length.
    # Solve for reads_per_tx: ceil(tx_length * target_cov / 200). Floor at 50
    # so very short transcripts still yield a paired-fragment minimum.
    local pol_out="$W/${tag}_polyester"
    rm -rf "$pol_out"; mkdir -p "$pol_out"
    "${RUN[@]}" Rscript - "$transcripts_out" "$pol_out" "$target_cov" \
        > "$W/${tag}_polyester.log" 2>&1 <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
fasta_path <- args[1]
outdir <- args[2]
target_cov <- as.numeric(args[3])
suppressPackageStartupMessages({ library(polyester); library(Biostrings) })
txs <- readDNAStringSet(fasta_path)
reads_per_tx <- pmax(50L, as.integer(ceiling(width(txs) * target_cov / 200)))
countmat <- matrix(reads_per_tx, ncol = 1L)
simulate_experiment_countmat(
    fasta_path, readmat = countmat, outdir = outdir,
    paired = TRUE, readlen = 100L,
    fraglen = 200, fragsd = 25,
    seed = 42L, strand_specific = TRUE,
    error_model = 'illumina5', bias = 'none'
)
cat(sprintf('polyester wrote %d transcripts, total reads = %d\n',
    length(txs), sum(reads_per_tx)))
RSCRIPT

    local r1="$pol_out/sample_01_1.fasta"
    local r2="$pol_out/sample_01_2.fasta"
    [[ -s "$r1" && -s "$r2" ]] || {
        echo "FATAL: polyester did not produce paired FASTAs in $pol_out" >&2
        tail -50 "$W/${tag}_polyester.log" >&2
        exit 1
    }

    # --- STAR alignment to genome (splice-aware) ---
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

    # Sanity 1: at least 50% of reads must map (transcripts are subset of the
    # genome, so mapping rate should be near 100% absent edge effects).
    local mapped total
    mapped=$("${RUN[@]}" samtools flagstat "$bam_out" 2>/dev/null \
        | awk '/^[0-9]+ \+ [0-9]+ primary mapped/{print $1; exit}')
    total=$("${RUN[@]}" samtools flagstat "$bam_out" 2>/dev/null \
        | awk '/^[0-9]+ \+ [0-9]+ primary$/{print $1; exit}')
    if ! [[ "$mapped" =~ ^[0-9]+$ && "$total" =~ ^[0-9]+$ ]] || (( total == 0 )) || (( mapped * 2 < total )); then
        echo "FATAL: only $mapped/$total reads mapped in $bam_out (need >=50%)" >&2
        exit 1
    fi

    # Sanity 2: BAM must carry intron-spanning split reads (CIGAR contains
    # an N operator). Zero N-cigar reads means the fixture is no better than
    # the old wgsim-from-genome and the genemark --ET / augustus hint paths
    # will not exercise their intron logic.
    local n_split
    n_split=$("${RUN[@]}" samtools view "$bam_out" 2>/dev/null \
        | awk '$6 ~ /N/ {n++} END {print n+0}')
    if (( n_split < 1 )); then
        echo "FATAL: $bam_out has zero N-CIGAR (intron-spanning) reads" >&2
        exit 1
    fi
    echo "OK: $tag fixture: $mapped/$total mapped, $n_split intron-spanning reads"

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
# Polyester samples 10x per-transcript coverage from the spliced transcripts.
# SAindexNbases 7 = min(14, log2(100_000)/2 - 1) ≈ 7.3 -> 7.
# ----------------------------------------------------------------------------
build_genome_fixture \
    100001 200000 10 7 \
    test_suite/mini-genome.fasta \
    test_suite/mini-transcripts.fasta \
    test_suite/mini-rnaseq.bam \
    test_suite/mini-repeats-species.fasta \
    100kb

# ----------------------------------------------------------------------------
# 1 Mb genemark fixture (X:100001-1100000)
# 10x per-transcript coverage (same as 100 kb; coverage is per-transcript, not
# per-genome-bp, so the same target gives roughly the same read budget per gene).
# SAindexNbases 9 = min(14, log2(1_000_000)/2 - 1) ≈ 8.96 -> 9 (10 overshoots
# the formula's cap and can segfault STAR's index builder on small genomes).
# Used by tests/rules/test_genemark.sh; gmes_petap needs >=1 Mb.
# ----------------------------------------------------------------------------
build_genome_fixture \
    100001 1100000 10 9 \
    test_suite/mini-genome-1mb.fasta \
    test_suite/mini-transcripts-1mb.fasta \
    test_suite/mini-rnaseq-1mb.bam \
    test_suite/mini-repeats-species-1mb.fasta \
    1mb

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

# Pin source-fixture mtimes to a fixed past date so every derived artifact
# (snapshot files, staged outputs, fresh rule outputs) is naturally newer.
# The glob also touches a handful of runtime-derived files left at the
# test_suite/ root (mini-genome.fasta.cidx, mini-transcripts.fasta.cidx/orf/
# annotations*). Those files are regenerated by tools that overwrite mtime,
# and snakemake does not track them as named inputs/outputs; touching them is harmless.
find test_suite -maxdepth 1 -name 'mini-*' -type f -exec touch -h -d '2020-01-01' {} +

echo "OK: fixtures rebuilt (100 kb default + 1 Mb genemark)"
ls -la test_suite/mini-*
