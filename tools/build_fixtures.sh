#!/usr/bin/env bash
# Reproducibly rebuild test_suite/mini-* fixtures from the bundled sources.
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
#     If you need a true upstream-release pin, replace SPROT_URL with a
#     previous_releases/release-YYYY_MM/... path that resolves.
#
#   - mini-repeats-species.fasta is sampled FROM the mini-genome itself
#     (no external Diptera RepeatModeler library is bundled). RepeatMasker
#     will match those windows at 100%, so up to ~5-10 kb of the mini-genome
#     will be masked when these are used. Phase-3 rule tests assert on rule
#     completion + output presence, NOT on prediction quality metrics, so
#     this masking aggressiveness is acceptable for fixture purposes.
#
# Outputs:
#   test_suite/mini-genome.fasta + .fai
#   test_suite/mini-rnaseq.bam   + .bai
#   test_suite/mini-transcripts.fasta
#   test_suite/mini-proteins.fasta + BLAST DB siblings
#   test_suite/mini-repeats-rna.fasta
#   test_suite/mini-repeats-species.fasta

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

# --- mini-genome ---------------------------------------------------------
# 1 Mb region of D. melanogaster X chromosome, well past the 5 kb telomeric
# zone. Contains the original 51 kb gene-rich window plus surrounding
# context for tools that require minimum genome sizes (e.g. GeneMark-ES
# self-training expects >=1 Mb).
TMP=$(mktemp --suffix .fasta)
bzcat test_suite/dmel-X-r5.53.fasta.bz2 > "$TMP"
"${RUN[@]}" samtools faidx "$TMP" 2>/dev/null
"${RUN[@]}" samtools faidx "$TMP" X:100001-1100000 2>/dev/null \
    | sed 's|^>X:100001-1100000|>X_mini|' > test_suite/mini-genome.fasta
rm -f "$TMP" "${TMP}.fai"
"${RUN[@]}" samtools faidx test_suite/mini-genome.fasta 2>/dev/null

# --- mini-rnaseq.bam -----------------------------------------------------
# 100k PE reads over 1 Mb is ~20x coverage at 100bp (sufficient for STAR
# alignment + augustus coverage hints).
rm -f "$W/mini_r1.fq" "$W/mini_r2.fq"
"${RUN[@]}" wgsim -S 42 -N 100000 -1 100 -2 100 -e 0.005 -r 0.001 -R 0 -X 0 \
    test_suite/mini-genome.fasta \
    "$W/mini_r1.fq" "$W/mini_r2.fq" > "$W/wgsim.mut" 2>"$W/wgsim.err"

rm -rf "$W/star_idx"
mkdir -p "$W/star_idx"
# SAindexNbases 10 for ~1 Mb (rule of thumb: min(14, log2(genome)/2 - 1)).
"${RUN[@]}" STAR --runMode genomeGenerate \
    --genomeDir "$W/star_idx" \
    --genomeFastaFiles test_suite/mini-genome.fasta \
    --genomeSAindexNbases 10 --runThreadN 4 \
    >"$W/star_index.log" 2>&1

"${RUN[@]}" STAR --runMode alignReads --runThreadN 4 \
    --genomeDir "$W/star_idx" \
    --readFilesIn "$W/mini_r1.fq" "$W/mini_r2.fq" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$W/star_aln_" \
    --alignIntronMax 50000 \
    >"$W/star_align.log" 2>&1

cp "$W/star_aln_Aligned.sortedByCoord.out.bam" test_suite/mini-rnaseq.bam
"${RUN[@]}" samtools index test_suite/mini-rnaseq.bam 2>/dev/null

# samtools flagstat lines/format has been stable since 1.3; first awk match
# on /primary mapped/ matches the "N + 0 primary mapped (P% : N/A)" row.
mapped=$("${RUN[@]}" samtools flagstat test_suite/mini-rnaseq.bam 2>/dev/null \
    | awk '/^[0-9]+ \+ [0-9]+ primary mapped/{print $1; exit}')
if ! [[ "$mapped" =~ ^[0-9]+$ ]] || (( mapped < 10000 )); then
    echo "FATAL: only $mapped mapped reads in mini-rnaseq.bam (need >=10000)" >&2
    exit 1
fi

# --- mini-transcripts.fasta + mini-annot.gff3 (mini-genome-local coords) ---
tar -xf test_suite/Drosophila_official_annotations_cleaned.tar -C "$W" \
    melanogaster/dmel-all-no-analysis-r5.53.gff3.gff3.clean.bz2
bunzip2 -kf "$W/melanogaster/dmel-all-no-analysis-r5.53.gff3.gff3.clean.bz2"
SRC_GFF=$W/melanogaster/dmel-all-no-analysis-r5.53.gff3.gff3.clean

awk -F'\t' -v OFS='\t' -v off=100000 '
    /^#/ { print; next }
    $1=="X" && $4 >= 100001 && $5 <= 1100000 && ($3=="gene"||$3=="mRNA"||$3=="exon"||$3=="CDS"||$3=="five_prime_UTR"||$3=="three_prime_UTR") {
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
    -outfmt 6 -max_target_seqs 500 -num_threads 4 -evalue 1e-10 \
    -out "$W/mini_blastx.tsv" 2>/dev/null

cut -f2 "$W/mini_blastx.tsv" | sort -u > "$W/prot_ids.txt"
"${RUN[@]}" blastdbcmd -db "$SPROT_FA" -entry_batch "$W/prot_ids.txt" \
    -out test_suite/mini-proteins.fasta 2>/dev/null

if [[ ! -s test_suite/mini-proteins.fasta.pdb ]]; then
    "${RUN[@]}" makeblastdb -in test_suite/mini-proteins.fasta -dbtype prot \
        -parse_seqids -hash_index -title "mini SwissProt subset" >/dev/null
fi

# --- mini-repeats-rna.fasta (20 entries from rnammer-SILVA) ---------------
bzcat databases/repeats/rnammer-SILVA.classified.nr95.renamed.fasta.bz2 \
    | awk 'BEGIN{n=0} /^>/{n++; if(n>20) exit} {print}' > test_suite/mini-repeats-rna.fasta

# --- mini-repeats-species.fasta (20 sampled windows from mini-genome) -----
python3 - <<'PY' > test_suite/mini-repeats-species.fasta
import random
random.seed(42)
seq = ""
with open("test_suite/mini-genome.fasta") as fh:
    for line in fh:
        if line.startswith(">"): continue
        seq += line.strip()
N = len(seq)
for i in range(20):
    L = random.randint(150, 500)
    start = random.randint(0, N - L)
    chunk = seq[start:start+L]
    print(f">species_repeat_{i+1}#LINE/L1")
    for j in range(0, len(chunk), 80):
        print(chunk[j:j+80])
PY

echo "OK: fixtures rebuilt"
ls -la test_suite/mini-*
