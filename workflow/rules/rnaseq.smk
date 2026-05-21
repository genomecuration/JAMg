# RNA-seq evidence pipeline. Three input modes selected by
# config["rnaseq"]["mode"]:
#   - "bam":   user-supplied coord-sorted BAM (passthrough; resort+reindex)
#   - "star":  align rnaseq.fastq_pairs[] with STAR
#   - "gsnap": align rnaseq.fastq_pairs[] with gsnap | samtools sort
# All three converge on {OUTDIR}/rnaseq/rnaseq.sorted.bam + augustus hints
# produced by augustus_RNAseq_hints.pl.

import os as _os

_GENOME_ABS = _os.path.abspath(config["genome"])
_RNASEQ_BAM_ABS = (
    _os.path.abspath(config["rnaseq"]["bam"])
    if config["rnaseq"].get("bam")
    else None
)
_RNASEQ_MODE = config["rnaseq"]["mode"]


def _aligned_bam_input(_wc):
    """Resolve the unsorted (or upstream-sorted) BAM for the current mode."""
    if _RNASEQ_MODE == "bam":
        return _RNASEQ_BAM_ABS
    if _RNASEQ_MODE == "star":
        return f"{OUTDIR}/rnaseq/star/Aligned.sortedByCoord.out.bam"
    if _RNASEQ_MODE == "gsnap":
        return f"{OUTDIR}/rnaseq/gsnap/aligned.sorted.bam"
    raise ValueError(f"Unknown rnaseq.mode: {_RNASEQ_MODE}")


rule rnaseq_align_star:
    # Built only when mode=star. STAR aligns rnaseq.fastq_pairs[] and writes
    # a coord-sorted BAM at the canonical Aligned.sortedByCoord.out.bam path.
    input:
        genome = config["genome"],
    output:
        bam = f"{OUTDIR}/rnaseq/star/Aligned.sortedByCoord.out.bam",
    params:
        genome_abs = _GENOME_ABS,
        max_intron = config["max_intron"],
        reads = lambda _wc: [
            fq for pair in (config["rnaseq"].get("fastq_pairs") or []) for fq in pair
        ],
        idx_dir = f"{OUTDIR}/rnaseq/star/idx",
        out_prefix = f"{OUTDIR}/rnaseq/star/Aligned_",
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 32000,
    shell:
        "mkdir -p {params.idx_dir} && "
        "STAR --runMode genomeGenerate --genomeDir {params.idx_dir} "
        "  --genomeFastaFiles {params.genome_abs} "
        "  --genomeSAindexNbases 8 --runThreadN {threads} && "
        "STAR --runMode alignReads --runThreadN {threads} "
        "  --genomeDir {params.idx_dir} "
        "  --readFilesIn {params.reads} "
        "  --outSAMtype BAM SortedByCoordinate "
        "  --outFileNamePrefix {params.out_prefix} "
        "  --alignIntronMax {params.max_intron} && "
        "mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam}"


rule rnaseq_align_gsnap:
    # Built only when mode=gsnap. gsnap streams SAM to samtools sort.
    input:
        genome = config["genome"],
    output:
        bam = f"{OUTDIR}/rnaseq/gsnap/aligned.sorted.bam",
    params:
        genome_abs = _GENOME_ABS,
        reads = lambda _wc: [
            fq for pair in (config["rnaseq"].get("fastq_pairs") or []) for fq in pair
        ],
        gmap_dir = f"{OUTDIR}/rnaseq/gsnap/gmapdb",
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "mkdir -p {params.gmap_dir} && "
        "gmap_build -D {params.gmap_dir} -d mini_genome {params.genome_abs} && "
        "gsnap -D {params.gmap_dir} -d mini_genome -t {threads} "
        "  --novelsplicing=1 --format=sam {params.reads} "
        "  | samtools sort -@ {threads} -o {output.bam} -"


rule rnaseq_sort_bam:
    # Mode-agnostic: produces the canonical sorted+indexed BAM that
    # downstream rules consume.
    input:
        bam = _aligned_bam_input,
    output:
        bam = f"{OUTDIR}/rnaseq/rnaseq.sorted.bam",
        bai = f"{OUTDIR}/rnaseq/rnaseq.sorted.bam.bai",
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 8000,
    shell:
        "mkdir -p $(dirname {output.bam}) && "
        # samtools sort works whether the input is sorted or not. If the
        # input is already coord-sorted (mode=bam, mode=star), sort is a
        # no-op pass-through write.
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && "
        "samtools index {output.bam}"


rule rnaseq_hints:
    # augustus_RNAseq_hints.pl writes coverage hints to <bam>.coverage.hints
    # and junction hints to <bam>.junctions.hints (when junctions are found).
    # Concatenate both into rnaseq.hints. With purely DNA-simulated reads
    # (no spliced reads) the junctions output may be empty; coverage hints
    # are always produced when there is mapped depth.
    input:
        bam = rules.rnaseq_sort_bam.output.bam,
        bai = rules.rnaseq_sort_bam.output.bai,
        genome = config["genome"],
    output:
        hints = f"{OUTDIR}/rnaseq/rnaseq.hints",
    params:
        genome_abs = _GENOME_ABS,
        bam_abs    = lambda _wc, input: _os.path.abspath(input.bam),
        hints_abs  = lambda _wc, output: _os.path.abspath(output.hints),
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 8000,
    shell:
        "cd $(dirname {output.hints}) && "
        "augustus_RNAseq_hints.pl -bam {params.bam_abs} -genome {params.genome_abs} "
        "  -cpus {threads}; "
        # Concatenate whatever hint files were produced. coverage.hints is
        # always present when there is read coverage; junctions.hints only
        # exists when spliced reads were found.
        "cat $(basename {params.bam_abs}).coverage.hints \\\n"
        "    $(basename {params.bam_abs}).junctions.hints 2>/dev/null \\\n"
        "    > {params.hints_abs} || "
        "cp $(basename {params.bam_abs}).coverage.hints {params.hints_abs}"
