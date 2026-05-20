# 4-pass RepeatMasker (contamination_check + general + RNA-specific + species-
# specific) -> merged GFF3 -> soft/hardmasked FASTA. contamination_check is
# mandatory pre-pass QC; its output is NOT in the merge.
#
# Each pass cd's into its own output dir to keep RepeatMasker's many sibling
# files isolated. Every external path referenced inside the shell is
# absolutized via os.path.abspath() at rule-parse time, so the cd does not
# break path resolution.

import os as _os

RG = f"{OUTDIR}/repeats/general/{GENOME_BASENAME}"
RR = f"{OUTDIR}/repeats/rna-specific/{GENOME_BASENAME}"
RS = f"{OUTDIR}/repeats/species-specific/{GENOME_BASENAME}"
RC = f"{OUTDIR}/repeats/contamination-check/{GENOME_BASENAME}"

# os.path.abspath() resolves config paths against the CWD snakemake is
# invoked from. The CLI runs snakemake from the repo root, so relative paths
# in the user's config (e.g. genome: test_suite/mini-genome.fasta) resolve
# correctly. Users invoking snakemake directly from a different CWD should
# write absolute paths in the config.
_GENOME_ABS  = _os.path.abspath(config["genome"])
_RNA_LIB_ABS = _os.path.abspath(config["repeats"]["rna_lib"])
_SPP_LIB_ABS = _os.path.abspath(config["repeats"]["species_lib"])


rule repeats_contamination_check:
    # RepeatMasker may report "no repetitive sequences detected" and emit only
    # the .out summary file, not the .gff. Touch the declared .gff afterwards
    # so the rule produces its expected output regardless. This pass is QC
    # only; its output is NOT merged into all_repeat_masks.gff3.
    input:
        genome = config["genome"],
    output:
        gff = RC + ".out.gff",
    params:
        cat        = config["repeats"]["species_category"],
        genome_abs = _GENOME_ABS,
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "mkdir -p $(dirname {output.gff}) && cd $(dirname {output.gff}) && "
        "ln -sf {params.genome_abs} . && "
        # famdb.py uses #!/usr/bin/env python3 which resolves to /usr/bin/python3
        # (Python 3.13, no h5py) in the container because /usr/bin precedes
        # /opt/conda/bin on PATH. Prepend conda so famdb.py gets conda python3
        # (3.12 + h5py) for the -species lookup.
        "(PATH=/opt/conda/bin:$PATH RepeatMasker -s -excln -pa {threads} -gff -xsmall -gccalc -frag 5000000000 "
        "   -e ncbi -is_only -species '{params.cat}' $(basename {params.genome_abs}) "
        " || true); "
        "touch $(basename {output.gff})"


rule repeats_general:
    input:
        genome = config["genome"],
    output:
        gff = RG + ".out.gff",
    params:
        cat        = config["repeats"]["species_category"],
        genome_abs = _GENOME_ABS,
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "mkdir -p $(dirname {output.gff}) && cd $(dirname {output.gff}) && "
        "ln -sf {params.genome_abs} . && "
        "(PATH=/opt/conda/bin:$PATH RepeatMasker -s -excln -pa {threads} -gff -xsmall -gccalc "
        "   -frag 500000 -e ncbi -species '{params.cat}' $(basename {params.genome_abs}) "
        " || true); "
        "touch $(basename {output.gff})"


rule repeats_rna_specific:
    input:
        genome = config["genome"],
    output:
        gff = RR + ".trim.out.gff",
    params:
        lib_abs    = _RNA_LIB_ABS,
        genome_abs = _GENOME_ABS,
        gff_abs    = lambda wc, output: _os.path.abspath(output.gff),
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        # `awk !/.../` consolidates the three grep -v passes into one and
        # avoids the pipefail-vs-empty-match grep failure mode.
        "mkdir -p $(dirname {output.gff}) && cd $(dirname {output.gff}) && "
        "ln -sf {params.genome_abs} . && "
        "(RepeatMasker -s -excln -pa {threads} -gff -xsmall -gccalc -frag 500000 "
        "   -no_is -nolow -lib {params.lib_abs} -e ncbi $(basename {params.genome_abs}) "
        " || true); "
        "if [ -s $(basename {params.genome_abs}).out ] "
        "   && grep -qE '^[[:space:]]*[0-9]' $(basename {params.genome_abs}).out; then "
        "    awk '!/ Unknown/ && !/ Simple_repeat/ && !/ Low_complexity/'"
        "        $(basename {params.genome_abs}).out > $(basename {params.genome_abs}).trim.out && "
        "    rmOutToGFF3.pl $(basename {params.genome_abs}).trim.out > {params.gff_abs}; "
        "else "
        "    : > {params.gff_abs}; "
        "fi"


rule repeats_species_specific:
    input:
        genome = config["genome"],
    output:
        gff = RS + ".trim.out.gff",
    params:
        lib_abs    = _SPP_LIB_ABS,
        genome_abs = _GENOME_ABS,
        gff_abs    = lambda wc, output: _os.path.abspath(output.gff),
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "mkdir -p $(dirname {output.gff}) && cd $(dirname {output.gff}) && "
        "ln -sf {params.genome_abs} . && "
        "(RepeatMasker -s -excln -pa {threads} -gff -xsmall -gccalc -frag 500000 "
        "   -lib {params.lib_abs} -e ncbi $(basename {params.genome_abs}) "
        " || true); "
        "if [ -s $(basename {params.genome_abs}).out ] "
        "   && grep -qE '^[[:space:]]*[0-9]' $(basename {params.genome_abs}).out; then "
        "    awk '!/ Unknown/ && !/ Simple_repeat/ && !/ Low_complexity/'"
        "        $(basename {params.genome_abs}).out > $(basename {params.genome_abs}).trim.out && "
        "    rmOutToGFF3.pl $(basename {params.genome_abs}).trim.out > {params.gff_abs}; "
        "else "
        "    : > {params.gff_abs}; "
        "fi"


rule repeats_merge:
    # contamination_check output is intentionally excluded from the merge but
    # kept as an input edge to force the QC pass to run before the merge.
    # makeblastdb writes ~9 sibling files: snakemake watches the touch-sentinel
    # instead so the DAG edge is well-defined.
    input:
        contam  = rules.repeats_contamination_check.output.gff,
        general = rules.repeats_general.output.gff,
        rna     = rules.repeats_rna_specific.output.gff,
        spp     = rules.repeats_species_specific.output.gff,
        genome  = config["genome"],
    output:
        merged       = f"{OUTDIR}/repeats/all_repeat_masks.gff3",
        soft         = f"{OUTDIR}/repeats/{GENOME_BASENAME}.softmasked",
        soft_fai     = f"{OUTDIR}/repeats/{GENOME_BASENAME}.softmasked.fai",
        hard         = f"{OUTDIR}/repeats/{GENOME_BASENAME}.hardmasked",
        hard_fai     = f"{OUTDIR}/repeats/{GENOME_BASENAME}.hardmasked.fai",
        blastdb_done = touch(f"{OUTDIR}/repeats/{GENOME_BASENAME}.softmasked.blastdb.done"),
        hints        = f"{OUTDIR}/repeats/all_repeat_masks.gff3.hints",
    container: "containers/jamg.sif"
    resources:
        mem_mb = 4000,
    shell:
        "grep -hv '^#' {input.general} {input.rna} {input.spp} | grep -v '^$' | "
        "  sort -k1,1 -k4,5n | awk '!seen[$1,$4,$5,$7]++' > {output.merged} && "
        "bedtools maskfasta -fi {input.genome} -fo {output.soft} -bed {output.merged} -soft && "
        "bedtools maskfasta -fi {input.genome} -fo {output.hard} -bed {output.merged} && "
        "samtools faidx {output.soft} && samtools faidx {output.hard} && "
        "makeblastdb -dbtype nucl -in {output.soft} -parse_seqids -hash_index "
        f"           -title '{GENOME_BASENAME}' && "
        # repeatmasker2hints.pl writes to <input>.hints which already matches
        # output.hints (declared as merged + '.hints'), so no rename is needed.
        # Guard against empty merged GFF (e.g. clean genome with no library hits):
        # the script exits 255 on empty input, so we touch the hints file instead.
        "if [ -s {output.merged} ]; then repeatmasker2hints.pl {output.merged}; "
        "else touch {output.hints}; fi"
