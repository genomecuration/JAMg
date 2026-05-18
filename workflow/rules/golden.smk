# Golden gene set: high-confidence PASA-derived gene models used to train
# the predictors (Augustus) and to seed EVM. Runs the v1
# bin/prepare_golden_genes_for_predictors.pl until Phase 4 lands the
# rewritten bin/prepare_golden_genes.pl + PerlLib/Golden/*.pm modules.

import os as _os

_GOLDEN_DIR        = f"{OUTDIR}/golden"
_GENOME_ABS        = _os.path.abspath(config["genome"])
_INTRON            = config["max_intron"]


rule golden_genes:
    input:
        softmasked     = rules.repeats_merge.output.soft,
        pasa_genome    = rules.pasa_compare_transdecoder.output.transdecoder,
        pasa_assembly  = rules.pasa_compare_transdecoder.output.pasa_assemblies,
    output:
        gff = f"{_GOLDEN_DIR}/final_golden_genes.gff3.nr.golden.gff3",
    params:
        genome_abs     = _GENOME_ABS,
        softmasked_abs = lambda _wc, input: _os.path.abspath(input.softmasked),
        pasa_genome_abs = lambda _wc, input: _os.path.abspath(input.pasa_genome),
        pasa_assembly_abs = lambda _wc, input: _os.path.abspath(input.pasa_assembly),
        golden_dir_abs = lambda _wc, output: _os.path.dirname(_os.path.abspath(output.gff)),
        intron         = _INTRON,
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "mkdir -p {params.golden_dir_abs} && cd {params.golden_dir_abs} && "
        "prepare_golden_genes_for_predictors.pl "
        "    --genome {params.genome_abs} "
        "    --softmasked {params.softmasked_abs} "
        "    --pasa_genome {params.pasa_genome_abs} "
        "    --pasa_assembly {params.pasa_assembly_abs} "
        "    --intron {params.intron} "
        "    --threads {threads}"


rule golden_hints:
    # gff2hints.pl writes <input>.hints which already matches output.hints
    # (declared as gff path + '.hints'); no rename needed.
    input:
        gff = rules.golden_genes.output.gff,
    output:
        hints = f"{_GOLDEN_DIR}/final_golden_genes.gff3.nr.golden.gff3.hints",
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "gff2hints.pl {input.gff}"
