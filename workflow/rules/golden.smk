# Golden gene set: high-confidence PASA-derived gene models used to train
# the predictors (Augustus) and to seed EVM. Runs the v1
# bin/prepare_golden_genes_for_predictors.pl until Phase 4 lands the
# rewritten bin/prepare_golden_genes.pl + PerlLib/Golden/*.pm modules.

import os as _os

_GOLDEN_DIR        = f"{OUTDIR}/golden"
_GENOME_ABS        = _os.path.abspath(config["genome"])
_INTRON            = config["max_intron"]


rule golden_genes:
    # bin/prepare_golden_genes_for_predictors.pl::check_for_options (line
    # 3415-3426) requires ALL FIVE PASA inputs to be supplied together
    # AND non-empty; passing only --pasa_genome + --pasa_assembly (the
    # initial naive choice) triggers pod2usage and exit 2.
    input:
        softmasked       = rules.repeats_merge.output.soft,
        pasa_gff         = rules.pasa_compare_transdecoder.output.transdecoder_gff,
        pasa_genome      = rules.pasa_compare_transdecoder.output.transdecoder,
        pasa_assembly    = rules.pasa_compare_transdecoder.output.assemblies_fasta,
        pasa_peptides    = rules.pasa_compare_transdecoder.output.transdecoder_pep,
        pasa_cds         = rules.pasa_compare_transdecoder.output.transdecoder_cds,
    output:
        gff = f"{_GOLDEN_DIR}/final_golden_genes.gff3.nr.golden.gff3",
    params:
        genome_abs         = _GENOME_ABS,
        softmasked_abs     = lambda _wc, input: _os.path.abspath(input.softmasked),
        pasa_gff_abs       = lambda _wc, input: _os.path.abspath(input.pasa_gff),
        pasa_genome_abs    = lambda _wc, input: _os.path.abspath(input.pasa_genome),
        pasa_assembly_abs  = lambda _wc, input: _os.path.abspath(input.pasa_assembly),
        pasa_peptides_abs  = lambda _wc, input: _os.path.abspath(input.pasa_peptides),
        pasa_cds_abs       = lambda _wc, input: _os.path.abspath(input.pasa_cds),
        golden_dir_abs     = lambda _wc, output: _os.path.dirname(_os.path.abspath(output.gff)),
        intron             = _INTRON,
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        # The v1 script defaults --gmap_dir to $RealBin/../databases/gmap/
        # which resolves to /opt/jamg/databases/gmap/ inside the SIF (read-
        # only). Point it at a writable subdirectory of the rule's workdir
        # so gmap_build can write its database where snakemake has bind-
        # mounted write access.
        #
        # The script's check_augustus derives augustus_dir =
        # dirname(dirname(augustus_exec)). With augustus at /opt/jamg/bin/
        # that gives /opt/jamg/, then it looks for /opt/jamg/scripts/...
        # which doesn't exist. Pass --augustus explicitly pointing at the
        # actual Augustus install location (jamg.sif puts it under
        # /opt/jamg/share/Augustus/ with bin/ + scripts/ subdirs).
        "mkdir -p {params.golden_dir_abs}/gmap && cd {params.golden_dir_abs} && "
        "prepare_golden_genes_for_predictors.pl "
        "    --genome {params.genome_abs} "
        "    --softmasked {params.softmasked_abs} "
        "    --pasa_gff {params.pasa_gff_abs} "
        "    --pasa_genome {params.pasa_genome_abs} "
        "    --pasa_assembly {params.pasa_assembly_abs} "
        "    --pasa_peptides {params.pasa_peptides_abs} "
        "    --pasa_cds {params.pasa_cds_abs} "
        "    --gmap_dir {params.golden_dir_abs}/gmap "
        "    --augustus /opt/jamg/share/Augustus "
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
