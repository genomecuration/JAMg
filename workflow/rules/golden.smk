# Golden gene set: high-confidence PASA-derived gene models used to train
# the predictors (Augustus) and to seed EVM. Wraps
# bin/prepare_golden_genes.pl + PerlLib/Golden/{Alignment,Filter,Augustus}.pm.

import os as _os

_GOLDEN_DIR        = f"{OUTDIR}/golden"
_GENOME_ABS        = _os.path.abspath(config["genome"])
_INTRON            = config["max_intron"]


rule golden_genes:
    # bin/prepare_golden_genes.pl::check_for_options requires ALL FIVE PASA
    # inputs together AND non-empty; passing only --pasa_genome +
    # --pasa_assembly triggers pod2usage and exit 2.
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
        repo_bin           = _os.path.abspath(_os.path.join(workflow.basedir, "..", "bin")),
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        # --gmap_dir points at a writable subdirectory of the rule's workdir
        # (the SIF's default /opt/jamg/databases/gmap/ is read-only). The
        # driver's check_augustus derives augustus_dir =
        # dirname(dirname(augustus_exec)) which gives /opt/jamg/, where the
        # auxiliary scripts/ subtree does not exist. Pass --augustus
        # explicitly pointing at the actual install location
        # (jamg.sif: /opt/jamg/share/Augustus/{bin,scripts}/).
        # bin/ lives on the host (bind-mounted into the SIF at the repo
        # root); /opt/jamg/bin/ holds the upstream-shipped helpers. Add
        # repo bin/ to PATH so prepare_golden_genes.pl resolves.
        "export PATH={params.repo_bin}:$PATH && "
        "mkdir -p {params.golden_dir_abs}/gmap && "
        "prepare_golden_genes.pl "
        "    --outdir {params.golden_dir_abs} "
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
