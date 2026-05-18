# PASA pre-EVM pipeline: setup -> align -> compare_transdecoder -> two
# compare-load passes. Produces pass1/pass2/pass3 sqlite snapshots; pass3
# is what 3j post-EVM consumes. Runs inside containers/pasa.sif (upstream
# Docker pasapipeline/pasapipeline:2.5.3); jamg.sif handles the seqclean
# + concat preparation step where pasa is not yet needed.

import os as _os

_PASA_DIR        = f"{OUTDIR}/pasa"
_GENOME_ABS      = _os.path.abspath(config["genome"])
_PASA_TEMPLATE   = _os.path.abspath(config["pasa"]["config_template"])
_PASA_MAX_INTRON = config["pasa"]["max_intron_length"]


def _abspath_or_none(value):
    """Absolutize a config path so the shell can cd before reading it.
    None/empty → None (the optional input was not configured)."""
    return _os.path.abspath(value) if value else None


_TRINITY_TDN_ABS = _abspath_or_none(config.get("trinity_denovo"))
_LONGREADS_ABS   = _abspath_or_none(config.get("longreads"))


def _pasa_db_path():
    """SQLite path inside pasa.tmp_dir (typically /dev/shm) for per-run scratch."""
    tmp = config["pasa"].get("tmp_dir", "/dev/shm")
    user = _os.environ.get("USER", "jamg")
    return f"{tmp}/{user}/pasa.sqlite"


rule pasa_setup_db:
    # Builds transcripts.fasta (concat of Trinity-GG, optional Trinity-TDN,
    # optional longreads), tdn.accs, transcripts.fasta.clean (via seqclean),
    # and renders alignAssembly.config from the template.
    input:
        trinity_gg = rules.tgg_trinity.output.fasta,
    output:
        transcripts        = f"{_PASA_DIR}/transcripts.fasta",
        transcripts_clean  = f"{_PASA_DIR}/transcripts.fasta.clean",
        tdn_accs           = f"{_PASA_DIR}/tdn.accs",
        config_rendered    = f"{_PASA_DIR}/alignAssembly.config",
    params:
        # Use sentinel string 'NONE' for unconfigured optionals; the shell
        # branches on `!= "NONE"` to distinguish "not configured" (skip
        # silently) from "configured but file missing" (fail loudly).
        trinity_tdn   = _TRINITY_TDN_ABS or "NONE",
        longreads     = _LONGREADS_ABS or "NONE",
        template_abs  = _PASA_TEMPLATE,
        db_path       = _pasa_db_path(),
        pasa_dir_abs  = lambda _wc, output: _os.path.dirname(_os.path.abspath(output.transcripts)),
        trinity_gg_abs = lambda _wc, input: _os.path.abspath(input.trinity_gg),
    container: "containers/pasa.sif"
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        "mkdir -p {params.pasa_dir_abs} && cd {params.pasa_dir_abs} && "
        # Build transcripts.fasta from up to three sources. Trinity-GG is
        # always present; Trinity-TDN and longreads are optional.
        "cat {params.trinity_gg_abs} > transcripts.fasta && "
        "if [ '{params.trinity_tdn}' != 'NONE' ]; then "
        "    if [ ! -s '{params.trinity_tdn}' ]; then "
        "        echo 'FATAL: trinity_denovo path \"{params.trinity_tdn}\" does not exist or is empty' >&2; "
        "        exit 1; "
        "    fi; "
        "    cat '{params.trinity_tdn}' >> transcripts.fasta; "
        "    accession_extractor.pl '{params.trinity_tdn}' > tdn.accs; "
        "else "
        "    : > tdn.accs; "
        "fi && "
        "if [ '{params.longreads}' != 'NONE' ]; then "
        "    if [ ! -s '{params.longreads}' ]; then "
        "        echo 'FATAL: longreads path \"{params.longreads}\" does not exist or is empty' >&2; "
        "        exit 1; "
        "    fi; "
        "    cat '{params.longreads}' >> transcripts.fasta; "
        "fi && "
        "seqclean transcripts.fasta && "
        "sed 's|<__DATABASE__>|{params.db_path}|' {params.template_abs} > alignAssembly.config"


rule pasa_align:
    # Launch_PASA_pipeline.pl -C (create db + align). No PASA pass-snapshot
    # here (that starts at compare_transdecoder), but we DO bzip2-snapshot
    # the sqlite into _PASA_DIR so a resumed run does not consume a /dev/shm
    # state that was reset between runs.
    input:
        cfg     = rules.pasa_setup_db.output.config_rendered,
        clean   = rules.pasa_setup_db.output.transcripts_clean,
        full    = rules.pasa_setup_db.output.transcripts,
        tdn     = rules.pasa_setup_db.output.tdn_accs,
        genome  = config["genome"],
    output:
        align_snapshot = f"{_PASA_DIR}/pasa.sqlite.align.bz2",
    params:
        genome_abs    = _GENOME_ABS,
        max_intron    = _PASA_MAX_INTRON,
        pasa_dir_abs  = lambda _wc, input: _os.path.dirname(_os.path.abspath(input.cfg)),
        db_path       = _pasa_db_path(),
    container: "containers/pasa.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "cd {params.pasa_dir_abs} && "
        "mkdir -p $(dirname {params.db_path}) && "
        # If a prior align ran and we have a snapshot but /dev/shm was
        # cleaned, restore the sqlite from the snapshot before invoking PASA.
        "if [ ! -s {params.db_path} ] && [ -s pasa.sqlite.align.bz2 ]; then "
        "    bunzip2 -c pasa.sqlite.align.bz2 > {params.db_path}; "
        "fi && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config -C "
        "  -g {params.genome_abs} "
        "  -T -t transcripts.fasta.clean -u transcripts.fasta "
        "  --TDN tdn.accs "
        "  --MAX_INTRON_LENGTH {params.max_intron} "
        "  --ALIGNERS gmap,minimap2 "
        "  --CPU {threads} && "
        "bzip2 -c {params.db_path} > pasa.sqlite.align.bz2"


rule pasa_compare_transdecoder:
    # Launch_PASA_pipeline.pl -R --TRANSDECODER. Emits pass1.bz2.
    input:
        align_snapshot = rules.pasa_align.output.align_snapshot,
        cfg            = rules.pasa_setup_db.output.config_rendered,
        clean          = rules.pasa_setup_db.output.transcripts_clean,
        full           = rules.pasa_setup_db.output.transcripts,
        tdn            = rules.pasa_setup_db.output.tdn_accs,
        genome         = config["genome"],
    output:
        pasa_assemblies = f"{_PASA_DIR}/pasa_assemblies.gff3",
        transdecoder    = f"{_PASA_DIR}/pasa.transdecoder.genome.gff3",
        polya           = f"{_PASA_DIR}/polyAsites.fasta",
        pass1           = f"{_PASA_DIR}/pasa.sqlite.pass1.bz2",
    params:
        genome_abs    = _GENOME_ABS,
        max_intron    = _PASA_MAX_INTRON,
        pasa_dir_abs  = lambda _wc, input: _os.path.dirname(_os.path.abspath(input.cfg)),
        db_path       = _pasa_db_path(),
    container: "containers/pasa.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        "cd {params.pasa_dir_abs} && "
        "mkdir -p $(dirname {params.db_path}) && "
        "if [ ! -s {params.db_path} ] && [ -s pasa.sqlite.align.bz2 ]; then "
        "    bunzip2 -c pasa.sqlite.align.bz2 > {params.db_path}; "
        "fi && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config -R "
        "  -g {params.genome_abs} "
        "  --MAX_INTRON_LENGTH {params.max_intron} "
        "  --ALIGNERS gmap,minimap2 "
        "  --TRANSDECODER --CPU {threads} "
        "  -T -t transcripts.fasta.clean -u transcripts.fasta "
        "  --TDN tdn.accs && "
        # PASA emits outputs with prefixes derived from the db basename; the
        # standardised names below are the actual file basenames produced
        # by Launch_PASA_pipeline.pl with our config. Snapshot the sqlite
        # state IMMEDIATELY after the compare; the next compare overwrites it.
        "bzip2 -c {params.db_path} > pasa.sqlite.pass1.bz2 && "
        "for f in *.assemblies.fasta.transdecoder.genome.gff3; do "
        "    [ -e \"$f\" ] && cp \"$f\" pasa.transdecoder.genome.gff3; "
        "done && "
        "for f in *.pasa_assemblies.gff3; do "
        "    [ -e \"$f\" ] && cp \"$f\" pasa_assemblies.gff3; "
        "done && "
        "for f in *.polyAsites.fasta; do "
        "    [ -e \"$f\" ] && cp \"$f\" polyAsites.fasta && break; "
        "done"


rule pasa_compare_load_1:
    # Launch_PASA_pipeline.pl -A -L --annots pasa_assemblies.gff3.
    # Emits pass2.bz2 + the first gene_structures_post_PASA_updates.<round>.gff3.
    input:
        pasa_assemblies = rules.pasa_compare_transdecoder.output.pasa_assemblies,
        cfg             = rules.pasa_setup_db.output.config_rendered,
        clean           = rules.pasa_setup_db.output.transcripts_clean,
        genome          = config["genome"],
    output:
        pass2     = f"{_PASA_DIR}/pasa.sqlite.pass2.bz2",
        updates_1 = f"{_PASA_DIR}/gene_structures_post_PASA_updates.round1.gff3",
    params:
        genome_abs    = _GENOME_ABS,
        pasa_dir_abs  = lambda _wc, input: _os.path.dirname(_os.path.abspath(input.cfg)),
        db_path       = _pasa_db_path(),
    container: "containers/pasa.sif"
    threads: THREADS
    resources:
        mem_mb = 8000,
    shell:
        "cd {params.pasa_dir_abs} && "
        "mkdir -p $(dirname {params.db_path}) && "
        "if [ ! -s {params.db_path} ] && [ -s pasa.sqlite.pass1.bz2 ]; then "
        "    bunzip2 -c pasa.sqlite.pass1.bz2 > {params.db_path}; "
        "fi && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config "
        "  -g {params.genome_abs} -t transcripts.fasta.clean "
        "  -A -L --annots pasa_assemblies.gff3 && "
        "bzip2 -c {params.db_path} > pasa.sqlite.pass2.bz2 && "
        # PASA emits multiple round files in a single compare-load; pick the
        # first numerically (round1.gff3 lexically before round2.gff3, valid
        # since round-N stays single-digit for our 2-pass pipeline).
        "first=$(ls *.gene_structures_post_PASA_updates.*.gff3 2>/dev/null | sort -V | head -1) && "
        "[ -n \"$first\" ] && cp \"$first\" gene_structures_post_PASA_updates.round1.gff3"


rule pasa_compare_load_2:
    # Launch_PASA_pipeline.pl -A -L --annots <previous gene_structures>.
    # Emits pass3.bz2 + the second gene_structures_post_PASA_updates.<round>.gff3.
    # pass3.bz2 is the snapshot §3j consumes via bunzip2 -fkc.
    input:
        updates_1 = rules.pasa_compare_load_1.output.updates_1,
        cfg       = rules.pasa_setup_db.output.config_rendered,
        clean     = rules.pasa_setup_db.output.transcripts_clean,
        genome    = config["genome"],
    output:
        pass3     = f"{_PASA_DIR}/pasa.sqlite.pass3.bz2",
        updates_2 = f"{_PASA_DIR}/gene_structures_post_PASA_updates.round2.gff3",
    params:
        genome_abs    = _GENOME_ABS,
        pasa_dir_abs  = lambda _wc, input: _os.path.dirname(_os.path.abspath(input.cfg)),
        db_path       = _pasa_db_path(),
    container: "containers/pasa.sif"
    threads: THREADS
    resources:
        mem_mb = 8000,
    shell:
        "cd {params.pasa_dir_abs} && "
        "mkdir -p $(dirname {params.db_path}) && "
        "if [ ! -s {params.db_path} ] && [ -s pasa.sqlite.pass2.bz2 ]; then "
        "    bunzip2 -c pasa.sqlite.pass2.bz2 > {params.db_path}; "
        "fi && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config "
        "  -g {params.genome_abs} -t transcripts.fasta.clean "
        "  -A -L --annots gene_structures_post_PASA_updates.round1.gff3 && "
        "bzip2 -c {params.db_path} > pasa.sqlite.pass3.bz2 && "
        # Pick the latest round file by version-sort; numeric tail (round2,
        # round3, ...) is robust to round10+ ordering issues.
        "last=$(ls *.gene_structures_post_PASA_updates.*.gff3 2>/dev/null | sort -V | tail -1) && "
        "[ -n \"$last\" ] && cp \"$last\" gene_structures_post_PASA_updates.round2.gff3"


rule pasa_hints:
    # Augustus hints from PASA's polyA sites + assembly GFF.
    # pasapolyA2hints.pl and gff2hints.pl live in jamg.sif.
    input:
        polya           = rules.pasa_compare_transdecoder.output.polya,
        pasa_assemblies = rules.pasa_compare_transdecoder.output.pasa_assemblies,
    output:
        polya_hints   = f"{_PASA_DIR}/polyAsites.hints",
        assembly_hints = f"{_PASA_DIR}/pasa_assemblies.gff3.hints",
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "pasapolyA2hints.pl {input.polya} > {output.polya_hints} && "
        "gff2hints.pl {input.pasa_assemblies}"
