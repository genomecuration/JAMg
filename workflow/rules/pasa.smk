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
        cln                = f"{_PASA_DIR}/transcripts.fasta.cln",
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
        # PASA's -C (create) refuses to run against an existing populated
        # sqlite ("table URL_templates already exists"). /dev/shm survives
        # across snakemake re-runs on the same compute node, so a stale
        # sqlite from a prior partial run can block this one. Always start
        # fresh; snakemake's job-level cache on pasa.sqlite.align.bz2
        # already prevents re-running this rule when its output is current.
        "rm -f {params.db_path} && "
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
        pasa_assemblies   = f"{_PASA_DIR}/pasa_assemblies.gff3",
        transdecoder      = f"{_PASA_DIR}/pasa.transdecoder.genome.gff3",
        # 4 additional pasa_asmbls_to_training_set.dbi products. Golden
        # genes (§3e) needs all 5 PASA inputs: --pasa_gff, --pasa_genome,
        # --pasa_assembly, --pasa_peptides, --pasa_cds (per the script's
        # check_for_options in bin/prepare_golden_genes_for_predictors.pl
        # line 3415-3426; missing any of them triggers pod2usage).
        transdecoder_gff  = f"{_PASA_DIR}/pasa.transdecoder.gff3",
        transdecoder_pep  = f"{_PASA_DIR}/pasa.transdecoder.pep",
        transdecoder_cds  = f"{_PASA_DIR}/pasa.transdecoder.cds",
        assemblies_fasta  = f"{_PASA_DIR}/pasa.assemblies.fasta",
        polya             = f"{_PASA_DIR}/polyAsites.fasta",
        pass1             = f"{_PASA_DIR}/pasa.sqlite.pass1.bz2",
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
        # /dev/shm persists across runs on the same compute node, so a stale
        # sqlite from a prior incomplete run could be present. Restore the
        # in-repo snapshot unconditionally (overwriting any stale state) so
        # -R resumes against the correct state.
        "rm -f {params.db_path} && "
        "bunzip2 -c pasa.sqlite.align.bz2 > {params.db_path} && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config -R "
        "  -g {params.genome_abs} "
        "  --MAX_INTRON_LENGTH {params.max_intron} "
        "  --ALIGNERS gmap,minimap2 "
        "  --TRANSDECODER --CPU {threads} "
        "  -T -t transcripts.fasta.clean -u transcripts.fasta "
        "  --TDN tdn.accs && "
        # PASA's --TRANSDECODER flag runs TransDecoder on the input
        # transcripts (produces transcripts.fasta.clean.transdecoder.gff3)
        # but does NOT call the script that converts the assembled-transcript
        # ORFs into genome coordinates. That requires
        # pasa_asmbls_to_training_set.dbi as a separate post-step; it
        # produces pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3,
        # which the §3i evm_tag_pasa_transdecoder rule consumes (renamed
        # below to pasa.transdecoder.genome.gff3).
        "$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi "
        "  --pasa_transcripts_fasta pasa.sqlite.assemblies.fasta "
        "  --pasa_transcripts_gff3 pasa.sqlite.pasa_assemblies.gff3 && "
        # PASA emits outputs with prefixes derived from the db basename; the
        # standardised names below are the actual file basenames produced
        # by Launch_PASA_pipeline.pl with our config. Snapshot the sqlite
        # state IMMEDIATELY after the compare; the next compare overwrites it.
        # Each rename loop: stop at the first match (break) so a stale
        # leftover file from a prior aborted run does not silently overwrite
        # the freshly-produced output. If no match exists, FATAL loudly
        # rather than letting `set -e` swallow the for-loop's exit code.
        "bzip2 -c {params.db_path} > pasa.sqlite.pass1.bz2 && "
        "for f in *.assemblies.fasta.transdecoder.genome.gff3; do "
        "    if [ -e \"$f\" ]; then cp \"$f\" pasa.transdecoder.genome.gff3; break; fi; "
        "done && "
        "test -s pasa.transdecoder.genome.gff3 || "
        "    {{ echo 'FATAL: pasa_asmbls_to_training_set.dbi produced no *.assemblies.fasta.transdecoder.genome.gff3' >&2; exit 1; }} && "
        # Rename the four additional pasa_asmbls_to_training_set.dbi
        # products into stable names (snakemake outputs declared above).
        # The .gff3, .pep, .cds are the transcript-coord TransDecoder
        # products; .assemblies.fasta is PASA's assembled transcripts.
        "for f in *.assemblies.fasta.transdecoder.gff3; do "
        "    if [ -e \"$f\" ]; then cp \"$f\" pasa.transdecoder.gff3; break; fi; "
        "done && "
        "test -s pasa.transdecoder.gff3 || "
        "    {{ echo 'FATAL: missing *.assemblies.fasta.transdecoder.gff3' >&2; exit 1; }} && "
        "for f in *.assemblies.fasta.transdecoder.pep; do "
        "    if [ -e \"$f\" ]; then cp \"$f\" pasa.transdecoder.pep; break; fi; "
        "done && "
        "test -s pasa.transdecoder.pep || "
        "    {{ echo 'FATAL: missing *.assemblies.fasta.transdecoder.pep' >&2; exit 1; }} && "
        "for f in *.assemblies.fasta.transdecoder.cds; do "
        "    if [ -e \"$f\" ]; then cp \"$f\" pasa.transdecoder.cds; break; fi; "
        "done && "
        "test -s pasa.transdecoder.cds || "
        "    {{ echo 'FATAL: missing *.assemblies.fasta.transdecoder.cds' >&2; exit 1; }} && "
        "for f in *.assemblies.fasta; do "
        # Match only the bare assemblies fasta, not *.assemblies.fasta.<ext>
        "    case \"$f\" in *.transdecoder*|*.fai|*.gmap|*.mm2) continue;; esac; "
        "    if [ -e \"$f\" ]; then cp \"$f\" pasa.assemblies.fasta; break; fi; "
        "done && "
        "test -s pasa.assemblies.fasta || "
        "    {{ echo 'FATAL: missing PASA *.assemblies.fasta' >&2; exit 1; }} && "
        "for f in *.pasa_assemblies.gff3; do "
        "    if [ -e \"$f\" ]; then cp \"$f\" pasa_assemblies.gff3; break; fi; "
        "done && "
        "test -s pasa_assemblies.gff3 || "
        "    {{ echo 'FATAL: PASA produced no *.pasa_assemblies.gff3' >&2; exit 1; }} && "
        "for f in *.polyAsites.fasta; do "
        "    if [ -e \"$f\" ]; then cp \"$f\" polyAsites.fasta; break; fi; "
        "done && "
        "test -e polyAsites.fasta || "
        "    {{ echo 'FATAL: PASA produced no *.polyAsites.fasta' >&2; exit 1; }}"


# NOTE: pasa_compare_load_1 / pasa_compare_load_2 were originally specified
# in plan §3d but moved out after testing. PASA's `-A -L --annots <gff>` mode
# updates an EXISTING gene annotation with PASA's transcript evidence (per
# PASA upstream wiki, PASA_genome_annotation.md inside pasa.sif: "The PASA
# software can update any preexisting set of protein-coding gene
# annotations"). Pre-EVM there is NO existing annotation; passing PASA's own
# assemblies as `--annots` produces a 0-byte gene_structures_post_PASA_
# updates output because PASA is asked to compare its own evidence against
# itself. Empirically confirmed 2026-05-19 on the mini-fixture.
#
# The two compare-load passes belong post-EVM (plan §3j step 2, "Post-EVM
# PASA compares (x2)") where EVM-combined gene models ARE the preexisting
# annotation. They will be implemented in ogs.smk when §3j lands.


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
        # pasapolyA2hints.pl lives under /opt/jamg/share/Augustus/scripts/
        # inside jamg.sif (Augustus auxiliary script, not in /opt/jamg/bin
        # which is on PATH). Its upstream shebang is broken (`#!/usr/bin
        # env`, missing slash); jamg.def %post patches it in place at SIF
        # build time. Invoke by absolute path.
        "/opt/jamg/share/Augustus/scripts/pasapolyA2hints.pl {input.polya} > {output.polya_hints} && "
        "gff2hints.pl {input.pasa_assemblies}"
