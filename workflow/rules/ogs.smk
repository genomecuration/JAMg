# Post-EVM Official Gene Set (OGS) build. Plan §3j.
#
# Pipeline shape:
#   1  Rescue Augustus-only genes that EVM dropped (trim_overlap_gff3.py)
#   2  Two PASA -A -L compare passes against the EVM-combined model
#      (the compare-load passes that were intentionally NOT placed in §3d;
#       see workflow/rules/pasa.smk's NOTE near pasa_compare_transdecoder).
#   3  ID rewrite + simple-feature post-processing (create_features_from_gff3.pl)
#   4  Validate + sort (fix_pasa_evm_gff.pl + gt gff3 -sort -tidy)
#   5  Final products (gffread → OGS.gff3, OGS.{mRNA,CDS,pep}.fasta, OGS.{gtf,bed})

import os as _os

_OGS_DIR        = f"{OUTDIR}/ogs"
_EVM_DIR_REF    = f"{OUTDIR}/evm"
_AUG_DIR_REF    = f"{OUTDIR}/augustus"
_PASA_DIR_REF   = f"{OUTDIR}/pasa"

_REPO_BIN_OGS   = _os.path.abspath(_os.path.join(workflow.basedir, "..", "bin"))
_GENOME_ABS_OGS = _os.path.abspath(config["genome"])
_OGS_CODE       = config.get("code", "JAMG")
_MANUAL_CUR     = config.get("manual_curations")  # optional path, default None


# ---------------------------------------------------------------------------
# §3j.1  ogs_rescue_augustus
# Augustus-only genes (not overlapping EVM models on the same strand) are
# rescued into a missed_genes.gff3 and concatenated with EVM.gff3.
# ---------------------------------------------------------------------------
rule ogs_rescue_augustus:
    input:
        evm      = rules.evm_run.output.evm_gff,
        augustus = rules.evm_tag_augustus.output.gff,
    output:
        missed   = f"{_OGS_DIR}/missed_genes.gff3",
        combined = f"{_OGS_DIR}/EVM.combined.gff3",
    params:
        ogs_dir      = _OGS_DIR,
        repo_bin     = _REPO_BIN_OGS,
        evm_abs      = lambda _wc, input: _os.path.abspath(input.evm),
        augustus_abs = lambda _wc, input: _os.path.abspath(input.augustus),
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "mkdir -p {params.ogs_dir} && cd {params.ogs_dir} && "
        "trim_overlap_gff3.py --file1 {params.augustus_abs} --file2 {params.evm_abs} "
        "    -o missed_genes.gff3 && "
        "cat missed_genes.gff3 {params.evm_abs} > EVM.combined.gff3"


# ---------------------------------------------------------------------------
# §3j.2a  ogs_pasa_compare_1
# First post-EVM PASA compare-load: snapshot pass1, emit
# gene_structures_post_PASA_updates.<round>.gff3 (round-1 update of the
# EVM-combined model).
# ---------------------------------------------------------------------------
rule ogs_pasa_compare_1:
    input:
        combined = rules.ogs_rescue_augustus.output.combined,
        cfg      = rules.pasa_setup_db.output.config_rendered,
        clean    = rules.pasa_setup_db.output.transcripts_clean,
        full     = rules.pasa_setup_db.output.transcripts,
        tdn      = rules.pasa_setup_db.output.tdn_accs,
        genome   = config["genome"],
        pass1    = rules.pasa_compare_transdecoder.output.pass1,
    output:
        snapshot = f"{_OGS_DIR}/pasa.sqlite.evm_pasa1.bz2",
        gff      = f"{_OGS_DIR}/EVM.combined.pasa1.gff3",
    params:
        genome_abs = _GENOME_ABS_OGS,
        max_intron = config.get("max_intron", 70000),
        ogs_dir    = _OGS_DIR,
        pasa_dir   = _PASA_DIR_REF,
        cfg_abs    = lambda _wc, input: _os.path.abspath(input.cfg),
        clean_abs  = lambda _wc, input: _os.path.abspath(input.clean),
        full_abs   = lambda _wc, input: _os.path.abspath(input.full),
        tdn_abs    = lambda _wc, input: _os.path.abspath(input.tdn),
        pass1_abs  = lambda _wc, input: _os.path.abspath(input.pass1),
    container: "containers/pasa.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        # PASA compare-load wants its CWD to hold transcripts.fasta /
        # transcripts.fasta.clean / tdn.accs / alignAssembly.config and
        # produces files there. Stage by symlink into ogs/ so writes don't
        # collide with pasa/.
        "mkdir -p {params.ogs_dir} && cd {params.ogs_dir} && "
        "ln -sf {params.cfg_abs}   alignAssembly.config && "
        "ln -sf {params.clean_abs} transcripts.fasta.clean && "
        "ln -sf {params.full_abs}  transcripts.fasta && "
        "ln -sf {params.tdn_abs}   tdn.accs && "
        # Fresh sqlite from the upstream pass1 snapshot. /dev/shm survives
        # scancel, so always overwrite.
        "rm -f /dev/shm/pasa.sqlite && "
        "bunzip2 -c {params.pass1_abs} > /dev/shm/pasa.sqlite && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config "
        "  -A -L --annots EVM.combined.gff3 "
        "  -g {params.genome_abs} "
        "  --MAX_INTRON_LENGTH {params.max_intron} "
        "  --CPU {threads} "
        "  -T -t transcripts.fasta.clean -u transcripts.fasta "
        "  --TDN tdn.accs && "
        # Snapshot + canonicalise the gene-structures output.
        "bzip2 -c /dev/shm/pasa.sqlite > pasa.sqlite.evm_pasa1.bz2 && "
        "for f in *.gene_structures_post_PASA_updates.*.gff3; do "
        "    if [ -s \"$f\" ]; then cp \"$f\" EVM.combined.pasa1.gff3; break; fi; "
        "done && "
        "test -s EVM.combined.pasa1.gff3 || "
        "    {{ echo 'FATAL: PASA compare_1 produced no non-empty gene_structures_post_PASA_updates' >&2; exit 1; }}"


# ---------------------------------------------------------------------------
# §3j.2b  ogs_pasa_compare_2
# Second compare-load, this time against the round-1 update. Snapshot
# pass2; emit the round-2 gene_structures_post_PASA_updates.<round>.gff3
# which is the OGS-candidate model.
# ---------------------------------------------------------------------------
rule ogs_pasa_compare_2:
    input:
        update_1 = rules.ogs_pasa_compare_1.output.gff,
        snap_1   = rules.ogs_pasa_compare_1.output.snapshot,
        cfg      = rules.pasa_setup_db.output.config_rendered,
        clean    = rules.pasa_setup_db.output.transcripts_clean,
        full     = rules.pasa_setup_db.output.transcripts,
        tdn      = rules.pasa_setup_db.output.tdn_accs,
        genome   = config["genome"],
    output:
        snapshot   = f"{_OGS_DIR}/pasa.sqlite.evm_pasa2.bz2",
        candidate  = f"{_OGS_DIR}/EVM.combined.pasa2.gff3",
    params:
        genome_abs   = _GENOME_ABS_OGS,
        max_intron   = config.get("max_intron", 70000),
        ogs_dir      = _OGS_DIR,
        snap_1_abs   = lambda _wc, input: _os.path.abspath(input.snap_1),
        update_1_abs = lambda _wc, input: _os.path.abspath(input.update_1),
        cfg_abs      = lambda _wc, input: _os.path.abspath(input.cfg),
        clean_abs    = lambda _wc, input: _os.path.abspath(input.clean),
        full_abs     = lambda _wc, input: _os.path.abspath(input.full),
        tdn_abs      = lambda _wc, input: _os.path.abspath(input.tdn),
    container: "containers/pasa.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        # Re-create the per-OGS-dir symlinks: a --forcerun ogs_pasa_compare_2
        # (or a partial clean of {ogs_dir}) wipes compare_1's symlinks and
        # PASA aborts because Launch_PASA_pipeline.pl reads
        # alignAssembly.config / transcripts.fasta / transcripts.fasta.clean /
        # tdn.accs by bare name from cwd. Snakemake doesn't track these
        # without `input:` entries; declare them, re-symlink unconditionally.
        "cd {params.ogs_dir} && "
        "ln -sf {params.cfg_abs}   alignAssembly.config && "
        "ln -sf {params.clean_abs} transcripts.fasta.clean && "
        "ln -sf {params.full_abs}  transcripts.fasta && "
        "ln -sf {params.tdn_abs}   tdn.accs && "
        "rm -f /dev/shm/pasa.sqlite && "
        "bunzip2 -c {params.snap_1_abs} > /dev/shm/pasa.sqlite && "
        "Launch_PASA_pipeline.pl -c alignAssembly.config "
        "  -A -L --annots {params.update_1_abs} "
        "  -g {params.genome_abs} "
        "  --MAX_INTRON_LENGTH {params.max_intron} "
        "  --CPU {threads} "
        "  -T -t transcripts.fasta.clean -u transcripts.fasta "
        "  --TDN tdn.accs && "
        "bzip2 -c /dev/shm/pasa.sqlite > pasa.sqlite.evm_pasa2.bz2 && "
        # The round-2 .gff3 file replaces the round-1 one (PASA reuses the
        # gene_structures_post_PASA_updates name with an incremented round).
        "for f in *.gene_structures_post_PASA_updates.*.gff3; do "
        "    if [ -s \"$f\" ]; then cp \"$f\" EVM.combined.pasa2.gff3; break; fi; "
        "done && "
        "test -s EVM.combined.pasa2.gff3 || "
        "    {{ echo 'FATAL: PASA compare_2 produced no non-empty gene_structures_post_PASA_updates' >&2; exit 1; }}"


# ---------------------------------------------------------------------------
# §3j.3  ogs_create_features
# Renames IDs to <code>.JAM*; injects novel_gene / novel_model markers;
# fixes phase on first CDS; strips Name= attributes; drops scaffolds with
# > N runs longer than 5.
# ---------------------------------------------------------------------------
rule ogs_create_features:
    input:
        candidate = rules.ogs_pasa_compare_2.output.candidate,
    output:
        renamed = f"{_OGS_DIR}/EVM.combined.pasa2.renamed.gff3",
    params:
        ogs_dir    = _OGS_DIR,
        repo_bin   = _REPO_BIN_OGS,
        genome_abs = _GENOME_ABS_OGS,
        code       = _OGS_CODE,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 8000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.ogs_dir} && "
        # create_features_from_gff3.pl emits <input>.renamed.gff3 (and a
        # number of sidecar files); --strip_name (not --strip) per the
        # script's GetOptions.
        "create_features_from_gff3.pl "
        "  -rename "
        "  -genome {params.genome_abs} "
        "  -gff EVM.combined.pasa2.gff3 "
        "  -simple "
        "  -delete_ns 5 "
        "  -strip_name "
        "  -fix_first_phase "
        "  -code {params.code} && "
        # Resolve the produced name to our declared canonical output.
        "if [ -s EVM.combined.pasa2.gff3.renamed.gff3 ]; then "
        "  cp EVM.combined.pasa2.gff3.renamed.gff3 EVM.combined.pasa2.renamed.gff3; "
        "else "
        "  for f in EVM.combined.pasa2.gff3*renamed*.gff3; do "
        "    [ -s \"$f\" ] && cp \"$f\" EVM.combined.pasa2.renamed.gff3 && break; "
        "  done; "
        "fi && "
        "test -s EVM.combined.pasa2.renamed.gff3 || "
        "    {{ echo 'FATAL: create_features_from_gff3.pl produced no renamed GFF' >&2; exit 1; }}"


# ---------------------------------------------------------------------------
# §3j.4  ogs_fix_validate
# fix_pasa_evm_gff.pl normalises a few PASA-EVM-specific glitches; gt
# gff3 -sort -tidy is the GenomeTools validator/canonicaliser.
# ---------------------------------------------------------------------------
rule ogs_fix_validate:
    input:
        renamed = rules.ogs_create_features.output.renamed,
    output:
        sorted = f"{_OGS_DIR}/EVM.combined.pasa2.renamed.sorted.gff3",
    params:
        ogs_dir  = _OGS_DIR,
        repo_bin = _REPO_BIN_OGS,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.ogs_dir} && "
        # fix_pasa_evm_gff.pl writes <input>.undefined.gff3 (per its v1
        # contract). gt then sorts + tidies + retains IDs.
        "fix_pasa_evm_gff.pl EVM.combined.pasa2.renamed.gff3 && "
        "test -s EVM.combined.pasa2.renamed.gff3.undefined.gff3 || "
        "    {{ echo 'FATAL: fix_pasa_evm_gff.pl produced no .undefined.gff3' >&2; exit 1; }} && "
        "gt gff3 -sort -tidy -retainids -force "
        "    -o EVM.combined.pasa2.renamed.sorted.gff3 "
        "    EVM.combined.pasa2.renamed.gff3.undefined.gff3"


# ---------------------------------------------------------------------------
# §3j.5  ogs_emit
# Final OGS products via gffread. Per plan: canonical filenames are
# OGS.{gff3,mRNA.fasta,CDS.fasta,pep.fasta,gtf,bed} under $OUTDIR/.
# The <code> prefix is applied internally by create_features (step 3) and
# carried through; canonical names at the workflow contract level do NOT
# include it. Symlinks let `rule all` find them at $OUTDIR/OGS.gff3.
# ---------------------------------------------------------------------------
rule ogs_emit:
    input:
        sorted = rules.ogs_fix_validate.output.sorted,
    output:
        gff   = f"{OUTDIR}/OGS.gff3",
        mrna  = f"{OUTDIR}/OGS.mRNA.fasta",
        cds   = f"{OUTDIR}/OGS.CDS.fasta",
        pep   = f"{OUTDIR}/OGS.pep.fasta",
        gtf   = f"{OUTDIR}/OGS.gtf",
        bed   = f"{OUTDIR}/OGS.bed",
    params:
        outdir     = OUTDIR,
        genome_abs = _GENOME_ABS_OGS,
        sorted_abs = lambda _wc, input: _os.path.abspath(input.sorted),
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 4000,
    shell:
        "cd {params.outdir} && "
        "cp {params.sorted_abs} OGS.gff3 && "
        # gffread one-shot: -y peptide, -x CDS, -w mRNA, -T GTF, --bed BED.
        "gffread {params.sorted_abs} -g {params.genome_abs} "
        "  -y OGS.pep.fasta "
        "  -x OGS.CDS.fasta "
        "  -w OGS.mRNA.fasta && "
        "gffread {params.sorted_abs} -T -o OGS.gtf && "
        "gffread {params.sorted_abs} --bed -o OGS.bed"


# ---------------------------------------------------------------------------
# Optional manual-curations sub-rule (plan §3j optional). When the config
# sets `manual_curations: <path>`, the GFF is tag-and-sorted and dropped
# into the EVM workdir to be picked up by evm_abinitio_cat via extra_gff.
# This rule is informational only — actual wiring through evm_abinitio_cat
# happens via config["evm"]["extra_gff"].
# ---------------------------------------------------------------------------
if _MANUAL_CUR:
    rule merge_manual_curations:
        input:
            manual = _MANUAL_CUR,
        output:
            tagged = f"{_OGS_DIR}/manual_curations.gff3.out.sorted",
        params:
            ogs_dir  = _OGS_DIR,
            repo_bin = _REPO_BIN_OGS,
        container: "containers/jamg.sif"
        threads: 1
        resources:
            mem_mb = 2000,
        shell:
            "export PATH={params.repo_bin}:$PATH && "
            "mkdir -p {params.ogs_dir} && cd {params.ogs_dir} && "
            "cp {input.manual} manual_curations.gff3 && "
            "add_source_gff.pl manual_curations.gff3 CURATIONS && "
            "sort_gff3.pl manual_curations.gff3.out"
