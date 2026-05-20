# EvidenceModeler (EVM) integration. Plan §3i.
#
# Pipeline shape:
#   1-7  per-source tag + sort the GFFs that feed EVM
#   8    cat the ab-initio predictions into a single sorted file
#   9    preflight: assert every source token has a weights row
#   10   stage repeats + symlink everything into the EVM workdir
#   11   run EVidenceModeler
#
# EVM expects ABINITIO_PREDICTION GFFs from {AUGUSTUS, GeneMarkHMM, GOLDEN,
# PASA}; TRANSCRIPT alignments from PASA; PROTEIN alignments from SWISSPROT.
# Optional `evm.extra_gff[]` items (e.g. liftoff projections, manual
# curations) flow through evm_extra_gff into the same cat step.

import os as _os

_EVM_DIR        = f"{OUTDIR}/evm"
_GOLDEN_DIR_OUT = f"{OUTDIR}/golden"
_GM_DIR_OUT     = f"{OUTDIR}/genemark"
_PASA_DIR_OUT   = f"{OUTDIR}/pasa"
_AUG_DIR_OUT    = f"{OUTDIR}/augustus"
_PROT_DIR_OUT   = f"{OUTDIR}/proteins"
_REP_DIR_OUT    = f"{OUTDIR}/repeats"

_REPO_BIN       = _os.path.abspath(_os.path.join(workflow.basedir, "..", "bin"))
_GENOME_ABS_EVM = _os.path.abspath(config["genome"])
_EVM_WEIGHTS    = _os.path.abspath(config.get("evm", {}).get("weights_file") or
                                   _os.path.join(workflow.basedir, "config", "evm_weights.example.txt"))
_EXTRA_GFF      = config.get("evm", {}).get("extra_gff") or []


# ---------------------------------------------------------------------------
# §3i.1  evm_tag_golden
# Symlinks golden_genes.gff3 → final_golden_genes.gff3.nr.golden.gff3
# (the §3e canonical name) and applies GOLDEN source-tag + ID prefix.
# Per plan: prevents ID collisions across predictor sources at EVM cat-time.
# ---------------------------------------------------------------------------
rule evm_tag_golden:
    input:
        golden = rules.golden_genes.output.gff,
    output:
        gff = f"{_GOLDEN_DIR_OUT}/golden_genes.gff3.out.sorted",
    params:
        golden_dir = _GOLDEN_DIR_OUT,
        repo_bin   = _REPO_BIN,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.golden_dir} && "
        "ln -sf final_golden_genes.gff3.nr.golden.gff3 golden_genes.gff3 && "
        "add_source_gff.pl golden_genes.gff3 GOLDEN && "
        "sort_gff3.pl golden_genes.gff3.out && "
        # Plan §3i.1 (review-#8 I-8-3): rewrite ID= and Parent= to prefix
        # `golden.`; comma-handling on Parent= is defensive (single-Parent
        # today, alt-splice tomorrow).
        "perl -i -pe 's{{ID=([^;\\s]+)}}{{q{{ID=golden.}} . $1}}ge' golden_genes.gff3.out.sorted && "
        "perl -i -pe 's{{Parent=([^;\\s]+)}}{{q{{Parent=}} . join(q{{,}}, map qq{{golden.$_}}, split /,/, $1)}}ge' golden_genes.gff3.out.sorted"


# ---------------------------------------------------------------------------
# §3i.2  evm_tag_genemark
# Genemark already carries its own source tag (GeneMarkHMM) from §3f's
# gtf_to_gff3_format.pl. Only sort.
# ---------------------------------------------------------------------------
rule evm_tag_genemark:
    input:
        gff = rules.genemark.output.gff3,
    output:
        gff = f"{_GM_DIR_OUT}/genemark.gff3.sorted",
    params:
        gm_dir   = _GM_DIR_OUT,
        repo_bin = _REPO_BIN,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.gm_dir} && "
        "sort_gff3.pl genemark.gff3"


# ---------------------------------------------------------------------------
# §3i.3  evm_tag_pasa_transdecoder
# ---------------------------------------------------------------------------
rule evm_tag_pasa_transdecoder:
    input:
        gff = rules.pasa_compare_transdecoder.output.transdecoder,
    output:
        gff = f"{_PASA_DIR_OUT}/pasa.transdecoder.genome.gff3.out.sorted",
    params:
        pasa_dir = _PASA_DIR_OUT,
        repo_bin = _REPO_BIN,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.pasa_dir} && "
        "add_source_gff.pl pasa.transdecoder.genome.gff3 PASA && "
        "sort_gff3.pl pasa.transdecoder.genome.gff3.out"


# ---------------------------------------------------------------------------
# §3i.4  evm_tag_augustus
# ---------------------------------------------------------------------------
rule evm_tag_augustus:
    input:
        gff = rules.augustus.output.gff,
    output:
        gff = f"{_AUG_DIR_OUT}/augustus_results.gff3.out.sorted",
    params:
        aug_dir  = _AUG_DIR_OUT,
        repo_bin = _REPO_BIN,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.aug_dir} && "
        "add_source_gff.pl augustus_results.gff3 AUGUSTUS && "
        "sort_gff3.pl augustus_results.gff3.out"


# ---------------------------------------------------------------------------
# §3i.5  evm_tag_proteins
# Canonicalises the output filename to proteins.gff3 (EVM's expected name).
# ---------------------------------------------------------------------------
rule evm_tag_proteins:
    input:
        gff = rules.proteins_to_gff.output.gff,
    output:
        gff = f"{_PROT_DIR_OUT}/proteins.gff3",
    params:
        prot_dir = _PROT_DIR_OUT,
        repo_bin = _REPO_BIN,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.prot_dir} && "
        "add_source_gff.pl swissprot.blastx.gff3 SWISSPROT && "
        "mv -f swissprot.blastx.gff3.out proteins.gff3"


# ---------------------------------------------------------------------------
# §3i.6  evm_tag_transcript_alignments
# PASA assemblies → transcript_alignments.gff3.out. Note the deliberate
# PASA source-token collision with §3i.3 (EVM weights file carries both
# rows: ABINITIO PASA and TRANSCRIPT PASA).
# ---------------------------------------------------------------------------
rule evm_tag_transcript_alignments:
    input:
        gff = rules.pasa_compare_transdecoder.output.pasa_assemblies,
    output:
        gff = f"{_PASA_DIR_OUT}/transcript_alignments.gff3.out",
    params:
        pasa_dir = _PASA_DIR_OUT,
        repo_bin = _REPO_BIN,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "cd {params.pasa_dir} && "
        "ln -sf pasa_assemblies.gff3 transcript_alignments.gff3 && "
        "add_source_gff.pl transcript_alignments.gff3 PASA && "
        "perl -i -pe 's{{Parent=([^;\\s]+)}}{{q{{Parent=}} . join(q{{,}}, split /,/, $1)}}ge' transcript_alignments.gff3.out"


# ---------------------------------------------------------------------------
# §3i.7  evm_extra_gff  (optional)
# Applies the tag-and-sort idiom to each user-supplied extra GFF. Outputs
# land beside the source path (caller controls placement).
# ---------------------------------------------------------------------------
def _extra_gff_outputs():
    return [f"{p['path']}.out.sorted" for p in _EXTRA_GFF]


if _EXTRA_GFF:
    rule evm_extra_gff:
        input:
            sources = [p["path"] for p in _EXTRA_GFF],
        output:
            gffs = _extra_gff_outputs(),
        params:
            entries  = _EXTRA_GFF,
            repo_bin = _REPO_BIN,
        container: "containers/jamg.sif"
        threads: 1
        resources:
            mem_mb = 2000,
        run:
            for p in params.entries:
                src = p["path"]
                tag = p["source_tag"]
                shell(
                    "export PATH={repo_bin}:$PATH && "
                    "add_source_gff.pl {src} {tag} && "
                    "sort_gff3.pl {src}.out".format(repo_bin=params.repo_bin, src=src, tag=tag)
                )


# ---------------------------------------------------------------------------
# §3i.8  evm_abinitio_cat
# Concatenates the 4 ab-initio (predictor) sources + any extras into one
# sorted abinitio_gene_predictions.gff3.sorted file, the EVM input.
# ---------------------------------------------------------------------------
rule evm_abinitio_cat:
    input:
        genemark      = rules.evm_tag_genemark.output.gff,
        augustus      = rules.evm_tag_augustus.output.gff,
        golden        = rules.evm_tag_golden.output.gff,
        pasa_td       = rules.evm_tag_pasa_transdecoder.output.gff,
        extra         = _extra_gff_outputs() if _EXTRA_GFF else [],
    output:
        gff = f"{_EVM_DIR}/abinitio_gene_predictions.gff3.sorted",
    params:
        evm_dir       = _EVM_DIR,
        repo_bin      = _REPO_BIN,
        genemark_abs  = lambda _wc, input: _os.path.abspath(input.genemark),
        augustus_abs  = lambda _wc, input: _os.path.abspath(input.augustus),
        golden_abs    = lambda _wc, input: _os.path.abspath(input.golden),
        pasa_td_abs   = lambda _wc, input: _os.path.abspath(input.pasa_td),
        extra_abs     = lambda _wc, input: " ".join(_os.path.abspath(p) for p in input.extra),
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "export PATH={params.repo_bin}:$PATH && "
        "mkdir -p {params.evm_dir} && cd {params.evm_dir} && "
        "cat {params.genemark_abs} {params.augustus_abs} {params.golden_abs} {params.pasa_td_abs} "
        "    {params.extra_abs} > abinitio_gene_predictions.gff3 && "
        "sort_gff3.pl abinitio_gene_predictions.gff3"


# ---------------------------------------------------------------------------
# §3i.9  evm_preflight
# Subset check: every source token in the cat'd ab-initio GFF must have a
# matching ABINITIO_PREDICTION row in evm_weights.txt. Extra weights rows
# are fine. Also asserts tab-delimited rows in the GFF.
# ---------------------------------------------------------------------------
rule evm_preflight:
    input:
        gff     = rules.evm_abinitio_cat.output.gff,
        weights = _EVM_WEIGHTS,
    output:
        sentinel = f"{_EVM_DIR}/.preflight.ok",
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 1000,
    shell:
        # Tab-delimiter precheck (review-#7 I-7-5). Skip comments and
        # whitespace-only lines (GFF3 allows blank record-separators that
        # may carry stray spaces under some emitters).
        "awk -F'\\t' '!/^#/ && /[^[:space:]]/ {{ if (NF < 9) {{ "
        "    printf(\"FATAL: non-tab-delimited row in %s:%d (got %d fields)\\n\", FILENAME, NR, NF); exit 1 "
        "}} }}' {input.gff} && "
        # Subset check (review-#7 C-7-2): predictions' source tokens ⊆ weights' ABINITIO_PREDICTION sources.
        "unknown=$(comm -23 "
        "    <(awk -F'\\t' '!/^#/ && NF>=2 {{print $2}}' {input.gff} | sort -u) "
        "    <(awk '$1 ~ /^ABINITIO/ {{print $2}}' {input.weights} | sort -u)) && "
        "if [ -n \"$unknown\" ]; then "
        "    echo 'FATAL: predictions GFF contains source tokens absent from EVM weights:' >&2 ; "
        "    echo \"$unknown\" >&2 ; "
        "    exit 1 ; "
        "fi && touch {output.sentinel}"


# ---------------------------------------------------------------------------
# §3i.10  evm_stage_repeats
# Symlink everything EVM needs into {OUTDIR}/evm/ + strip comments from
# the merged repeats GFF (EVM doesn't tolerate header lines).
# ---------------------------------------------------------------------------
rule evm_stage_repeats:
    input:
        preflight    = rules.evm_preflight.output.sentinel,
        repeats      = rules.repeats_merge.output.merged,
        golden       = rules.evm_tag_golden.output.gff,
        pasa_td      = rules.evm_tag_pasa_transdecoder.output.gff,
        genemark     = rules.evm_tag_genemark.output.gff,
        augustus     = rules.evm_tag_augustus.output.gff,
        proteins     = rules.evm_tag_proteins.output.gff,
        tr_align     = rules.evm_tag_transcript_alignments.output.gff,
        abinitio_cat = rules.evm_abinitio_cat.output.gff,
    output:
        repeats_local = f"{_EVM_DIR}/all_repeat_masks.gff3",
    params:
        evm_dir          = _EVM_DIR,
        repeats_abs      = lambda _wc, input: _os.path.abspath(input.repeats),
        golden_abs       = lambda _wc, input: _os.path.abspath(input.golden),
        pasa_td_abs      = lambda _wc, input: _os.path.abspath(input.pasa_td),
        genemark_abs     = lambda _wc, input: _os.path.abspath(input.genemark),
        augustus_abs     = lambda _wc, input: _os.path.abspath(input.augustus),
        proteins_abs     = lambda _wc, input: _os.path.abspath(input.proteins),
        tr_align_abs     = lambda _wc, input: _os.path.abspath(input.tr_align),
        abinitio_cat_abs = lambda _wc, input: _os.path.abspath(input.abinitio_cat),
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 1000,
    shell:
        "cd {params.evm_dir} && "
        "grep -v '^#' {params.repeats_abs} > all_repeat_masks.gff3 && "
        "ln -sf {params.golden_abs}       . && "
        "ln -sf {params.pasa_td_abs}      . && "
        "ln -sf {params.genemark_abs}     . && "
        "ln -sf {params.augustus_abs}     . && "
        "ln -sf {params.proteins_abs}     . && "
        "ln -sf {params.tr_align_abs}     ."
        # abinitio_gene_predictions.gff3.sorted is already in {evm_dir}
        # (evm_abinitio_cat writes there directly); no symlink needed.


# ---------------------------------------------------------------------------
# §3i.11  evm_run
# Invokes EVidenceModeler. Output is consumed by §3j (ogs.smk).
# ---------------------------------------------------------------------------
rule evm_run:
    input:
        sentinel      = rules.evm_stage_repeats.output.repeats_local,
        weights       = _EVM_WEIGHTS,
        abinitio_cat  = rules.evm_abinitio_cat.output.gff,
        proteins      = rules.evm_tag_proteins.output.gff,
        tr_align      = rules.evm_tag_transcript_alignments.output.gff,
    output:
        evm_gff = f"{_EVM_DIR}/EVM.gff3",
    params:
        evm_dir     = _EVM_DIR,
        genome_abs  = _GENOME_ABS_EVM,
        weights_abs = _EVM_WEIGHTS,
        sample_id   = lambda _wc: _os.path.splitext(_os.path.basename(_GENOME_ABS_EVM))[0],
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 32000,
    shell:
        # /opt/jamg/bin/EVidenceModeler is a symlink to
        # /opt/jamg/share/EVidenceModeler/EVidenceModeler, but EVM uses
        # $FindBin::Bin to locate its sibling EvmUtils/ + PerlLib/, and
        # FindBin returns the symlink's directory, not the target. Invoke
        # by the real path so EvmUtils/ resolves alongside it.
        "cd {params.evm_dir} && "
        "/opt/jamg/share/EVidenceModeler/EVidenceModeler "
        "  --weights {params.weights_abs} "
        "  --segmentSize 5000000 --overlapSize 80000 "
        "  --sample_id {params.sample_id} "
        "  --genome {params.genome_abs} "
        "  --gene_predictions abinitio_gene_predictions.gff3.sorted "
        "  --protein_alignments proteins.gff3 "
        "  --transcript_alignments transcript_alignments.gff3.out "
        "  --repeats all_repeat_masks.gff3 "
        "  --CPU {threads} && "
        # EVM names its output <sample_id>.EVM.gff3; canonicalise to EVM.gff3.
        "cp {params.sample_id}.EVM.gff3 EVM.gff3"
