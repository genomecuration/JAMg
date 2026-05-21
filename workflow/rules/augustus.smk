# Augustus gene prediction with chunked ParaFly execution. Hint pile sources
# (per plan §3h, NOT including GeneMark — its output goes straight to EVM):
#   - repeats merged hints (RM)
#   - rnaseq hints (JR + RCOV)
#   - PASA assembly hints (PASA) + polyA hints
#   - golden hints (GLD)
#   - protein hints (HU from blastx)
#
# augustus_preflight (§3h.0) gates the main rule with 3 fail-fast checks
# so a misconfigured species or extrinsic.cfg surfaces in seconds rather
# than after the chunked ParaFly run.

import os as _os

_AUG_DIR        = f"{OUTDIR}/augustus"
_GENOME_ABS     = _os.path.abspath(config["genome"])
_EXTRINSIC_ABS  = _os.path.abspath(config["augustus"]["extrinsic_cfg"])


rule augustus_preflight:
    # Surfaces three misconfigurations before any compute is spent:
    #   1. species not in Augustus's library;
    #   2. extrinsic.cfg [SOURCES] line drifted from the v2 token set;
    #   3. extrinsic.cfg fails to parse under Augustus 3.5.
    # Plan §3h.0 incorporated reviewer findings #7 (C-7-1 species=help,
    # I-7-3 awk-after-[SOURCES], I-7-1 err-keyword grep, I-7-5).
    output:
        sentinel = touch(f"{_AUG_DIR}/.preflight.ok"),
    params:
        species   = config["species"],
        extrinsic = _EXTRINSIC_ABS,
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 1000,
    shell:
        # (1) species exists in Augustus's species library. Check the species
        # directory directly: augustus reads $AUGUSTUS_CONFIG_PATH/species/<name>/
        # at load time, so the directory's existence is the authoritative test.
        # Previously `augustus --species=help | grep -q ...` was used, but the
        # combination of `set -euo pipefail` + grep -q + augustus's 118-line
        # output triggers a SIGPIPE race: grep -q exits on first match, the
        # pipe closes, augustus exits with SIGPIPE status, pipefail propagates
        # the non-zero exit even though grep matched. The dir check is immune.
        "if [ ! -d \"$AUGUSTUS_CONFIG_PATH/species/{params.species}\" ]; then "
        "    echo \"FATAL: species '{params.species}' not in Augustus library.\"; "
        "    echo \"Looked at: $AUGUSTUS_CONFIG_PATH/species/{params.species}\"; "
        "    echo 'Available species:'; "
        "    ls -1 \"$AUGUSTUS_CONFIG_PATH/species/\" | head -40; "
        "    echo 'Either pick a built-in species OR set augustus.optimise: true to train one.'; "
        "    exit 1; "
        "fi && "
        # (2) extrinsic.cfg [SOURCES] line matches the v2 token set.
        # HU stays (Augustus protein-hint source, used by gff2hints.pl on
        # BLASTX output); PASA is NOT in this line (PASA is a transcript-
        # alignments source, not an Augustus hints source).
        "awk '/^\\[SOURCES\\]/ {{found=1; next}} found {{print; exit}}' {params.extrinsic} "
        "  | grep -qE '^M JR RCOV XNT RM HU GLD E\\s*$' || "
        "  {{ echo 'FATAL: extrinsic.cfg [SOURCES] line does not match v2 token set.'; "
        "     echo 'Expected: M JR RCOV XNT RM HU GLD E'; exit 1; }} && "
        # (3) extrinsic.cfg parses cleanly under Augustus 3.5. Capture stderr
        # only (stdout → /dev/null), `|| true` swallows the non-zero exit
        # from feeding augustus an empty /dev/null input. Then grep stderr
        # for known parse-error keywords; if absent, the parse passed.
        "augustus_err=$(augustus --extrinsicCfgFile={params.extrinsic} "
        "                         --species={params.species} /dev/null 2>&1 1>/dev/null || true); "
        "if echo \"$augustus_err\" | grep -qE 'extrinsicCfgFile|parse error|cfg|configuration'; then "
        "    echo \"FATAL: extrinsic.cfg parse error under Augustus 3.5:\"; "
        "    echo \"$augustus_err\"; exit 1; "
        "fi"


rule augustus:
    # Concatenate all hint sources into one hintfile, partition the genome
    # + hints into per-chunk files via bin/run_split_augustus.py, execute
    # the per-chunk augustus invocations under ParaFly, merge results.
    # GeneMark output is NOT in the hint pile (goes to EVM directly).
    input:
        preflight           = rules.augustus_preflight.output.sentinel,
        softmasked          = rules.repeats_merge.output.soft,
        repeats_hints       = rules.repeats_merge.output.hints,
        rnaseq_hints        = rules.rnaseq_hints.output.hints,
        pasa_assembly_hints = rules.pasa_hints.output.assembly_hints,
        polya_hints         = rules.pasa_hints.output.polya_hints,
        golden_hints        = rules.golden_hints.output.hints,
        proteins_hints      = rules.proteins_hints.output.hints,
    output:
        gff = f"{_AUG_DIR}/augustus_results.gff3",
    params:
        species         = config["species"],
        extrinsic       = _EXTRINSIC_ABS,
        aug_dir_abs     = lambda _wc, output: _os.path.dirname(_os.path.abspath(output.gff)),
        softmasked_abs  = lambda _wc, input: _os.path.abspath(input.softmasked),
        # absolutise every hint input so the `cd {aug_dir_abs}` below
        # doesn't break the relative paths snakemake hands us.
        repeats_hints_abs       = lambda _wc, input: _os.path.abspath(input.repeats_hints),
        rnaseq_hints_abs        = lambda _wc, input: _os.path.abspath(input.rnaseq_hints),
        pasa_assembly_hints_abs = lambda _wc, input: _os.path.abspath(input.pasa_assembly_hints),
        polya_hints_abs         = lambda _wc, input: _os.path.abspath(input.polya_hints),
        golden_hints_abs        = lambda _wc, input: _os.path.abspath(input.golden_hints),
        proteins_hints_abs      = lambda _wc, input: _os.path.abspath(input.proteins_hints),
        repo_bin                = _os.path.abspath(_os.path.join(workflow.basedir, "..", "bin")),
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 16000,
    shell:
        # Use the host repo's bin/run_split_augustus.py (shebang
        # /opt/conda/bin/python3 to access BioPython) instead of the
        # SIF-baked copy that points at /usr/bin/python3 (no Bio).
        "export PATH={params.repo_bin}:$PATH && "
        "mkdir -p {params.aug_dir_abs} && cd {params.aug_dir_abs} && "
        # Pile all hint files into one. Per plan: repeats + rnaseq + pasa
        # assemblies + polyA + golden + proteins. NOT GeneMark.
        "cat {params.repeats_hints_abs} {params.rnaseq_hints_abs} "
        "    {params.pasa_assembly_hints_abs} {params.polya_hints_abs} "
        "    {params.golden_hints_abs} {params.proteins_hints_abs} "
        "    > genome_augustus.all.hintfile && "
        # run_split_augustus.py partitions the softmasked genome + the
        # hintfile into {threads} balanced chunks and emits per-chunk
        # augustus invocations into augustus.cmds.
        "run_split_augustus.py run "
        "    -f {params.softmasked_abs} "
        "    -s {params.species} "
        "    -n {threads} "
        "    -t genome_augustus.all.hintfile "
        "    -c {params.extrinsic} "
        "    -u on "
        "    -o augustus.cmds && "
        # ParaFly executes the chunked augustus commands in parallel.
        # Stays here (per plan, review-#2 I12: ParaFly NOT replaced by
        # snakemake parallelism — different abstraction levels).
        "ParaFly -c augustus.cmds -CPU {threads} -v -shuffle && "
        # Merge per-chunk GFF3 outputs. uniqueGeneId=true (run_split_
        # augustus.py default) means chunk gene IDs do not collide.
        # Write basename (not {output.gff}) because CWD is already
        # {params.aug_dir_abs} — same invariant as genemark.smk:118.
        "cat run/result.* > augustus_results.gff3"
