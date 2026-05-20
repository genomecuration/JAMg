# GeneMark-ES gene prediction. Two rules:
#   - genemark_preflight: verifies the licence and tarball, extracts the
#     archive, ungzips the key, ABI-checks gmes_petap.pl under SIF Perl.
#   - genemark: concatenates RNA-seq intron evidence + golden hint introns,
#     runs gmes_petap.pl --ET, converts the GTF to canonical GFF3.
#
# GeneMark is bind-mounted from the host via collect_binds (it walks
# config["genemark"]["path"] + config["genemark"]["key"] and adds the
# parent dir of each archive to --apptainer-args). HOME is overridden to
# the workdir holding .gm_key so gmes_petap.pl finds $HOME/.gm_key.

import os as _os

_GM_DIR        = f"{OUTDIR}/genemark"
_GENOME_ABS    = _os.path.abspath(config["genome"])
_GM_PATH_ABS   = _os.path.abspath(_os.path.expanduser(config["genemark"]["path"]))
_GM_KEY_ABS    = _os.path.abspath(_os.path.expanduser(config["genemark"]["key"]))


rule genemark_preflight:
    # (1) gzipped key present; (2) tarball present; (3) extract tarball
    # (idempotent); (4) ungzip key into a workdir-private $HOME; (5) ABI
    # check via gmes_petap.pl --help.
    output:
        sentinel = touch(f"{_GM_DIR}/.preflight.ok"),
    params:
        gm_tarball = _GM_PATH_ABS,
        gm_key_gz  = _GM_KEY_ABS,
        workdir    = f"{_GM_DIR}/extracted",
        gm_home    = f"{_GM_DIR}/home",
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 1000,
    shell:
        '[[ -s "{params.gm_key_gz}" ]] || '
        '  {{ echo "FATAL: GeneMark gzipped key not found at {params.gm_key_gz}"; '
        '     echo "Set config.genemark.key to the .gz path, or download gm_key_64.gz from"; '
        '     echo "  http://exon.gatech.edu/GeneMark/license_download.cgi"; '
        '     exit 1; }} && '
        '[[ -s "{params.gm_tarball}" ]] || '
        '  {{ echo "FATAL: GeneMark tarball not found at {params.gm_tarball}"; '
        '     echo "Set config.genemark.path to gmes_linux_64_4.tar.gz, or download from"; '
        '     echo "  http://exon.gatech.edu/GeneMark/license_download.cgi"; '
        '     exit 1; }} && '
        'mkdir -p "{params.workdir}" "{params.gm_home}" && '
        'if [[ ! -x "{params.workdir}/gmes_linux_64_4/gmes_petap.pl" ]]; then '
        '    tar -xzf "{params.gm_tarball}" -C "{params.workdir}"; '
        'fi && '
        '[[ -x "{params.workdir}/gmes_linux_64_4/gmes_petap.pl" ]] || '
        '  {{ echo "FATAL: gmes_petap.pl missing after extract"; exit 1; }} && '
        'if [[ ! -s "{params.gm_home}/.gm_key" ]]; then '
        '    gunzip -c "{params.gm_key_gz}" > "{params.gm_home}/.gm_key" && '
        '    chmod 600 "{params.gm_home}/.gm_key"; '
        'fi && '
        # ABI check: invoke gmes_petap.pl with no args (it prints "Usage:"
        # and exits non-zero) and confirm the "Usage:" line appears. The
        # subshell disables pipefail because gmes_petap.pl's exit=1 would
        # otherwise propagate through the pipe even when grep matches.
        # A missing Perl module would surface as "Can't locate <module>"
        # without ever reaching the usage banner.
        '( set +o pipefail; '
        '  HOME="{params.gm_home}" "{params.workdir}/gmes_linux_64_4/gmes_petap.pl" 2>&1 '
        '    | grep -q "^Usage:" ) || '
        '  {{ echo "FATAL: gmes_petap.pl could not load Perl modules under SIF Perl"; exit 1; }}'


rule genemark:
    # Run GeneMark-ET with RNA-seq intron evidence. Produces genemark.gtf
    # and a canonical GFF3 via gtf_to_gff3_format.pl.
    #
    # ancient() on inputs is INTENDED to suppress snakemake's in-session "Input
    # files updated by another job" cascade. Without it, when rnaseq_hints or
    # repeats_merge run in the same session, snakemake would re-trigger genemark
    # even though the pre-staged genemark output (test_suite/fixtures/genemark.gff3,
    # a 100 kb subset of a 1 Mb prediction) is already present. gmes_petap --ES
    # cannot converge on the 100 kb default fixture (needs >=1 Mb training data),
    # so re-running would fail. Whether ancient() fully suppresses this cascade in
    # snakemake 9.21 has been verified empirically for the test_evm and test_ogs
    # scenarios (pre-staged upstream); the test_full_dag / regen_snapshot scenario
    # (upstream rules fresh) is not yet fully verified at the time of writing.
    # Missing output still triggers the rule normally. Force a re-run after
    # upstream evidence changes with `snakemake --forcerun genemark`.
    input:
        preflight    = ancient(rules.genemark_preflight.output.sentinel),
        rnaseq_hints = ancient(rules.rnaseq_hints.output.hints),
        softmasked   = ancient(rules.repeats_merge.output.soft),
    output:
        gtf  = f"{_GM_DIR}/genemark.gtf",
        gff3 = f"{_GM_DIR}/genemark.gff3",
    params:
        gm_home_abs = lambda _wc, output: _os.path.abspath(f"{_GM_DIR}/home"),
        workdir_abs = lambda _wc, output: _os.path.abspath(f"{_GM_DIR}/extracted"),
        softmasked_abs = lambda _wc, input: _os.path.abspath(input.softmasked),
        rnaseq_abs     = lambda _wc, input: _os.path.abspath(input.rnaseq_hints),
        gm_dir_abs     = lambda _wc, output: _os.path.dirname(_os.path.abspath(output.gtf)),
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 8000,
    shell:
        "cd {params.gm_dir_abs} && "
        # Build an intron-only evidence file from the rnaseq hints. GFF3 is
        # strictly tab-delimited; use awk -F'\t' so field 3 is matched
        # exactly. POSIX ERE does not define \t and the SIF's GNU grep 3.11
        # treats it as a literal two-char sequence, so the previous
        # grep -E '...\tintron\t...' returned zero matches inside the SIF.
        "awk -F'\\t' '!/^#/ && $3==\"intron\" && $9 !~ /noncanonical=true/' "
        "    {params.rnaseq_abs} "
        "    | sort -k1,1 -k4,4n > all_evidence_introns.gff3; "
        "[ -s all_evidence_introns.gff3 ] || "
        "    {{ echo 'FATAL: no intron evidence for genemark --ET (rnaseq_hints broken or fixture lacks splice information)' >&2; exit 1; }} && "
        "HOME='{params.gm_home_abs}' "
        "{params.workdir_abs}/gmes_linux_64_4/gmes_petap.pl "
        "    --ET all_evidence_introns.gff3 --et_score 10 "
        "    --soft_mask 1 --max_mask 10000 "
        "    --cores {threads} "
        "    --sequence {params.softmasked_abs} && "
        "[ -s genemark.gtf ] || "
        "    {{ echo 'FATAL: gmes_petap.pl did not produce genemark.gtf' >&2; exit 1; }} && "
        # Convert GTF to canonical GFF3 with source=GeneMarkHMM.
        # Write the basename (not {output.gff3}) because the shell has
        # already `cd`'d into {params.gm_dir_abs}, which IS the directory
        # containing {output.gff3}.
        "gtf_to_gff3_format.pl genemark.gtf {params.softmasked_abs} GeneMarkHMM "
        "    | grep -v '^# ' > genemark.gff3"
