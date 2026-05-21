# Trinity genome-guided assembly. Consumes the sorted RNA-seq BAM from
# rnaseq.smk and produces Trinity-GG.fasta. Runs inside trinity.sif
# (upstream Docker trinityrnaseq/trinityrnaseq:2.15.2) because Trinity is
# NOT in jamg.sif.

import os as _os


rule tgg_trinity:
    input:
        bam = rules.rnaseq_sort_bam.output.bam,
        bai = rules.rnaseq_sort_bam.output.bai,
    output:
        fasta = f"{OUTDIR}/tgg/Trinity-GG.fasta",
    params:
        max_intron = config["max_intron"],
        workdir    = f"{OUTDIR}/tgg/trinity_out_dir",
        bam_abs    = lambda _wc, input: _os.path.abspath(input.bam),
        fasta_abs  = lambda _wc, output: _os.path.abspath(output.fasta),
    container: "containers/trinity.sif"
    threads: THREADS
    resources:
        mem_mb = 32000,
    shell:
        "mkdir -p $(dirname {output.fasta}) && "
        "rm -rf {params.workdir} && "
        # Trinity --max_memory wants Gb (e.g. "8G"); snakemake's
        # resources.mem_mb is in megabytes. Round-down to whole Gb,
        # minimum 1G.
        "mem_gb=$(({resources.mem_mb} / 1024)); [ $mem_gb -lt 1 ] && mem_gb=1; "
        "Trinity --genome_guided_bam {params.bam_abs} "
        "        --genome_guided_max_intron {params.max_intron} "
        "        --genome_guided_min_coverage 2 "
        "        --max_memory ${{mem_gb}}G "
        "        --CPU {threads} "
        "        --full_cleanup "
        "        --output {params.workdir} && "
        # --full_cleanup renames {workdir}/Trinity-GG.fasta to
        # {workdir}.Trinity-GG.fasta (dot suffix, sibling of workdir)
        # and removes {workdir}/. Source the renamed path.
        "cp {params.workdir}.Trinity-GG.fasta {params.fasta_abs}"
