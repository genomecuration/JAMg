# Protein evidence: blastx softmasked genome against the SwissProt subset
# -> blast2gff.py -> gff2hints.pl. Consumed downstream by EVM (3i) and the
# Augustus hint pile (3h).

import os as _os

_PROT_DIR     = f"{OUTDIR}/proteins"
_GENOME_ABS   = _os.path.abspath(config["genome"])
_SWISSPROT_DB = _os.path.abspath(config["proteins"]["swissprot_db"])


rule proteins_blastx:
    # blastx softmasked genome against the supplied SwissProt subset.
    # The softmasked genome from repeats_merge is used so blastx hits to
    # repetitive regions are suppressed.
    input:
        softmasked     = rules.repeats_merge.output.soft,
        softmasked_fai = rules.repeats_merge.output.soft_fai,
    output:
        tsv = f"{_PROT_DIR}/swissprot.blastx.tsv",
    params:
        db_abs     = _SWISSPROT_DB,
        max_intron = config["max_intron"],
    container: "containers/jamg.sif"
    threads: THREADS
    resources:
        mem_mb = 8000,
    shell:
        "mkdir -p $(dirname {output.tsv}) && "
        "blastx -query {input.softmasked} -db {params.db_abs} "
        "       -outfmt 6 -num_threads {threads} -evalue 1e-10 "
        "       -max_target_seqs 5 -max_intron_length {params.max_intron} -lcase_masking "
        "       -out {output.tsv}"


rule proteins_to_gff:
    # Convert blastx tab output to GFF3 via blast2gff.py. The script hardcodes
    # source=BLAST; the SWISSPROT source tag is applied later in §3i's
    # evm_tag_proteins step via add_source_gff.pl.
    input:
        tsv = rules.proteins_blastx.output.tsv,
    output:
        gff = f"{_PROT_DIR}/swissprot.blastx.gff3",
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "blast2gff.py -i {input.tsv} -o {output.gff}"


rule proteins_hints:
    # Emit augustus hints from the blastx GFF. gff2hints.pl writes to
    # <input>.hints which already matches the declared output.
    input:
        gff = rules.proteins_to_gff.output.gff,
    output:
        hints = f"{_PROT_DIR}/swissprot.blastx.gff3.hints",
    container: "containers/jamg.sif"
    threads: 1
    resources:
        mem_mb = 2000,
    shell:
        "gff2hints.pl {input.gff}"
