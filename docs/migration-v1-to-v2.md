# JAMg v1 → v2 migration

JAMg v2.0.0 is a breaking release. v1 is preserved at the `v1.x-final` tag; v2 work happens on branch `v2-revamp` until merged into `master` at v2.0.0 release.

## What is changing

| Area | v1 | v2 |
|---|---|---|
| Orchestration | Interactive shell + ~90 `bin/*.pl` orchestrators wired together by hand | `workflow/` Snakemake DAG; resumable, parallel, slurm-aware |
| Tool delivery | ~40 third-party tools vendored under `3rd_party/` and built locally | Single Apptainer/Singularity SIF (`containers/jamg.sif`) holding the bioconda tool stack + in-tree Perl glue |
| Env management | `env.sh` wizard generating per-genome `env.source` files | `pixi.toml` (host-dev) + per-genome `workflow/config/<your>.yaml` |
| Entry point | `bin/*.pl` invoked directly with sourced environment | `bin/jamg run --config <config>.yaml` (Click CLI wrapping `snakemake`) |
| Endpoint | Ad-hoc; tutorial stops at OGS but plumbing scattered | OGS GFF3 + `OGS.{mRNA,CDS,pep}.fasta` under `<outdir>/`, single canonical contract |
| Tests | None | Per-script unit tests + per-rule integration tests + end-to-end smoke on a 51 kb mini-genome |

## Predictor set reduction

| Predictor | v1 status | v2 status |
|---|---|---|
| Augustus | Live; HMM with evidence hints | Live; HMM with evidence hints |
| GeneMark-ES | Live; EVM input + Augustus hint source | Live; **EVM input only** (no longer routed to Augustus hints) |
| PASA TransDecoder | Live; ORF prediction from transcript alignments | Live; unchanged |
| SNAP | Live; bundled, called from `prepare_golden_genes_for_predictors.pl` | **Dropped.** All output emission, `gff2zff`, and the `fathom`-based validator are removed |
| geneid | Live; bundled, called from `prepare_golden_genes_for_predictors.pl` | **Dropped.** Output emission and helpers removed |
| GlimmerHMM | Live; bundled, called from `prepare_golden_genes_for_predictors.pl` | **Dropped.** `gb2glimmer` and emission removed |
| HHblits domain-exon search | `prepare_domain_exon_annotation.pl` + hh-suite + multi-GB DB stubs | **Dropped.** Orchestrator deleted; `databases/hhblits/` retired |
| GSNAP as RNA aligner | Live (`bin/align_rnaseq_gsnap.pl`) | **Retained as an option** via `config["rnaseq"]["mode"]: gsnap`; alternative to STAR or pre-aligned BAM |
| GSNAP as DNA aligner | Live (`bin/align_dnaseq_gsnap.pl`) | **Dropped** (DNA-seq alignment not in JAMg scope) |

The golden-gene structure validator that previously called SNAP's `fathom` is replaced by a pure-Perl canonical-splice + completeness check in `PerlLib/Golden/Filter.pm`.

## Container runtime move

v2 ships a single Apptainer/Singularity image (`containers/jamg.sif`) containing the bioconda tool stack and the in-tree Perl glue. The host needs only:

- Apptainer (or Singularity) on the PATH.
- `pixi` for host-side dev tooling (snakemake, click, pytest, ruff, pyright, shellcheck, perlcritic). Optional; users can also `apptainer exec jamg.sif jamg run ...`.
- A GeneMark-ES tarball + license-key archive (license-encumbered; not redistributable). Default config paths: `~/software/gmes_linux_64_4.tar.gz` and `~/software/gm_key_64.gz`. The workflow extracts both at runtime into `<outdir>/genemark/`; no manual extraction needed.

No `env.sh`, no `LD_LIBRARY_PATH`, no `PERL5LIB` to source. Per-rule shells run inside the SIF via Snakemake's `--software-deployment-method apptainer`.

## Top-level deletions in v2

Preserved at tag `v1.x-final`; gone from `v2-revamp`:

```
3rd_party/                                              entire tree
env.sh                                                  replaced by per-genome config.yaml
env.source                                              ditto
docs/procedure.asciidoc, docs/procedure.html            replaced by docs/procedure_v2.md
docs/tutorial.asciidoc,  docs/tutorial.html             replaced by README + procedure_v2
docs/creating_and_using_hhblits_databases.{asciidoc,html}
docs/index.{asciidoc,html}                              superseded by README + procedure_v2
docs/build.sh                                           obsolete asciidoc renderer
databases/hhblits/                                      HHblits DB stubs
```

`bin/` orphans deleted (full canonical list lives in `tests/no-legacy/forbidden_files.txt`):

```
bin/SNAP_output_to_gff3.pl
bin/optimize_geneid.sh
bin/purge.geneid.real.gff.pl
bin/zff2hintzff.pl
bin/align_rnaseq_gsnap.pl                superseded by rules/rnaseq.smk gsnap branch
bin/align_dnaseq_gsnap.pl                DNA-seq alignment not in scope
bin/align_rnaseq_cleanup.pl              gsnap log scraper; superseded
bin/bam_to_jbrowse.pl                    JBrowse plumbing out of scope
bin/bigwig_to_jbrowse.pl
bin/sra_to_jbrowse.pl
bin/tab_to_jbrowse.pl
bin/prepare_domain_exon_annotation.pl    HHblits orchestrator
bin/prepare_evm_calls_from_maker_output.pl
bin/download_ncbi_assembly_raw_data.pl
bin/prepare_trinity_genome_assembly_pbs.pl   superseded by rules/tgg.smk
bin/JAMG_TGG_cmds.pl                         superseded by rules/tgg.smk
bin/prepare_golden_genes_for_predictors.pl   refactored into PerlLib/Golden/ + bin/prepare_golden_genes.pl
```

## Recovering v1

```
git checkout v1.x-final           # detached HEAD at the pre-revamp tip
```

The `3rd_party/` tree, `env.sh`, and the full v1 `bin/` set are intact at that tag.

## Out of scope for v2.0.0

The following are deliberately not part of v2 and will be tracked for v2.1+ if needed:

- JBrowse / Nextcloud / Galaxy / `scpresume` plumbing
- snpEff database build + variant annotation
- Salmon decoy index build
- JCVI gene-based synteny
- JAMp functional annotation
- liftoff projection (the workflow accepts liftoff GFFs via `evm.extra_gff[]` but does not produce them)
- Manual curation merge as a workflow stage (user-supplied input only)
- Cactus pangenome / minimap2 whole-genome synteny
- Variant calling
- MariaDB-backed PASA (SQLite in `/dev/shm` is the v2 default)
- Multi-genome parallel dispatch
- ARM/aarch64 SIF builds
