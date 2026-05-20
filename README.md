# JAMg v2 — Just Annotate My Genome

Eukaryotic genome-annotation meta-pipeline. Snakemake DAG wrapping
Augustus, GeneMark-ES, PASA (TransDecoder), RepeatMasker / RepeatModeler,
Trinity-GG, EvidenceModeler, and a small in-tree Perl glue layer (golden
gene-set construction, hint emission, Augustus chunking).

## Quick start

```bash
# 1. Build the four Apptainer images. RM_LIB_HOST is your RepeatMasker
#    Libraries directory (host-supplied; license-encumbered).
make sifs RM_LIB_HOST=/path/to/RepeatMasker/Libraries

# 2. Author a per-genome config YAML. test_suite/mini-config.yaml is
#    a working example against the bundled D. melanogaster X-fragment
#    fixtures.
cp test_suite/mini-config.yaml my-genome.yaml
$EDITOR my-genome.yaml

# 3. Run the DAG. Produces $OUTDIR/OGS.gff3 plus the matching FASTA /
#    GTF / BED products.
bin/jamg run --config my-genome.yaml --cores 20
```

## Outputs (under `outdir:` from the config)

| File | What |
|---|---|
| `OGS.gff3` | Official Gene Set, GFF3 |
| `OGS.mRNA.fasta` | mRNA sequences |
| `OGS.CDS.fasta`  | CDS sequences  |
| `OGS.pep.fasta`  | Protein sequences |
| `OGS.gtf` | GTF view of OGS.gff3 |
| `OGS.bed` | BED view of OGS.gff3 |

## Repository layout

| Path | Purpose |
|---|---|
| `bin/` | Perl + Python helpers. `bin/jamg` is the CLI entry point. |
| `PerlLib/` | Perl modules (`Golden::Filter`, `Golden::Alignment`, `Golden::Augustus`; GFF / FASTA / BioPerl glue). |
| `containers/` | Apptainer build recipes + (gitignored) built `.sif` files. |
| `workflow/` | Snakemake DAG (`Snakefile`, `rules/*.smk`, `config/`). |
| `tests/` | Unit, integration, e2e, and policy tests. |
| `test_suite/` | Mini-fixtures (1 Mb dmel-X fragment + RNA-seq + transcripts) for smoke. |
| `tools/` | Build and fixture-construction helpers. |
| `docs/` | Reference docs: see below. |

## Documentation

- [`docs/procedure_v2.asciidoc`](docs/procedure_v2.asciidoc) — pipeline walkthrough.
- [`docs/containers.md`](docs/containers.md) — Apptainer images: build, rebuild, version pinning.
- [`docs/migration-v1-to-v2.md`](docs/migration-v1-to-v2.md) — what changed since v1.x.

## License

CSIRO 2012-, see LICENSE.
