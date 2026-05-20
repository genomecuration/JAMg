# Containers

JAMg v2 ships four Apptainer images that together carry every binary
dependency. The host needs Apptainer (or Singularity ≥ 3.8) and a writable
build cache.

## Image inventory

| Image | Base | Carries | Built via |
|---|---|---|---|
| `containers/jamg-base.sif` | `ubuntu:24.04` + bioconda | bioconda toolchain, system Perl + apt-installed BioPerl, Click + Snakemake | `make sif-base` |
| `containers/jamg.sif` | `jamg-base.sif` (`Bootstrap: localimage`) | Augustus master, EVidenceModeler (alpapan fork), TransDecoder, GeneMark glue, RepeatMasker + Libraries (vendored via git-lfs at `containers/rm_libs/`), JAMg `bin/` + `PerlLib/` | `make sif` |
| `containers/trinity.sif` | upstream Docker | Trinity (de novo + GG modes) | `make sif-trinity` |
| `containers/pasa.sif` | upstream Docker | PASA pipeline, Launch_PASA_pipeline.pl, pasa_asmbls_to_training_set.dbi | `make sif-pasa` |

`trinity.sif` and `pasa.sif` are pulled unmodified from upstream and have
no `%post` step; their glibc 2.31 conflicts with apptainer's libfakeroot
GLIBC_2.33+ requirement, so any in-build modification fails at faked
startup.

## Building

```bash
# RepeatMasker libraries (license-encumbered) are vendored under
# containers/rm_libs/ via git-lfs. A clean clone followed by 'git lfs pull'
# is sufficient; jamg.def %files mounts the directory into the SIF build at
# /rm_lib_host and the in-SIF post-install.sh copies it into place.
make sifs
```

Individual targets: `make sif-base`, `make sif-trinity`, `make sif-pasa`,
then `make sif`. The four images are independent except
for `jamg.sif` layering on `jamg-base.sif`.

## Rebuilding for a tool upgrade

1. Edit `containers/environment.yml` (bioconda packages) or
   `containers/post-install.sh` (in-place pip/cpan installs).
2. `rm -f containers/jamg-base.sif containers/jamg.sif`
3. `make sifs`

Build logs land in `containers/build-base.log` and `containers/build.log`.

## Version pinning

Bioconda versions are pinned in `containers/environment.yml`. Augustus,
EVidenceModeler, TransDecoder, and GeneMark glue install at fixed
upstream commits enumerated in `containers/jamg.def`'s `%post`.

## Running without Apptainer

`containers/Dockerfile` provides a Docker fallback for the `jamg.sif`
layer. The PASA and Trinity images come from their upstream Docker tags
directly. The bundled `bin/jamg` defaults to apptainer; pass
`--container-engine docker` to switch.

## PATH conventions inside the SIF

`jamg.def`'s `%environment` puts `/opt/jamg/bin:/usr/bin:/opt/conda/bin:...`
on PATH so that `#!/usr/bin/env perl` and `env python3` shebangs resolve
to system Perl + system Python (the apt versions carry BioPerl + Click +
PyYAML + jsonschema). Bioconda tools that exist only at `/opt/conda/bin`
still resolve via fallthrough.

Host repo's `bin/` is bind-mounted into the SIF at the repo root.
Snakemake rules that need a host script add `$REPO/bin` to PATH inline.
