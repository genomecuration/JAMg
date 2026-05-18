# JAMg containers

This directory holds the recipe and supporting files for the JAMg v2 Apptainer
SIF image (`jamg.sif`).

## Files

| File | Purpose |
|---|---|
| `jamg.def` | Apptainer definition file (primary build recipe) |
| `Dockerfile` | Docker equivalent of `jamg.def` (same base image, 1:1 translation) |
| `environment.yml` | Conda/mamba environment spec (rendered from `pixi.toml`; do not edit by hand) |
| `post-install.sh` | Builds exonerate, cdbtools, and ParaFly from source; configures RepeatMasker |
| `source-pins.toml` | Exact commit SHAs for source-compiled tools |
| `trf-license.txt` | TRF license acceptance record (required for RepeatMasker configure) |

## Building the SIF

From the repo root:

```bash
# Render environment.yml from pixi.toml (after pixi install -e sif):
pixi workspace export conda-environment -e sif containers/environment.yml

# Build the SIF (requires Apptainer 1.x):
make sif
```

`make sif` automatically adds `--fakeroot` when the host supports it.
If fakeroot is unavailable, the build falls back to a plain build (may need root).

## Verifying the SIF

```bash
bash tests/container/test_tools_present.sh containers/jamg.sif
```

## Updating the tool stack

1. Edit `[feature.sif.dependencies]` in `pixi.toml`.
2. Run `pixi workspace export conda-environment -e sif containers/environment.yml`.
3. Run `make sif` to rebuild.
4. Run `bash tests/container/test_tools_present.sh` to verify.

For tools compiled from source (exonerate, cdbfasta, ParaFly), update the SHA
in `source-pins.toml` and the matching variable in `post-install.sh`.
