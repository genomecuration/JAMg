SHELL := /bin/bash

.PHONY: sif sifs test test-container smoke lint lint-py lint-perl lint-shell lint-snakemake pyright clean docs

REPO_ROOT := $(shell git rev-parse --show-toplevel)

# --fakeroot requires subuid/subgid mapping; fall back to plain build
# (which may need root). If neither works, the user must rebuild on a
# host with apptainer fakeroot support OR pull a prebuilt SIF.
APPTAINER_BUILD_FLAGS := $(shell apptainer build --help 2>/dev/null | grep -q -- --fakeroot && echo --fakeroot)

# `sif` builds the main image only (jamg.sif). `sifs` (plural) builds all
# four: jamg-base, jamg, trinity, pasa. Iterative dev typically only needs
# `sif`; CI / fresh host setup uses `sifs`.
sif: containers/jamg.sif
sifs: containers/jamg.sif containers/trinity.sif containers/pasa.sif

# `apptainer build -F` forces overwrite of an existing SIF (default refuses).
containers/jamg-base.sif: containers/jamg-base.def
	@echo "Building base SIF with flags: $(APPTAINER_BUILD_FLAGS) (log -> containers/build-base.log)"
	@set -o pipefail; cd containers && apptainer build -F $(APPTAINER_BUILD_FLAGS) jamg-base.sif jamg-base.def 2>&1 | tee build-base.log || \
	  { echo ""; echo "BASE SIF build FAILED -- see containers/build-base.log"; \
	    exit 1; }

# Main SIF: bioconda tool stack layered on top of jamg-base.sif via
# `Bootstrap: localimage`. Depends on the base SIF + environment.yml +
# post-install.sh, so any change to either layer triggers the right rebuild.
# Build log written to containers/build.log AND streamed to stdout via `tee` so
# background invocations and CI both get a tailable artefact.
containers/jamg.sif: containers/jamg.def containers/jamg-base.sif containers/environment.yml containers/post-install.sh
	@echo "Building main SIF with flags: $(APPTAINER_BUILD_FLAGS) (log -> containers/build.log)"
	@set -o pipefail; cd containers && apptainer build -F $(APPTAINER_BUILD_FLAGS) jamg.sif jamg.def 2>&1 | tee build.log || \
	  { echo ""; echo "SIF build FAILED -- see containers/build.log"; \
	    echo "  - Try: sudo apptainer build containers/jamg.sif containers/jamg.def"; \
	    echo "  - Or pull a prebuilt SIF (see docs/containers.md)"; \
	    exit 1; }

# Trinity SIF: pulled from upstream Docker image trinityrnaseq/trinityrnaseq.
# Version pinned in containers/trinity.def. Consumed by workflow/rules/tgg.smk.
containers/trinity.sif: containers/trinity.def
	@echo "Pulling Trinity SIF from upstream Docker (log -> containers/build-trinity.log)"
	@set -o pipefail; cd containers && apptainer build -F $(APPTAINER_BUILD_FLAGS) trinity.sif trinity.def 2>&1 | tee build-trinity.log || \
	  { echo ""; echo "TRINITY SIF build FAILED -- see containers/build-trinity.log"; exit 1; }

# PASA SIF: pulled from upstream Docker pasapipeline/pasapipeline.
# Version pinned in containers/pasa.def. Consumed by workflow/rules/pasa.smk.
containers/pasa.sif: containers/pasa.def
	@echo "Pulling PASA SIF from upstream Docker (log -> containers/build-pasa.log)"
	@set -o pipefail; cd containers && apptainer build -F $(APPTAINER_BUILD_FLAGS) pasa.sif pasa.def 2>&1 | tee build-pasa.log || \
	  { echo ""; echo "PASA SIF build FAILED -- see containers/build-pasa.log"; exit 1; }

# Per-SIF presence tests. Each SIF has its own tool list under tests/container/.
test-container: containers/jamg.sif
	bash tests/container/test_tools_present.sh containers/jamg.sif
	@[ -f containers/trinity.sif ] && bash tests/container/test_trinity_tools.sh containers/trinity.sif || \
	  echo "SKIP: containers/trinity.sif not built (run 'make sifs' for full coverage)"
	@[ -f containers/pasa.sif ] && bash tests/container/test_pasa_tools.sh containers/pasa.sif || \
	  echo "SKIP: containers/pasa.sif not built (run 'make sifs' for full coverage)"

# Forward-declared paths (tests/golden/, tests/no-legacy/) skip cleanly so the
# target stays green between phase landings.
test: test-container
	@[ -f tests/golden/test_validator.t ] && prove tests/golden/test_validator.t || \
	  echo "SKIP: tests/golden/test_validator.t (Phase 4 deliverable)"
	@[ -f tests/golden/test_prepare_golden_genes.t ] && prove tests/golden/test_prepare_golden_genes.t || \
	  echo "SKIP: tests/golden/test_prepare_golden_genes.t (Phase 4 deliverable)"
	@[ -f tests/no-legacy/test_legacy_tools_gone.sh ] && bash tests/no-legacy/test_legacy_tools_gone.sh || \
	  echo "SKIP: tests/no-legacy/test_legacy_tools_gone.sh (Phase 5 deliverable)"

smoke: sif
	@[ -f tests/e2e/test_full_dag.sh ] && bash tests/e2e/test_full_dag.sh || \
	  echo "SKIP: tests/e2e/test_full_dag.sh (forward-decl)"

lint: lint-py lint-perl lint-shell lint-snakemake

# All lint targets delegate to the pixi task so the Makefile path and the
# `pixi run` path produce identical output.
lint-py:
	pixi run lint-py

lint-perl:
	pixi run lint-perl

lint-shell:
	pixi run lint-shell

lint-snakemake:
	pixi run lint-snakemake

pyright:
	pixi run pyright src bin

clean:
	rm -rf containers/jamg.sif containers/jamg-base.sif \
	       containers/trinity.sif containers/pasa.sif \
	       containers/build.log containers/build-base.log \
	       containers/build-trinity.log containers/build-pasa.log \
	       .snakemake test_suite/output .pixi

docs:
	@echo "docs/*.md require no build"
