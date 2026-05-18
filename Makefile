SHELL := /bin/bash

.PHONY: sif test test-container smoke lint lint-py lint-perl lint-shell lint-snakemake pyright clean docs

REPO_ROOT := $(shell git rev-parse --show-toplevel)

# --fakeroot requires subuid/subgid mapping; fall back to plain build
# (which may need root). If neither works, the user must rebuild on a
# host with apptainer fakeroot support OR pull a prebuilt SIF.
APPTAINER_BUILD_FLAGS := $(shell apptainer build --help 2>/dev/null | grep -q -- --fakeroot && echo --fakeroot)

sif: containers/jamg.sif

# Base SIF: debian-slim + apt bootstrap + statically linked micromamba.
# Rebuilt only when containers/jamg-base.def changes.
# `-F` forces overwrite of an existing SIF (apptainer build refuses to clobber
# by default -- reviewer I1).
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

# test-container: just the container presence test (Phase 1 deliverable).
# Always runs.
test-container: sif
	bash tests/container/test_tools_present.sh containers/jamg.sif

# test: full unit + integration test suite. Phase 2+ adds tests/golden/ + the
# Snakemake rule tests + tests/no-legacy/. For now, only run what exists.
# Skip-with-notice for forward-decl paths so this target stays green between
# phase landings (review B2).
test: test-container
	@[ -f tests/golden/test_validator.t ] && prove tests/golden/test_validator.t || \
	  echo "SKIP: tests/golden/test_validator.t (Phase 4 deliverable)"
	@[ -f tests/golden/test_prepare_golden_genes.t ] && prove tests/golden/test_prepare_golden_genes.t || \
	  echo "SKIP: tests/golden/test_prepare_golden_genes.t (Phase 4 deliverable)"
	@[ -f tests/no-legacy/test_legacy_tools_gone.sh ] && bash tests/no-legacy/test_legacy_tools_gone.sh || \
	  echo "SKIP: tests/no-legacy/test_legacy_tools_gone.sh (Phase 5 deliverable)"

smoke: sif                             # slow: full DAG on mini-genome (~30 min)
	@[ -f tests/e2e/test_full_dag.sh ] && bash tests/e2e/test_full_dag.sh || \
	  echo "SKIP: tests/e2e/test_full_dag.sh (Phase 6 deliverable)"

lint: lint-py lint-perl lint-shell lint-snakemake

# All lint targets delegate to the pixi task so the Makefile path and the
# `pixi run` path produce identical output (review B4). The pixi tasks handle
# phase-aware path existence checks internally.
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
	       containers/build.log containers/build-base.log \
	       .snakemake test_suite/output .pixi

docs:
	# rendered docs are markdown; nothing to build, but kept for symmetry
	@echo "docs/*.md require no build"
