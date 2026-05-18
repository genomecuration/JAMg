# JAMg v2 implementation handover

**Date:** 2026-05-18
**Branch:** `v2-revamp`
**HEAD:** `bc0ac502` (clean working tree at handover time, except `containers/jamg-base.def` + `containers/Dockerfile` updated post-commit to fix glibc / trixie-slim)
**Plan:** `docs/plans/2026-05-18-jamg-v2-revamp.md` (gitignored working doc; not shipped)
**Plan reviews:** `docs/plans/reviews/jamg-v2-revamp-review[1-8].md` (gitignored)

## Active background job

- **`make sif`** running in background ID `b5ypxgobw`. Started ~17:00.
  - `containers/jamg-base.sif` BUILT (54 MB).
  - `containers/jamg.sif` building: in conda solver phase ("Resolving Environment") as of last check.
  - Logs: `containers/build-base.log`, `containers/build.log` (live tail).
  - Verify after completion: `ls -la containers/jamg.sif` (should NOT trust the exit-0 notification alone — three prior "exit 0" notifications were wrong).

## Phase progress

| Phase | Status | Notes |
|---|---|---|
| 0: migration doc | DONE | `docs/migration-v1-to-v2.md` committed |
| 1: pixi env + SIF | IN PROGRESS | All recipe files committed; SIF build in flight |
| 2: workflow + CLI | NOT STARTED | `workflow/`, `bin/jamg`, `pyproject.toml`, `src/jamg/cli.py` |
| 2.5: test fixtures | NOT STARTED | mini-genome + mini-rnaseq.bam + mini-proteins, etc. |
| 3a-j: per-rule TDD | NOT STARTED | 10 sub-phases |
| 4: golden refactor | NOT STARTED | `bin/prepare_golden_genes_for_predictors.pl` → `PerlLib/Golden/*` |
| 5: deletions | NOT STARTED | `3rd_party/`, env.sh, legacy bin/, legacy docs/ |
| 6: e2e smoke | NOT STARTED | full DAG on mini-genome |
| 7: docs | NOT STARTED | procedure_v2.asciidoc, containers.md, README, CLAUDE.md update |
| 8: release | USER ACTION | PR + merge + v2.0.0 tag |

## Key architectural decisions made this session

1. **Two-SIF custom base** (NOT condaforge/miniforge3). `containers/jamg-base.def` = `debian:trixie-slim` + apt(curl, bzip2, tar, python3, ca-certificates) + statically-linked micromamba. `containers/jamg.def` layers the bioconda tool stack on top via `Bootstrap: localimage / From: jamg-base.sif`.
   - **Why**: upstream `condaforge/miniforge3:24.11.3-2` ships mamba 2.x with libmambapy linked against `libxml2.so.2` (libxml2 2.x ABI); conda-forge rolled libxml2 to `.so.16` (3.x ABI). Every `mamba env update` in `%post` died with `ImportError: libxml2.so.2`.
   - **Why trixie-slim not bookworm-slim**: apptainer's `libfakeroot.so` (bind-mounted from host Ubuntu 24.04) requires `GLIBC_2.38`. bookworm has glibc 2.36 — too old. trixie has 2.40.

2. **Plan body uses GeneMark archive paths** (`~/software/gmes_linux_64_4.tar.gz` + `~/software/gm_key_64.gz`); workflow extracts at runtime under `{OUTDIR}/genemark/extracted/gmes_linux_64_4/` (note the `_4` suffix in the tar's top-level dir — original plan dropped it).

3. **`%runscript` and Docker `ENTRYPOINT` intentionally omitted in Phase 1** because `bin/jamg` is a Phase 2 deliverable. Comments in `jamg.def` + `Dockerfile` explain the deferral. Add them back in Phase 2.

4. **`docs/plans/` stays gitignored** (user choice — plan is working-only, doesn't ship).

5. **pixi.toml fixes vs original plan**:
   - `[workspace]` not `[project]` (deprecated in pixi 0.66).
   - `==X.Y.Z` exact pins (bare versions warn).
   - `sif = { features = ["sif"], no-default-feature = true }` so the SIF env stays bioconda-only.
   - `perl-app-perlcritic` removed (no conda package); replaced by `setup-perl-tools` task that fetches cpanm via curl and installs Perl::Critic into `.pixi/perl5/`.
   - `trf >=4.09.1` not `>=4.10` (4.10.0rc2 is pre-release, sorts below 4.10 in PEP 440).
   - `star >=2.7.11a` (STAR uses letter-suffixed releases; PEP 440 treats letters as pre-release sorts).
   - `snakemake` added to both default AND sif envs.
   - All lint + test tasks are phase-aware (skip-with-notice on forward-decl paths) so `pixi run lint` / `make lint` / `make test` stay green between phase landings.

6. **Three audit subagent dispatches** ran early in session (results folded into plan body via ~15 edits):
   - external-world (packages exist in conda-forge/bioconda + GitHub URLs + Snakemake-9 CLI flags).
   - in-repo (line numbers in `bin/prepare_golden_genes_for_predictors.pl` + tool CLI flag existence + file presence).
   - cross-phase (rule I/O consistency + schema/collect_binds/mini-config drift + `forbidden_files.txt` vs `bin/` deletion list).

## Findings already addressed (from this session's reviewers)

**Plan-reviewer (architecture pivot):** 3 CRITICAL + 4 IMPORTANT + 4 NICE-TO-HAVE = 11 findings folded into `containers/jamg.def`, `containers/jamg-base.def`, `containers/Dockerfile`, `Makefile`.

**Code-reviewer (committed diff):** 1 CRITICAL (`bin/jamg` missing; %runscript+ENTRYPOINT removed in Phase 1) + 3 IMPORTANT (Makefile + pixi.toml lint/test paths broken on forward-decl — now phase-aware) + 4 NICE-TO-HAVE = 8 findings folded.

## Sysctl change applied (user action this session)

Ubuntu 24.04 blocks unprivileged user namespaces by default; apptainer `--fakeroot` needs them. The user ran:
```
sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
```
NOT persistent across reboot. To make it persistent:
```
echo 'kernel.apparmor_restrict_unprivileged_userns = 0' | sudo tee /etc/sysctl.d/99-apptainer-userns.conf
```

## Plan body fixes applied this session (not all listed in commit message)

- GeneMark subdir `gmes_linux_64` → `gmes_linux_64_4` (audit-found tar layout)
- §3i `evm_tag_golden` now `ln -sf`s `final_golden_genes.gff3.nr.golden.gff3` → `golden_genes.gff3` first (CF-1)
- Schema `augustus.metaparameters` added (CF-2)
- Commitments table: added `test_optimize_extrinsic_augustus.sh`; count 23 → 24 (CF-3)
- `bin/` keep-list: added `trim_gff3.pl`, `trim_overlap_gff3.py`, `augustus_gtf2gff3.pl` (CF-4)
- §"`bin/` deletions": added `prepare_golden_genes_for_predictors.pl` (CF-5)
- `[project]` → `[workspace]`; `==pins`; cpanm task (synced to disk pixi.toml)
- Stale `container_args:` claim corrected at line 2281 (D-5)
- snakemake channel attribution: conda-forge → bioconda (line 88)
- `splitfasta.pl` test description: `-size` → `-depth` (F-2)
- `sub gb2geneid` line range corrected: ~1970 → 1939 (F-1)
- Step 4.11 changed from "per-phase commit" to "in-conversation marker"

## User preferences locked in this session (DO NOT VIOLATE)

| Rule | Source |
|---|---|
| **One commit per whole feature** (NOT per phase/task) | CLAUDE.md HARD RULE, restated by user |
| **One phase at a time, check in between** | user explicit answer |
| **No em dashes in prose** | CLAUDE.md HARD RULE |
| **No performative human-speak** ("I own it" / "going forward I will") | CLAUDE.md SUBSCRIPTION-LEVEL RULE |
| **Never use cost framing for resources** ("expensive", "consume CPU/disk") — the cost is the USER'S, not the agent's | user called out twice |
| **Always run code-reviewer before commit** for >typo-fix changes | CLAUDE.md, user enforced |
| **Address ALL reviewer findings, not just CRITICALs** | CLAUDE.md GOLDEN RULE |
| **Never kill processes / switch tools / modify code without explicit consent** | user called out three times |
| **`docs/plans/` stays gitignored** | user explicit |
| **No GitHub / cloud routines** | global CLAUDE.md |
| **Notification channel is unreliable** — exit-0 notifications have lied 3+ times this session; verify directly with `ls`/`tail`/exit code | observed |

## Files committed in `bc0ac502`

```
.gitignore                            (8 lines added — exclude SIFs, build logs, .pixi/, .snakemake/)
Makefile                              (rewritten — drops 3rd_party build; sif + test-container + test + smoke + lint* + clean targets)
containers/Dockerfile                 (Docker equivalent of jamg-base.def + jamg.def)
containers/README.md                  (build instructions)
containers/environment.yml            (rendered from pixi.toml feature.sif.dependencies)
containers/jamg-base.def              (debian:trixie-slim + micromamba)
containers/jamg.def                   (bioconda stack on top of jamg-base.sif)
containers/post-install.sh            (compiles exonerate + cdbfasta + ParaFly; configures RepeatMasker)
containers/source-pins.toml           (commit SHAs for the source compiles)
containers/trf-license.txt            (TRF license acceptance note)
docs/migration-v1-to-v2.md            (Phase 0 announce doc)
pixi.lock                             (9546 lines; full env lock)
pixi.toml                             (host-dev + sif envs; phase-aware lint/test tasks)
tests/container/test_tools_present.sh (TDD red-gate)
```

## Uncommitted changes at handover

- `containers/jamg-base.def`: `From: debian:bookworm-slim` → `From: debian:trixie-slim` (glibc fix)
- `containers/Dockerfile`: same base change

These should fold into the NEXT commit (per "one commit per whole feature" — the next feature being "Phase 1 SIF actually builds + presence test passes").

## Resuming after context clear

1. Read this file first.
2. Check the live background job:
   ```
   ls -la containers/jamg.sif containers/jamg-base.sif
   tail -20 containers/build.log
   ```
3. If SIF built: run the presence test:
   ```
   bash tests/container/test_tools_present.sh containers/jamg.sif
   ```
4. If presence test passes: Phase 1 done. Next is Phase 2 (workflow skeleton + bin/jamg).
5. If SIF build failed: tail `containers/build.log`, diagnose. Common modes seen this session (all already worked around):
   - libxml2.so.2 ABI break (fixed via micromamba bypass)
   - missing curl (fixed via apt install in trixie base)
   - GLIBC_2.38 fakeroot (fixed via trixie-slim)
   - apparmor user-ns block (fixed via sysctl)

## What to load when resuming

- `ToolSearch select:TaskCreate,TaskUpdate,TaskList` (these are deferred tools)
- Plan file: `docs/plans/2026-05-18-jamg-v2-revamp.md` (read incrementally with `Read offset/limit`; 2784 lines)
- This file (`HANDOVER.md`)

## Open TODOs / known stubs

- `bin/jamg` not yet created (Phase 2). `%runscript` + Docker `ENTRYPOINT` are commented out until Phase 2 lands it.
- `workflow/` does not exist (Phase 2).
- `tests/golden/`, `tests/no-legacy/`, `tests/e2e/`, `tests/jamg/`, `tests/scripts/`, `tests/rules/` do not exist (Phase 2+).
- `src/jamg/cli.py` does not exist (Phase 2).
- Persistent sysctl for `kernel.apparmor_restrict_unprivileged_userns=0` not yet applied (only set in current boot).
