# JAMg v2-revamp  -  handover for next agent

## Read first

- `docs/procedure_v2.asciidoc`  -  operational reference; covers test architecture + Layout B Slurm dispatch.
- `docs/plans/2026-05-20-fast-rule-tests.md`  -  feature plan (all tasks complete; later issues uncovered, see below).
- `tmp/temp.md`  -  user's personal v1-pipeline recipe (authoritative for "what does v1 do at step X?"; shell history, not a spec).
- This file  -  current state of the EVM-blocker fix + container rewiring work.

## Repo state

```
JAMg branch: v2-revamp
JAMg HEAD:   eeab866e (working tree has extensive uncommitted changes; see below)
```

Submodules (registered, working trees populated, parent pointer staged):

| Submodule | Path | Local HEAD | Pushed to origin? |
|---|---|---|---|
| alpapan/EVidenceModeler | `containers/src/EVidenceModeler` | `905ebe2` on `fix/negative-offset-warn` (= origin/master, ff-merged) | **YES** (master) |
| alpapan/PASApipeline | `containers/src/PASApipeline` | `dc5497d` on `fix/gcc14-implicit-int` (1 commit ahead of origin/master `cc1f8d7`) | **NO**  -  needs explicit user push approval |

Uncommitted in JAMg (large diff; staged-as-add for submodules + rm_libs LFS, modified for code/docs):

| File | What changed |
|---|---|
| `.gitmodules` | EVM + PASApipeline submodule entries |
| `.gitattributes` | LFS pattern `containers/rm_libs/**` |
| `Makefile` | RM_LIB_HOST removed; rm_libs file prereq; `--bind` removed |
| `containers/jamg.def` | %setup removed; %files adds rm_libs + EVM src |
| `containers/pasa.def` | Fully rewritten: `From: jamg-base.sif`, %files overlay, %post apt+bioconda+make |
| `containers/post-install.sh` | EVM block uses /opt/jamg-src cp (no git clone); python3 shim; RM section comment + FATAL message rewritten |
| `containers/source-pins.toml` | `[evidencemodeler]` is now `submodule = ...` |
| `containers/src/EVidenceModeler` | submodule (89 LFS-tracked files at `containers/rm_libs/`) |
| `containers/src/PASApipeline` | submodule |
| `containers/rm_libs/**` | 89 LFS-tracked files (667 MB of RepeatMasker libs) |
| `Makefile`, `tools/regen_snapshot.sh`, `tools/stage_fixtures.sh`, `tests/rules/test_evm.sh`, `tests/rules/test_ogs.sh`, `tests/e2e/test_full_dag.sh` | `--mtime-only` flag purged |
| `README.md`, `docs/containers.md`, `docs/procedure_v2.asciidoc`, 7 `tests/rules/test_*.sh`, `tests/golden/test_prepare_golden_genes.t` | Stale `RM_LIB_HOST=...` references replaced |
| `src/jamg/cli.py` | `--mtime-only` option + handling removed (22 lines deleted) |
| `tools/stage_fixtures.sh` | `-h` flag on touch (don't follow symlinks); mtime "now" not 2038 |
| `test_suite/mini-genome.fasta` | mtime restored to current (symlink-target poisoning fixed) |

SIF state:

| SIF | Built when | Verified contents |
|---|---|---|
| `containers/jamg-base.sif` | 2026-05-20 16:47 (unchanged) | Debian trixie, perl 5.40, BioPerl, micromamba at /usr/local/bin/micromamba, /opt/conda empty |
| `containers/jamg.sif` | 2026-05-20 22:44 (current) | EVM patch at evidence_modeler.pl:3517 (warn, not confess); python3 shim at /opt/jamg/bin/python3 → /opt/conda/bin/python3.12; RepeatMasker libs at /opt/conda/share/RepeatMasker/Libraries/; no /rm_lib_host stub |
| `containers/pasa.sif` | 2026-05-20 23:16 (current) | PASA at /opt/PASApipeline (fork content + compiled bin/), blat+pblat from bioconda at /opt/conda/bin/, fasta symlink at /usr/local/bin/fasta → /usr/bin/fasta36, PASAHOME=/opt/PASApipeline |
| `containers/trinity.sif` | unchanged | trinityrnaseq/trinityrnaseq upstream Docker (per containers/trinity.def) |

## What works (verified)

- EVM patch confirmed inside jamg.sif (`apptainer exec containers/jamg.sif grep "WARNING: offset to augment" /opt/jamg/share/EVidenceModeler/EvmUtils/evidence_modeler.pl` returns line 3517).
- jamg.sif RepeatMasker stack: famdb.py runs cleanly via `env python3` (shim resolves it to conda's python with h5py 3.16.0).
- pasa.sif PASA stack: `Launch_PASA_pipeline.pl --version` returns 2.5.3; compiled binaries (pasa, slclust, cdbfasta, cdbyank, mdust, psx, seqclean, seqclean.psx, trimpoly, cln2qual) all present at `/opt/PASApipeline/bin/`; aligners gmap/blat/pblat/minimap2/fasta/samtools all on PATH inside the SIF.
- Submodule sub-submodules: PASApipeline's `pasa-plugins/{cdbtools,seqclean,slclust,transdecoder}` populated locally (transdecoder via `git submodule update --init --recursive`).
- EVM fix branch already on alpapan/EVidenceModeler origin/master (905ebe2).

## What's blocked

`tests/rules/test_evm.sh`, `tests/rules/test_ogs.sh`, `tests/e2e/test_full_dag.sh`, and `tools/regen_snapshot.sh` all fail in snakemake 9.21 even with the EVM patch in place. Two distinct issues:

### Issue 1: snakemake 9.21 `ancient()` does not suppress the in-session cascade

`workflow/rules/genemark.smk:84-87` marks all three `genemark` rule inputs as `ancient()`. The rule comment (lines 72-83) claims this prevents snakemake's "Input files updated by another job" cascade. **It does not, in snakemake 9.21.**

Symptom: pre-staged `genemark.gff3` + `genemark.gtf` exist at far-future mtime (2038-01-15), `.preflight.ok` is pre-staged too  -  but when upstream rules (`rnaseq_hints`, `repeats_merge`, `genemark_preflight`) run in the same session and produce their outputs, snakemake reschedules `genemark` with reason "Input files updated by another job". The 100 kb fixture is too small for GeneMark-ES self-training (`error, input sequence size is too small data/training.fna: 70478`), so the re-execution fails.

`--rerun-triggers mtime` does not help (previously plumbed as `--mtime-only`; both removed because they had no effect on this cascade and were misleading).

### Issue 2: snakemake clock-skew detector deletes new outputs

When `evm_tag_genemark` (or any rule whose input is at mtime 2038) produces an output at current-time, snakemake's clock-skew detector fires: `Output ... has older modification time (2026-05-20 ...) than input ... (2038-01-15 ...). This could indicate a clock skew problem ... Removing output files of failed job since they might be corrupted`. snakemake deletes the rule's just-produced output and marks the rule failed.

The 2038 strategy in `tests/rules/test_evm.sh:48`, `tests/rules/test_ogs.sh:45`, `tests/e2e/test_full_dag.sh:57`, and `tools/regen_snapshot.sh:25` still uses `touch -d '2038-01-15'` for the genemark fixture. `tools/stage_fixtures.sh` was already converted to current-time mtime; the others have not been converted because they pre-stage outputs of rules whose inputs would then be in the future.

### Diagnosis-but-not-fix

`--consider-ancient genemark=preflight,rnaseq_hints,softmasked` (snakemake 9 CLI flag) is the modern equivalent of in-rule `ancient()` and may suppress the cascade where the in-rule decorator does not. Not tried yet; would need plumbing through `src/jamg/cli.py` (no passthrough mechanism currently).

## Next actions (ordered)

1. **Diagnose the snakemake 9 cascade.** Try `--consider-ancient genemark=preflight,rnaseq_hints,softmasked` directly via a hand-rolled snakemake invocation (bypassing `bin/jamg`); if it suppresses the cascade, plumb it through `src/jamg/cli.py`. If it does not, the architectural assumption in `workflow/rules/genemark.smk:72-83` is wrong and the workflow needs restructuring (e.g., make genemark.gff3 a config-supplied input, not a rule output).
2. **Resolve the 2038 mtime issue in `regen_snapshot.sh`, `test_evm.sh`, `test_ogs.sh`, `test_full_dag.sh`.** Either drop the 2038 touch (use current mtime) AND fix the cascade so pre-staging still works, or switch to a different pre-staging mechanism (e.g., `--touch <target>` to mark outputs as up-to-date in snakemake metadata).
3. **Push PASApipeline `fix/gcc14-implicit-int` to origin (alpapan/PASApipeline).** Local commit `dc5497d` is the GCC 14 implicit-int fix for `pasa-plugins/seqclean/psx/psx.c`. Requires explicit per-invocation user approval before `git push origin fix/gcc14-implicit-int` or `git push origin master` (after ff-merge).
4. **Run integration tests once 1+2 are resolved**: `tests/rules/test_evm.sh`, `tests/rules/test_ogs.sh`, `tests/e2e/test_full_dag.sh`. Each must produce its declared assertion (EVM.gff3 ≥ 1 EVM gene row, OGS GFF non-empty, full DAG completes).
5. **Dispatch `feature-dev:code-reviewer` on the full diff** (multi-file, ~hundreds of lines once committed). Required by the "code-review before commit" rule for diffs >30 lines.
6. **Address every reviewer finding** (CRITICAL through Question for Author  -  see standing rules below).
7. **Single commit** via pathspec form (no `git add -A`). Author: Alexie Papanicolaou.

## Standing rules (project + user instructions)

- **Never push without explicit per-invocation approval.** Prior approval for a different branch/commit does not carry over.
- **Never commit with `git add -A` / `git add .`**  -  pathspec form only.
- **Never wipe uncommitted changes without explicit approval**  -  including `git checkout <path>`, `git reset --hard`, `git clean`, etc. See the global subscription-level rule.
- **Code-review before any first execution of a new/modified script** (shell, srun, sbatch, snakemake rule shell, apptainer exec). Self-review for diffs ≤30 lines (`git diff` + targeted greps); dispatch `feature-dev:code-reviewer` for larger diffs.
- **Address every reviewer finding**  -  CRITICAL, IMPORTANT, NICE-TO-HAVE, Suggestion, Question. Open every review-return status with the full severity count.
- **Agent prompts terse** (≤500 words); reviewers Write findings to disk only if they have the Write tool (`feature-dev:code-reviewer` does NOT  -  ask for inline findings, ≤15 per dispatch).
- **One commit per whole feature**, not per task.
- **No em dashes in prose** (project HARD RULE). En dashes for numerical ranges OK. Hyphens for compound words OK.
- **No phase gating**  -  execute multi-step plans autonomously without "ready for next phase?" pauses.
- **Brainstorm before committing to a fix direction.** Don't shortcut with A/B/C framing; investigate root cause first.
- **No performative human-speak**  -  no "I own it" / "going forward" / forward self-promises / contrition-theatre.

## Architecture invariants

These are NOT historical decisions; they are how the system currently IS:

- **EVM source** is consumed via the `alpapan/EVidenceModeler` submodule at `containers/src/EVidenceModeler`. The SIF build copies the working tree to `/opt/jamg-src/EVidenceModeler` via `containers/jamg.def`'s `%files`, then `containers/post-install.sh` `cp -a`'s into `$JAMG/share/EVidenceModeler/`. No `git clone` at SIF-build time. The fork carries a `warn`-not-`confess` patch at `evidence_modeler.pl:3517` for negative-offset noncoding-score cases.
- **PASApipeline source** is consumed via the `alpapan/PASApipeline` submodule at `containers/src/PASApipeline`. The SIF build copies the working tree to `/opt/PASApipeline` via `containers/pasa.def`'s `%files`, then `%post` runs apt/cpanm/bioconda installs and `make` to compile pasa_cpp + plugins. PASAHOME = /opt/PASApipeline (NOT /usr/local/src/PASApipeline). The fork carries upstream cherry-picks + a GCC 14 implicit-int fix at `pasa-plugins/seqclean/psx/psx.c`.
- **pasa.sif bootstraps from `jamg-base.sif`**, not from the upstream `pasapipeline:2.5.3` Docker image. Reason: trixie's glibc-2.40 is compatible with apptainer's libfakeroot, so `%post` works under `--fakeroot`. The upstream image's glibc-2.31 makes `%post` impossible.
- **RepeatMasker libraries** are vendored at `containers/rm_libs/` via `git-lfs` (`.gitattributes` pattern `containers/rm_libs/** filter=lfs ...`). `containers/jamg.def`'s `%files` mounts them at `/rm_lib_host`; `containers/post-install.sh` copies them into `/opt/conda/share/RepeatMasker/Libraries/`; the staging dir is `rm -rf /rm_lib_host` in `containers/jamg.def`'s `%post`. **No `--bind`-mount** in the Makefile; `git lfs pull` on a fresh clone is sufficient.
- **`python3` shim** at `/opt/jamg/bin/python3` → `/opt/conda/bin/python3`. The SIF's `%environment` puts `/usr/bin` before `/opt/conda/bin` (so `env perl` finds system perl 5.40 with apt BioPerl); the shim puts conda's python3 in front of system python3 for scripts whose shebang is `#!/usr/bin/env python3` and that need h5py (RepeatMasker's famdb.py, util/RM2Bed.py).
- **`fasta` binary** is `/usr/local/bin/fasta` → `/usr/bin/fasta36`. The Debian fasta3 apt package ships `fasta36` only; PASA expects `fasta`.
- **No `--mtime-only` flag.** Removed from `src/jamg/cli.py` and all callers. snakemake's native `--rerun-triggers mtime` was what `--mtime-only` translated to; that flag does NOT suppress the in-session cascade in snakemake 9, so the wrapper was misleading.
- **`stage_fixtures.sh` uses current-time mtime + `-h`** flag on `touch`. `-h` prevents following symlinks (snapshot/repeats/*/mini-genome.fasta is a symlink to test_suite/mini-genome.fasta; without `-h`, touch poisons the source fixture's mtime). Current-time mtime avoids the clock-skew detector firing on downstream rule outputs.

## Slurm + I/O guidance

| Work | Where |
|---|---|
| `apptainer build` (any SIF) | **lazebnik directly**  -  it's I/O, sbatching adds NFS round-trips |
| `git clone` / `git submodule update` / `git lfs pull` / file mv | lazebnik directly |
| Tarball extract, SIF copy, sha256 check, log read, `apt-cache` inside SIF | lazebnik directly |
| `tools/build_fixtures.sh` (downloads + polyester + extracts) | lazebnik directly |
| `tools/regen_snapshot.sh` (drives snakemake → predictor compute) | sbatch (20 CPU, 48 GB mem, partition main) |
| `tests/rules/test_*.sh` (each invokes snakemake → predictor runs) | sbatch (20 CPU, 32 GB mem) |
| `apptainer exec ... bin/jamg run ...` (the snakemake driver itself is light, but rules under it call gcc/PASA/Augustus/etc.) | sbatch |
| PASA / Trinity / STAR / BLAST / RepeatMasker / Augustus / GeneMark / EVM (any direct invocation) | sbatch |

After every sbatch: schedule a +1 min check first, then a +5 min check if still running. Single long sleep misses fast failures.

```bash
# Pattern: short check first
( sleep 60 && sacct -j "$JID" -X --format=JobID,State,ExitCode,Elapsed --noheader -P && \
  tail -30 "$TMP/<prefix>-<jobname>-$JID.err" ) > "$TMP/sanity-1min-$JID.log" 2>&1 &
# Then long
( sleep 300 && ... ) > "$TMP/sanity-5min-$JID.log" 2>&1 &
```

Do NOT set `--time` on sbatch; partition default is INFINITE.

`$TMP=/mnt/lazebnik/home/30042108/tmp` is NFS-shared. Never write to `/tmp` (per-host ramdisk, invisible across nodes).

`/mnt/lazebnik/home/30042108/.apptainer/cache` is NFS-shared too; one `apptainer cache clean --force` covers all nodes. `/mnt/lazebnik/home/30042108/tmp/build-temp-*` directories accumulate per failed build at ~5 GB each; clean periodically via `rm -rf /mnt/lazebnik/home/30042108/tmp/build-temp-*` (only when no build is active).

## Lessons embedded as forward-looking rules

- **Validate apt package availability on the EXACT base before submitting a build.** `apptainer exec --writable-tmpfs containers/jamg-base.sif bash -c 'apt-get update >/dev/null 2>&1 && for p in pkg1 pkg2; do apt-cache madison "$p" | head -1; done'` takes seconds; submitting a 5-min build that fails on `Package 'blat' has no installation candidate` is wasteful.
- **bioconda packages need conda-forge** (zlib, libpng, libgcc-ng deps). Always list `-c conda-forge -c bioconda` (conda-forge first).
- **Apptainer `%files` is host-side; `%post` enters the container.** On a base with old glibc (pasapipeline:2.5.3 → glibc-2.31), `%post` fails under `--fakeroot` because libfakeroot needs GLIBC_2.33+. `%files`-only overlays still work. For full from-source builds with `%post`, bootstrap from a modern base (jamg-base.sif = trixie = glibc-2.40+).
- **`find ... -exec touch -d ... {} +` without `-h` follows symlinks** and modifies the targets, not the symlinks. Always pass `-h` to `touch` (or `-not -type l` to `find`) when staging directories that may contain symlinks to source fixtures.
- **snakemake 9's `ancient()` decorator in a rule does NOT fully suppress the in-session "Input files updated by another job" cascade.** Use `--consider-ancient RULE=INPUTITEMS` at invocation time, or restructure the rule.
- **snakemake's clock-skew detector deletes just-produced outputs** when their inputs have future mtime. Avoid setting input mtimes to 2038 if any downstream rule will produce new output during the run.
- **GCC 14 (trixie default) escalates implicit-int and implicit-function-declaration from warnings to hard errors.** Legacy K&R-style C (e.g., PASA's `pasa-plugins/seqclean/psx/psx.c`) needs explicit return types. Patch the source rather than passing `-Wno-error=implicit-int`  -  modernize, don't silence.
- **The Debian fasta3 apt package ships `fasta36` only, not `fasta`.** Symlink if PASA-like tools expect the plain name.
- **The `blat` binary is bioconda-only on Debian trixie.** Install via `micromamba install -c conda-forge -c bioconda blat`.

## Where things are documented

- `tools/temp.md`: user's v1 hand-driven pipeline recipe (cdna of run commands).
- `docs/procedure_v2.asciidoc`: v2 Snakemake reference.
- `docs/containers.md`: container build invocation + image table.
- `containers/jamg-base.def`, `containers/jamg.def`, `containers/pasa.def`, `containers/trinity.def`: SIF definitions.
- `containers/post-install.sh`: jamg.sif source-builds (Augustus, EVidenceModeler, TransDecoder, cdbfasta, ParaFly, RepeatMasker config, python3 shim).
- `containers/source-pins.toml`: build pins (informational; submodule SHAs are authoritative).
- `workflow/rules/*.smk`: Snakemake rules.
- `workflow/profiles/slurm/`: Slurm executor profile.
- `~/.claude/projects/-mnt-lazebnik-home-30042108-software-gits-JAMg/memory/`: per-user memory files (slurm-vs-I/O rule, code-review dispatch scope, v1-recipe pointer, etc.).
