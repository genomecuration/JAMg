# JAMg v2-revamp  -  handover for next agent

## Read first

- `docs/procedure_v2.asciidoc`  -  operational reference; covers test architecture + Layout B Slurm dispatch.
- `docs/plans/2026-05-20-fast-rule-tests.md`  -  feature plan (all tasks complete; later issues uncovered, see below).
- `tmp/temp.md`  -  user's personal v1-pipeline recipe (authoritative for "what does v1 do at step X?"; shell history, not a spec).
- This file  -  current state of the EVM-blocker fix + container rewiring work.

## Repo state

```
JAMg branch: v2-revamp
JAMg HEAD:   3de6d229 (working tree clean; 1 commit ahead of origin)
```

Submodules (registered, parent pointer committed):

| Submodule | Path | Local HEAD | Pushed to origin? |
|---|---|---|---|
| alpapan/EVidenceModeler | `containers/src/EVidenceModeler` | `905ebe2` on `master` | **YES** (master) |
| alpapan/PASApipeline | `containers/src/PASApipeline` | `dc5497d` on `master` (ff-merged from `fix/gcc14-implicit-int`); 1 commit ahead of origin/master `cc1f8d7` | **NO**  -  needs explicit user push approval |

What landed in commit `3de6d229`:

| File / area | Change |
|---|---|
| `.gitmodules` | EVM + PASApipeline submodule entries |
| `.gitattributes` | LFS pattern `containers/rm_libs/**` |
| `Makefile` | RM_LIB_HOST removed; rm_libs file prereq; `--bind` removed |
| `containers/jamg.def` | %setup removed; %files adds rm_libs + EVM src |
| `containers/pasa.def` | Fully rewritten: `From: jamg-base.sif`, %files overlay, %post apt+bioconda+make |
| `containers/post-install.sh` | EVM block uses /opt/jamg-src cp (no git clone); python3 shim; RM section comment + FATAL message rewritten |
| `containers/source-pins.toml` | `[evidencemodeler]` is now `submodule = ...` |
| `containers/src/EVidenceModeler` | submodule pointer (905ebe2) |
| `containers/src/PASApipeline` | submodule pointer (dc5497d) |
| `containers/rm_libs/**` | 89 LFS-tracked files (667 MB of RepeatMasker libs) |
| `src/jamg/cli.py` | `--mtime-only` option + handling removed (22 lines deleted) |
| `tools/stage_fixtures.sh` | `-h` flag on touch (don't follow symlinks); mtime "now" not 2038 |
| `tools/regen_snapshot.sh`, `tests/rules/test_evm.sh`, `tests/rules/test_ogs.sh`, `tests/e2e/test_full_dag.sh` | `--mtime-only` flag purged |
| `README.md`, `docs/containers.md`, `docs/procedure_v2.asciidoc`, 7 `tests/rules/test_*.sh`, `tests/golden/test_prepare_golden_genes.t` | Stale `RM_LIB_HOST=...` references replaced |
| `HANDOVER.md` | this file (force-added; was gitignored) |

Not in `3de6d229` (intentionally; require user direction):

- `test_suite/mini-genome.fasta`: mtime was poisoned to 2038 by the latent symlink-follow bug in `stage_fixtures.sh`. Restored to current mtime on disk. File is gitignored (regeneratable via `tools/build_fixtures.sh`); not part of this commit. Future runs will not re-poison because the `-h` flag now prevents it.
- Push of either submodule's master to origin.
- Push of `genomecuration/JAMg` `v2-revamp` to origin.

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

## Test infrastructure fix (applied 2026-05-21)

The blocker on `tests/rules/test_evm.sh`, `tests/rules/test_ogs.sh`, `tests/e2e/test_full_dag.sh`, and `tools/regen_snapshot.sh` is RESOLVED for the pre-staged scenarios. The original diagnosis (snakemake 9 `ancient()` bug) was incorrect; the actual root cause was the `touch -d '2038-01-15'` strategy combined with `stage_fixtures.sh`'s non-atomic `find -exec touch`.

### Actual root cause

**Primary (always fires):** snakemake's clock-skew detector deletes the just-produced output of `evm_tag_genemark` because its input `genemark.gff3` was at `2038-01-15` (set by `touch -d '2038-01-15'` in test_evm.sh, test_ogs.sh, test_full_dag.sh, regen_snapshot.sh). Per snakemake's source: when an output is produced with an mtime older than any input, the output is deleted and the job fails. Empirical reproduction: see `$TMP/red-test_evm.log`.

**Secondary (state-dependent):** `tools/stage_fixtures.sh:44`'s `find ... -exec touch -h {} +` uses the `+` batching form. On large staged subtrees (pasa/ has many sqlite checkpoint files), multiple `touch` invocations within one `find` walk produce a multi-second mtime spread; from a fresh `rm -rf` baseline the spread is sub-second and harmless, but with stale `test_suite/output` it grows to 5-7 seconds, sufficient to make some staged inputs newer than other staged outputs and trigger the spurious "Updated input files" cascade.

`ancient()` in `genemark.smk:84-87` IS working correctly in snakemake 9.21 for the test_evm and test_ogs scenarios (genemark is correctly excluded from the rerun list whenever its output is pre-staged). Whether it suppresses the in-session "Input files updated by another job" cascade in the test_full_dag / regen_snapshot full-DAG scenarios (where rnaseq_hints, repeats_merge, genemark_preflight run fresh) is the go/no-go gate not yet empirically tested at the time of this commit; see Next actions below.

### Fix applied

Source-fixture mtimes are now the temporal floor. `tools/build_fixtures.sh` pins all `test_suite/mini-*` files to `2020-01-01` at the end of every fixture rebuild. Every derived artifact (snapshot, staged outputs, fresh rule outputs) is naturally newer. The non-atomic `find -exec touch` in `stage_fixtures.sh` and the `touch -d '2038-01-15'` lines in the 4 test/regen scripts are deleted. The `ancient()` decorator in `genemark.smk` is retained (the rule comment is corrected to drop the misdiagnosis claim).

Bootstrap: the existing on-disk fixtures were re-stamped to `2020-01-01` via a one-shot `touch -h -d '2020-01-01' test_suite/mini-*`. Future `build_fixtures.sh` runs maintain the invariant automatically.

Empirically verified: `bash tests/rules/test_evm.sh` exits 0 and produces `EVM.gff3` with 4 EVM gene rows.

## Next actions (ordered)

The mtime fix unmasked three pre-existing blockers downstream of it. Two of the three integration tests now bail past the mtime hurdle and hit a different wall:

1. **Augustus species `fly` missing under Augustus 3.5 inside `containers/jamg.sif`.** `tests/e2e/test_full_dag.sh` reached `augustus_preflight` (slurm jobid 10056) and failed with `FATAL: species 'fly' not in Augustus 3.5 library`. The species set bundled with apt Augustus 3.5 differs from Augustus 3.4 (which the rule was originally written against), and `fly` is not in the default library now. Either (a) install the missing `fly` species pack in `containers/jamg.def`'s `%post`, or (b) set `augustus.optimise: true` in `test_suite/mini-config.yaml` to train one. Out of scope for the mtime commit.
2. **PASA `.cln` symlink missing in `workflow/rules/ogs.smk:91-95`** (and likely lines 156-160 for `ogs_pasa_compare_2`). `tests/rules/test_ogs.sh` reached `ogs_pasa_compare_1` and failed with `ERROR: I cannot locate the .cln file generated by seqclean, expecting transcripts.fasta.cln`. The shell command symlinks `alignAssembly.config`, `transcripts.fasta.clean`, `transcripts.fasta`, `tdn.accs` into `test_suite/output/ogs/` but not `transcripts.fasta.cln`. One-line fix: add `"ln -sf {params.cln_abs} transcripts.fasta.cln && "` and a matching `cln = ...` entry in `params:`. Out of scope for the mtime commit.
3. **Full-DAG gate test is inconclusive.** `tests/e2e/test_full_dag.sh` failed at `augustus_preflight` BEFORE reaching `genemark`. So whether `ancient()` actually suppresses the in-session "Input files updated by another job" cascade in the full-DAG scenario (which HANDOVER originally claimed was broken) remains empirically unverified. Snakemake DID plan `genemark` with `count: 1` in the DAG (possibly due to wiped `.snakemake/` provenance metadata, separately from `ancient()`). Re-run the gate once (1) is fixed: pass through genemark + into PASA = `ancient()` confirmed working. Fail at genemark with cascade = restore `touch -d '2038-01-15'` in `regen_snapshot.sh` AND `test_full_dag.sh` only, and write a follow-up plan for `--consider-ancient` plumbing through `src/jamg/cli.py`.
4. **Snapshot refresh.** After (1), (2) and (3) pass, run `tools/snapshot_dag_outputs.sh` to refresh `test_suite/fixtures/snapshot/` with the now-working `evm/` subdir (the existing snapshot was built before EVM passed and lacks `evm/`, which blocked `test_ogs.sh` until I manually copied `test_suite/output/evm/` into the snapshot to unblock the regression run; that copy is a stopgap, the snapshot needs to be regenerated cleanly).
5. **Push PASApipeline `fix/gcc14-implicit-int` to origin (alpapan/PASApipeline).** Local commit `dc5497d` is the GCC 14 implicit-int fix for `pasa-plugins/seqclean/psx/psx.c`. Requires explicit per-invocation user approval before `git push origin fix/gcc14-implicit-int` or `git push origin master` (after ff-merge).

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
- **No `--mtime-only` flag.** Removed from `src/jamg/cli.py` and all callers. snakemake's native `--rerun-triggers mtime` was what `--mtime-only` translated to. The wrapper was misleading.
- **Source fixtures are the temporal floor.** `tools/build_fixtures.sh` pins all `test_suite/mini-*` files to `2020-01-01` on every rebuild. Every derived artifact (snapshot, staged outputs, fresh rule outputs) is naturally newer. `tools/stage_fixtures.sh` is now `cp -a` only (the snapshot mtimes are preserved); no touch needed inside it. The `touch -d '2038-01-15'` on pre-staged genemark fixtures has been removed from test_evm.sh, test_ogs.sh, test_full_dag.sh, and regen_snapshot.sh because that pattern fired snakemake's clock-skew detector inside `evm_tag_genemark`.

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
