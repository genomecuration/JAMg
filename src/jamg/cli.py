"""JAMg v2 CLI. Thin wrapper around Snakemake + apptainer."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import click
import yaml

from . import __version__

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
SIF = REPO / "containers" / "jamg.sif"


def _detect_runtime() -> str:
    for tool in ("apptainer", "singularity"):
        if shutil.which(tool):
            return tool
    # When jamg runs INSIDE the SIF (apptainer exec ... jamg ...), the host's
    # apptainer is not on the in-container PATH. Detect that case via the
    # container-runtime env vars apptainer/singularity export to %environment.
    if os.environ.get("APPTAINER_CONTAINER") or os.environ.get("APPTAINER_NAME"):
        return "apptainer"
    if os.environ.get("SINGULARITY_CONTAINER") or os.environ.get("SINGULARITY_NAME"):
        return "singularity"
    raise click.ClickException("Neither apptainer nor singularity found on PATH.")


def collect_binds(config_dict: dict[str, Any]) -> list[str]:
    """Walk the config and return sorted bind-mount paths for apptainer.

    Every external input or output the workflow reads or writes must have
    its parent dir bind-mounted into the SIF, otherwise rules fail with
    "No such file or directory" inside the container.

    Inputs that already exist: the file's parent (or the dir itself) is
    bound. Inputs that don't exist are silently dropped (they may be
    optional, may be schema-required but not yet staged for dry-run).
    Output paths (currently just `outdir`) must be handled differently:
    the directory itself typically does not exist on first run, but the
    workflow needs to write into it, so the nearest existing ancestor is
    bound instead.
    """
    scalar_input_keys = [
        ("genome",),
        ("rnaseq", "bam"),
        ("repeats", "rna_lib"),
        ("repeats", "species_lib"),
        ("proteins", "swissprot_db"),
        ("pasa", "config_template"),
        ("genemark", "path"),
        ("genemark", "key"),
        ("trinity_denovo",),
        ("longreads",),
        ("manual_curations",),
        ("evm", "weights_file"),
        ("augustus", "extrinsic_cfg"),
        ("augustus", "metaparameters"),
    ]
    output_keys = [("outdir",)]
    binds: set[str] = {"/dev/shm", str(SIF.parent.resolve()), str(WORKFLOW.resolve())}

    def _walk(node: Any, path: tuple[str, ...]) -> Any:
        for k in path:
            if not isinstance(node, dict) or k not in node:
                return None
            node = node[k]
        return node

    def _add_input(value: Any) -> None:
        if not value:
            return
        p = Path(os.path.expanduser(str(value))).resolve()
        if p.exists():
            binds.add(str(p if p.is_dir() else p.parent))

    def _add_output(value: Any) -> None:
        """Bind the nearest existing ancestor of an output path.

        Output paths (outdir) typically do not exist on first run, so
        `Path.exists()` would drop them. Walk up until an ancestor exists
        (root is the ultimate guarantee) and bind that, so the container
        can mkdir the output tree inside it.
        """
        if not value:
            return
        p = Path(os.path.expanduser(str(value))).resolve()
        ancestor = p if p.is_dir() else p.parent
        while not ancestor.exists() and ancestor != ancestor.parent:
            ancestor = ancestor.parent
        binds.add(str(ancestor))

    for path in scalar_input_keys:
        _add_input(_walk(config_dict, path))

    for path in output_keys:
        _add_output(_walk(config_dict, path))

    fastq_pairs = _walk(config_dict, ("rnaseq", "fastq_pairs")) or []
    for pair in fastq_pairs:
        if isinstance(pair, (list, tuple)):
            for fq in pair:
                _add_input(fq)
        else:
            _add_input(pair)

    extra_gff = _walk(config_dict, ("evm", "extra_gff")) or []
    for item in extra_gff:
        if isinstance(item, dict):
            _add_input(item.get("path"))
        else:
            _add_input(item)

    return sorted(binds)


@click.group()
@click.version_option(__version__)
def cli() -> None:
    pass


@cli.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--dry-run", is_flag=True)
@click.option("--unlock", is_flag=True, help="Release a stale Snakemake lock before running.")
@click.option("--cores", default=1, type=int, help="local: per-job CPU; slurm: ignored.")
@click.option(
    "--jobs",
    default=None,
    type=int,
    help="slurm: concurrent submissions cap. Required when --executor=slurm.",
)
@click.option(
    "--until",
    "until_rule",
    default=None,
    help="Stop after the named rule and its dependencies are produced.",
)
@click.option("--executor", default="local", type=click.Choice(["local", "slurm"]))
@click.option(
    "--container-runtime",
    default=None,
    type=click.Choice(["apptainer", "singularity"]),
    help="Override autodetect.",
)
def run(
    config_path: str,
    dry_run: bool,
    unlock: bool,
    cores: int,
    jobs: int | None,
    until_rule: str | None,
    executor: str,
    container_runtime: str | None,
) -> None:
    runtime = container_runtime or _detect_runtime()
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}
    binds = collect_binds(cfg)
    bind_args = " ".join(f"--bind {b}" for b in binds)

    # Detect whether we're running inside an apptainer/singularity SIF.
    # When inside, every rule's tools come from the SIF we're already in,
    # so we MUST NOT pass --software-deployment-method (which would try
    # to dispatch nested apptainer calls per rule, breaking the run).
    inside_sif = bool(
        os.environ.get("APPTAINER_CONTAINER")
        or os.environ.get("APPTAINER_NAME")
        or os.environ.get("SINGULARITY_CONTAINER")
        or os.environ.get("SINGULARITY_NAME")
    )

    cmd = [
        "snakemake",
        "--snakefile",
        str(WORKFLOW / "Snakefile"),
        "--configfile",
        config_path,
    ]
    if not inside_sif:
        # Snakemake 9: --software-deployment-method selects the runtime;
        # --apptainer-args is passed through to apptainer; the cache-prefix
        # flag is still --singularity-prefix in 9.x (apptainer rename not
        # adopted upstream).
        cmd += [
            "--software-deployment-method",
            "apptainer",
            "--singularity-prefix",
            str(REPO / ".snakemake" / "apptainer"),
            "--apptainer-args",
            bind_args,
        ]
    if executor == "slurm":
        if jobs is None:
            raise click.ClickException("--jobs N required when --executor=slurm")
        cmd += [
            "--executor",
            "slurm",
            "--profile",
            str(WORKFLOW / "profiles" / "slurm"),
            "--jobs",
            str(jobs),
        ]
    else:
        cmd += ["--cores", str(cores)]
    if dry_run:
        cmd.append("--dry-run")
    if unlock:
        cmd.append("--unlock")
    if until_rule:
        # Pass the rule name as a positional target. Snakemake's --until takes
        # nargs='+', so "--until X X" eats both tokens as --until args and
        # leaves no positional target, meaning `rule all` is still the
        # resolution root. A positional target alone gives the same DAG
        # behavior (build the rule + its dependencies) without that bug.
        cmd.append(until_rule)
    sys.exit(subprocess.call(cmd, env={**os.environ, "JAMG_CONTAINER_RUNTIME": runtime}))


if __name__ == "__main__":
    cli()
