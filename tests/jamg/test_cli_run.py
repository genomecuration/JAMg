"""Tests for `jamg run` argv construction (Phase 2.3 Click CLI).

We monkeypatch `subprocess.call` (the EXTERNAL boundary — the only mock allowed
per the v2 plan §4.10 and global CLAUDE.md test policy) to capture the argv
list that `jamg run` builds for `snakemake`, without actually invoking it.

Cases (4):
  1. --executor local --dry-run => --dry-run + --cores + no --executor slurm
  2. --executor slurm --jobs 5 => --executor slurm + --jobs 5 + --profile
  3. --executor slurm without --jobs => ClickException ('--jobs N required')
  4. --unlock => --unlock token appears in the argv

We use Click's CliRunner so SystemExit is captured cleanly.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

# Make src/ importable
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

import jamg.cli as cli_mod  # noqa: E402
from jamg.cli import cli  # noqa: E402


def _write_min_config(tmp_path: Path) -> Path:
    """Write a minimal config file the CLI can yaml.safe_load."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "genome: /nonexistent/genome.fa\n"
        "outdir: /nonexistent/out\n"
    )
    return cfg


@pytest.fixture
def capture_subprocess(monkeypatch):
    """Replace subprocess.call (third-party boundary) with a no-op capture.

    The mock returns 0 immediately; the test asserts on the captured argv,
    not on snakemake's behavior.
    """
    captured: dict[str, list] = {"argv": [], "env": {}}

    def fake_call(cmd, env=None, **kwargs):
        captured["argv"] = list(cmd)
        captured["env"] = dict(env) if env else {}
        return 0

    monkeypatch.setattr(cli_mod.subprocess, "call", fake_call)
    # The CLI also calls _detect_runtime() which probes shutil.which.
    # Force it to return 'apptainer' so the test doesn't need apptainer on PATH.
    monkeypatch.setattr(cli_mod, "_detect_runtime", lambda: "apptainer")
    # The CLI also reads APPTAINER_* env vars to detect inside-SIF mode.
    # Ensure those are NOT set (so we exercise the deployment-method branch).
    for var in ("APPTAINER_CONTAINER", "APPTAINER_NAME",
                "SINGULARITY_CONTAINER", "SINGULARITY_NAME"):
        monkeypatch.delenv(var, raising=False)
    return captured


def test_case1_local_dry_run(tmp_path, capture_subprocess):
    cfg = _write_min_config(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--config", str(cfg), "--executor", "local",
                                  "--cores", "4", "--dry-run"])
    assert result.exit_code == 0, result.output
    argv = capture_subprocess["argv"]
    assert "snakemake" in argv[0]
    assert "--cores" in argv and "4" in argv
    assert "--dry-run" in argv
    # local => no slurm executor token
    assert "slurm" not in argv


def test_case2_slurm_with_jobs(tmp_path, capture_subprocess):
    cfg = _write_min_config(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--config", str(cfg), "--executor", "slurm",
                                  "--jobs", "5"])
    assert result.exit_code == 0, result.output
    argv = capture_subprocess["argv"]
    assert "--executor" in argv
    # The token after --executor must be 'slurm'
    assert argv[argv.index("--executor") + 1] == "slurm"
    assert "--jobs" in argv and "5" in argv
    assert "--profile" in argv


def test_case3_slurm_without_jobs_raises(tmp_path, capture_subprocess):
    cfg = _write_min_config(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--config", str(cfg), "--executor", "slurm"])
    # ClickException -> exit_code != 0
    assert result.exit_code != 0, "slurm without --jobs should error"
    assert "--jobs" in (result.output + (result.stderr_bytes or b"").decode("utf8", "ignore"))


def test_case4_unlock_appears(tmp_path, capture_subprocess):
    cfg = _write_min_config(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--config", str(cfg), "--executor", "local",
                                  "--cores", "1", "--unlock"])
    assert result.exit_code == 0, result.output
    assert "--unlock" in capture_subprocess["argv"]
