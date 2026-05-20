"""Tests for bin/run_split_augustus.py::prepare_augustus_commands.

We test the pure-function `prepare_augustus_commands` directly — it builds the
augustus argv list per chunk. No subprocess; no Augustus. Inputs are simple
dicts of (chunk_id -> seq/hints file objects).
"""

from __future__ import annotations

import importlib.util
import sys
import types
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "run_split_augustus.py"


def _load_module() -> types.ModuleType:
    """Load run_split_augustus.py as a module without invoking its argparse main.

    Stub `Bio` if biopython isn't installed: `prepare_augustus_commands` does
    not touch SeqIO so a sentinel module is enough to satisfy the top-level
    import. Biopython is an external boundary (the I/O-level callers exercise
    real biopython via the SIF) and stubbing it here keeps this unit test
    runnable on the host perl/python env per the pure-function contract.
    """
    if "Bio" not in sys.modules:
        try:  # use real biopython if present (SIF perl env supplies it)
            import Bio  # noqa: F401
        except ModuleNotFoundError:
            bio_mod = types.ModuleType("Bio")
            seqio_mod = types.ModuleType("Bio.SeqIO")
            sys.modules["Bio"] = bio_mod
            sys.modules["Bio.SeqIO"] = seqio_mod
            bio_mod.SeqIO = seqio_mod

    spec = importlib.util.spec_from_file_location("run_split_augustus", str(SCRIPT))
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["run_split_augustus"] = mod
    spec.loader.exec_module(mod)
    return mod


@dataclass
class FakeFile:
    """Stub for an argparse FileType-opened file. Augustus cmd uses only `.name`."""

    name: str


def test_case1_basic_cmd_construction(tmp_path: Path) -> None:
    """Case 1: minimum-args invocation produces a list of cmds, one per chunk,
    each beginning with 'augustus' and including --species + --uniqueGeneId."""
    mod = _load_module()
    rundir = tmp_path
    seq_files = {0: FakeFile(name=str(tmp_path / "genome.fa.00"))}
    hint_files: dict[int, FakeFile] = {}  # no hints
    cmds = mod.prepare_augustus_commands(
        UTR="off",
        gff3="on",
        species="drosophila",
        uniqueGeneId="true",
        genemodel="complete",
        alternatives="false",
        extrinsicCfgFile=None,
        hint_files=hint_files,
        seq_files=seq_files,
        rundir=str(rundir),
        softmasking="1",
    )
    assert isinstance(cmds, list) and len(cmds) == 1
    cmd = cmds[0]
    assert cmd.startswith("augustus "), f"cmd must start with 'augustus ': {cmd}"
    assert "--species=drosophila" in cmd
    assert "--uniqueGeneId=true" in cmd
    assert "--gff3=on" in cmd
    assert "--UTR=off" in cmd
    assert "--genemodel=complete" in cmd
    assert "--softmasking=1" in cmd
    # hints/extrinsic not in cmd when no hints
    assert "--hintsfile" not in cmd
    assert "--extrinsicCfgFile" not in cmd


def test_case2_utr_on_propagates(tmp_path: Path) -> None:
    """Case 2: UTR='on' propagates into every chunk's cmd."""
    mod = _load_module()
    seq_files = {0: FakeFile(name=str(tmp_path / "a.fa.00")),
                 1: FakeFile(name=str(tmp_path / "a.fa.01"))}
    cmds = mod.prepare_augustus_commands(
        UTR="on",
        gff3="on",
        species="dmel",
        uniqueGeneId="true",
        genemodel="complete",
        alternatives="false",
        extrinsicCfgFile=None,
        hint_files={},
        seq_files=seq_files,
        rundir=str(tmp_path),
        softmasking="1",
    )
    assert len(cmds) == 2
    for cmd in cmds:
        assert "--UTR=on" in cmd, f"--UTR=on missing in: {cmd}"


def test_case3_uniqueGeneId_false_and_hints_chunk(tmp_path: Path) -> None:
    """Case 3: uniqueGeneId='false' propagates; hint file appears for the
    chunks that have hints (and only for those)."""
    mod = _load_module()
    seq_files = {0: FakeFile(name=str(tmp_path / "g.fa.00")),
                 1: FakeFile(name=str(tmp_path / "g.fa.01"))}
    # chunk 0 has hints; chunk 1 does not
    hint_files = {0: FakeFile(name=str(tmp_path / "h.gff.00"))}
    cmds = mod.prepare_augustus_commands(
        UTR="off",
        gff3="on",
        species="dmel",
        uniqueGeneId="false",
        genemodel="partial",
        alternatives="true",
        extrinsicCfgFile="/some/extrinsic.cfg",
        hint_files=hint_files,
        seq_files=seq_files,
        rundir=str(tmp_path),
        softmasking="1",
    )
    assert len(cmds) == 2
    assert "--uniqueGeneId=false" in cmds[0]
    assert "--hintsfile=" + str(tmp_path / "h.gff.00") in cmds[0]
    assert "--extrinsicCfgFile=/some/extrinsic.cfg" in cmds[0]
    # chunk 1: no hints
    assert "--hintsfile" not in cmds[1]
    assert "--extrinsicCfgFile" not in cmds[1]
