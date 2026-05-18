"""Tests for jamg.cli.collect_binds.

Eight cases enumerated in the Phase 2 plan (Step 2.3.6). All use pytest's
`tmp_path` fixture for real-path testing; no `unittest.mock`, no patched
filesystem. The function under test is the only externally-visible Python
function the workflow depends on for path-binding correctness.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from jamg.cli import SIF, WORKFLOW, collect_binds

DEFAULTS = {"/dev/shm", str(SIF.parent.resolve()), str(WORKFLOW.resolve())}


def _full_cfg(tmp_path: Path) -> dict:
    """Materialize every scalar key the function walks to a real tmp file."""
    files = {}
    for name in (
        "genome.fa",
        "rnaseq.bam",
        "rna.lib",
        "species.lib",
        "swissprot.fasta",
        "pasa.config",
        "gmes.tar.gz",
        "gm_key.gz",
        "trinity.fasta",
        "longreads.fasta",
        "manual.gff3",
        "weights.txt",
        "extrinsic.cfg",
        "metaparams.cfg",
    ):
        p = tmp_path / name
        p.write_text("")
        files[name] = p
    return {
        "genome": str(files["genome.fa"]),
        "outdir": str(tmp_path),
        "rnaseq": {"bam": str(files["rnaseq.bam"])},
        "repeats": {
            "rna_lib": str(files["rna.lib"]),
            "species_lib": str(files["species.lib"]),
        },
        "proteins": {"swissprot_db": str(files["swissprot.fasta"])},
        "pasa": {"config_template": str(files["pasa.config"])},
        "genemark": {
            "path": str(files["gmes.tar.gz"]),
            "key": str(files["gm_key.gz"]),
        },
        "trinity_denovo": str(files["trinity.fasta"]),
        "longreads": str(files["longreads.fasta"]),
        "manual_curations": str(files["manual.gff3"]),
        "evm": {
            "weights_file": str(files["weights.txt"]),
            "extra_gff": [],
        },
        "augustus": {
            "extrinsic_cfg": str(files["extrinsic.cfg"]),
            "metaparameters": str(files["metaparams.cfg"]),
        },
    }


def test_full_config_resolves_every_scalar_key(tmp_path: Path) -> None:
    """Case 1: every scalar key plus defaults appear in the bind list."""
    cfg = _full_cfg(tmp_path)
    binds = set(collect_binds(cfg))
    assert DEFAULTS.issubset(binds)
    expected_parent = str(tmp_path.resolve())
    assert expected_parent in binds


def test_empty_config_returns_only_defaults() -> None:
    """Case 2: empty dict yields only the three default binds."""
    assert set(collect_binds({})) == DEFAULTS


def test_null_optional_path_omitted() -> None:
    """Case 3: optional key set to None contributes no bind."""
    binds = set(collect_binds({"trinity_denovo": None}))
    assert binds == DEFAULTS


def test_tilde_expansion(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Case 4: ~ paths expand via os.path.expanduser to the user's home."""
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    sw = fake_home / "software"
    sw.mkdir()
    tarball = sw / "gmes.tar.gz"
    tarball.write_text("")
    monkeypatch.setenv("HOME", str(fake_home))
    binds = set(collect_binds({"genemark": {"path": "~/software/gmes.tar.gz"}}))
    assert str(sw.resolve()) in binds


def test_non_existent_path_dropped(tmp_path: Path) -> None:
    """Case 5: input paths that don't exist are silently dropped."""
    cfg = {"genome": str(tmp_path / "absent.fa")}
    binds = set(collect_binds(cfg))
    assert binds == DEFAULTS


def test_outdir_binds_nearest_existing_ancestor(tmp_path: Path) -> None:
    """outdir is an output: even when it does not yet exist, the nearest
    existing ancestor must be bound so the container can mkdir into it."""
    nested = tmp_path / "results" / "run_001" / "stage"
    binds = set(collect_binds({"outdir": str(nested)}))
    assert str(tmp_path.resolve()) in binds


def test_duplicate_paths_deduped(tmp_path: Path) -> None:
    """Case 6: two config keys pointing into the same dir collapse to one."""
    f1 = tmp_path / "a.fa"
    f2 = tmp_path / "b.fa"
    f1.write_text("")
    f2.write_text("")
    cfg = {"genome": str(f1), "proteins": {"swissprot_db": str(f2)}}
    binds = collect_binds(cfg)
    parent = str(tmp_path.resolve())
    assert binds.count(parent) == 1
    assert parent in binds


def test_evm_extra_gff_walked(tmp_path: Path) -> None:
    """Case 7: evm.extra_gff[].path values contribute binds."""
    d1 = tmp_path / "liftoff"
    d1.mkdir()
    d2 = tmp_path / "curations"
    d2.mkdir()
    a = d1 / "A.gff3"
    b = d2 / "B.gff3"
    a.write_text("")
    b.write_text("")
    cfg = {
        "evm": {
            "extra_gff": [
                {"path": str(a), "source_tag": "LIFTOFF"},
                {"path": str(b), "source_tag": "MANUAL"},
            ]
        }
    }
    binds = set(collect_binds(cfg))
    assert str(d1.resolve()) in binds
    assert str(d2.resolve()) in binds


def test_rnaseq_fastq_pairs_walked(tmp_path: Path) -> None:
    """Case 8: rnaseq.fastq_pairs[][] walks pairs and bare strings alike.

    collect_binds runs BEFORE snakemake validates the config against the
    schema, so it must tolerate inputs the schema would later reject.
    The bare-string branch is therefore a defensive codepath, not a
    schema-supported config shape. The schema constrains items to length-
    2 string arrays; the bare-string case is exercised here to lock in
    the defensive behaviour against accidental future regressions.
    """
    p1 = tmp_path / "lane1"
    p2 = tmp_path / "lane2"
    p1.mkdir()
    p2.mkdir()
    pairs = []
    for d in (p1, p2):
        r1 = d / "R1.fq"
        r2 = d / "R2.fq"
        r1.write_text("")
        r2.write_text("")
        pairs.append([str(r1), str(r2)])
    binds = set(collect_binds({"rnaseq": {"fastq_pairs": pairs}}))
    assert str(p1.resolve()) in binds
    assert str(p2.resolve()) in binds

    single = tmp_path / "single.fq"
    single.write_text("")
    binds_single = set(collect_binds({"rnaseq": {"fastq_pairs": [str(single)]}}))
    assert str(tmp_path.resolve()) in binds_single
