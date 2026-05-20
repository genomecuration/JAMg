"""Tests for bin/blast2gff.py.

Script invokes argparse with required -i / -o. It does NOT filter by
percent identity (no -p flag); column 3 is read but used only as the
GFF score field. Strand is determined by sstart<send (`+`) or sstart>send (`-`).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "blast2gff.py"


def run_script(input_path: Path, output_path: Path) -> subprocess.CompletedProcess:
    """Invoke the script as a subprocess (script is the unit; we test its behavior)."""
    return subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(input_path), "-o", str(output_path)],
        capture_output=True,
        text=True,
        check=False,
    )


def test_case1_three_hits_produce_three_gff_lines(tmp_path: Path) -> None:
    """Case 1: synthetic 3-row BLAST input -> 3 GFF lines with BLAST source, exon type."""
    blast = tmp_path / "blast.tsv"
    blast.write_text(
        "queryA\thit1\t95.5\t120\t5\t0\t1\t120\t1000\t1120\t1e-50\t200\n"
        "queryB\thit2\t99.0\t200\t2\t0\t1\t200\t5000\t5200\t1e-99\t400\n"
        "queryC\thit3\t88.0\t100\t12\t0\t1\t100\t3000\t2901\t1e-30\t150\n"
    )
    out = tmp_path / "out.gff3"
    proc = run_script(blast, out)
    assert proc.returncode == 0, f"script failed: {proc.stderr}"
    lines = [ln for ln in out.read_text().splitlines() if ln.strip()]
    assert len(lines) == 3, "expected 3 GFF lines"
    for ln in lines:
        cols = ln.split("\t")
        assert len(cols) == 9, f"GFF must have 9 tabs, got: {ln!r}"
        assert cols[1] == "BLAST"
        assert cols[2] == "exon"


def test_case2_strand_inference_from_sstart_send(tmp_path: Path) -> None:
    """Case 2: sstart<send => '+' strand; sstart>send => '-' strand and coords ordered."""
    blast = tmp_path / "blast.tsv"
    # row1: sstart 1000 < send 1120 => '+'; row2: sstart 3000 > send 2901 => '-'
    blast.write_text(
        "qP\thitP\t95.5\t120\t5\t0\t1\t120\t1000\t1120\t1e-50\t200\n"
        "qM\thitM\t88.0\t100\t12\t0\t1\t100\t3000\t2901\t1e-30\t150\n"
    )
    out = tmp_path / "out.gff3"
    proc = run_script(blast, out)
    assert proc.returncode == 0
    lines = out.read_text().splitlines()
    plus_row = next(ln for ln in lines if ln.startswith("qP"))
    minus_row = next(ln for ln in lines if ln.startswith("qM"))
    p = plus_row.split("\t")
    m = minus_row.split("\t")
    assert p[6] == "+", f"plus row strand: {plus_row}"
    assert p[3] == "1000" and p[4] == "1120"
    assert m[6] == "-", f"minus row strand: {minus_row}"
    assert int(m[3]) < int(m[4]), "minus row start<end after reorder"
    assert m[3] == "2901" and m[4] == "3000"


def test_case3_empty_input_produces_empty_output(tmp_path: Path) -> None:
    """Case 3: empty input file -> empty output file, exit 0."""
    blast = tmp_path / "empty.tsv"
    blast.write_text("")
    out = tmp_path / "out.gff3"
    proc = run_script(blast, out)
    assert proc.returncode == 0, f"empty input should exit 0; stderr={proc.stderr}"
    assert out.exists()
    assert out.read_text() == "", "empty input -> empty output"
