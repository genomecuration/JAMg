"""Tests for bin/trim_overlap_gff3.py.

Script CLI: -f1 <query> -f2 <ref> -o <out>. From file1, drops gene records
that overlap (same seqid, overlapping span, same strand) any gene in file2.
Writes survivors (gene + its child lines) to output.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "trim_overlap_gff3.py"


def write_gff(path: Path, *records: str) -> None:
    """Write a GFF3 file from a list of gene-records (each record is one or more
    GFF3 lines joined with '\n'). Records are separated by an empty line."""
    body = "##gff-version 3\n" + "\n\n".join(records) + "\n"
    path.write_text(body)


def run(file1: Path, file2: Path, out: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "-f1", str(file1), "-f2", str(file2), "-o", str(out)],
        capture_output=True,
        text=True,
        check=False,
    )


def gene_count(out: Path) -> int:
    return sum(1 for ln in out.read_text().splitlines() if "\tgene\t" in ln)


def test_case1_same_strand_overlap_excluded(tmp_path: Path) -> None:
    """Case 1: same-strand overlap -> file1 gene excluded."""
    f1 = tmp_path / "f1.gff"
    f2 = tmp_path / "f2.gff"
    out = tmp_path / "out.gff"
    write_gff(
        f1,
        "ctgA\tmaker\tgene\t100\t300\t.\t+\t.\tID=gQ\n"
        "ctgA\tmaker\tmRNA\t100\t300\t.\t+\t.\tID=mQ;Parent=gQ",
    )
    write_gff(
        f2,
        "ctgA\tmaker\tgene\t200\t400\t.\t+\t.\tID=gR",
    )
    proc = run(f1, f2, out)
    assert proc.returncode == 0, proc.stderr
    assert gene_count(out) == 0, f"expected 0 surviving genes, got: {out.read_text()}"


def test_case2_opposite_strand_kept(tmp_path: Path) -> None:
    """Case 2: spatial overlap but opposite strand -> file1 gene kept."""
    f1 = tmp_path / "f1.gff"
    f2 = tmp_path / "f2.gff"
    out = tmp_path / "out.gff"
    write_gff(
        f1,
        "ctgA\tmaker\tgene\t100\t300\t.\t+\t.\tID=gQ\n"
        "ctgA\tmaker\tmRNA\t100\t300\t.\t+\t.\tID=mQ;Parent=gQ",
    )
    write_gff(
        f2,
        "ctgA\tmaker\tgene\t200\t400\t.\t-\t.\tID=gR",
    )
    proc = run(f1, f2, out)
    assert proc.returncode == 0, proc.stderr
    assert gene_count(out) == 1, f"expected 1 surviving gene, got: {out.read_text()}"


def test_case3_different_seqid_kept(tmp_path: Path) -> None:
    """Case 3: overlap on coords but different seqid -> file1 gene kept."""
    f1 = tmp_path / "f1.gff"
    f2 = tmp_path / "f2.gff"
    out = tmp_path / "out.gff"
    write_gff(
        f1,
        "ctgA\tmaker\tgene\t100\t300\t.\t+\t.\tID=gQ",
    )
    write_gff(
        f2,
        "ctgB\tmaker\tgene\t100\t300\t.\t+\t.\tID=gR",
    )
    proc = run(f1, f2, out)
    assert proc.returncode == 0, proc.stderr
    assert gene_count(out) == 1


def test_case4_no_overlap_kept(tmp_path: Path) -> None:
    """Case 4: disjoint ranges -> file1 gene kept."""
    f1 = tmp_path / "f1.gff"
    f2 = tmp_path / "f2.gff"
    out = tmp_path / "out.gff"
    write_gff(
        f1,
        "ctgA\tmaker\tgene\t100\t300\t.\t+\t.\tID=gQ",
    )
    write_gff(
        f2,
        "ctgA\tmaker\tgene\t500\t800\t.\t+\t.\tID=gR",
    )
    proc = run(f1, f2, out)
    assert proc.returncode == 0, proc.stderr
    assert gene_count(out) == 1
