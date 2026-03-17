"""Tests for workflow/select_hgs.py"""

import os
import sys
import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import select_hgs


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_fasta(path: Path, sequences):
    """sequences: list of (header, seq) tuples"""
    with open(path, "w") as fh:
        for hdr, seq in sequences:
            fh.write(f">{hdr}\n{seq}\n")


# ---------------------------------------------------------------------------
# analyze_fasta
# ---------------------------------------------------------------------------

def test_analyze_fasta_counts_sequences(tmp_path):
    f = tmp_path / "hg1.fasta"
    write_fasta(f, [("Hsap_g1", "ACGT"), ("Mmus_g2", "TTTT"), ("Drer_g3", "AAAA")])
    nseq, has_soi, n_sps = select_hgs.analyze_fasta(f, None)
    assert nseq == 3


def test_analyze_fasta_counts_species(tmp_path):
    f = tmp_path / "hg1.fasta"
    write_fasta(f, [("Hsap_g1", "ACGT"), ("Hsap_g2", "TTTT"), ("Mmus_g1", "AAAA")])
    _, _, n_sps = select_hgs.analyze_fasta(f, None)
    assert n_sps == 2


def test_analyze_fasta_detects_soi_present(tmp_path):
    f = tmp_path / "hg1.fasta"
    write_fasta(f, [("Hsap_g1", "ACGT"), ("Mmus_g1", "TTTT")])
    _, has_soi, _ = select_hgs.analyze_fasta(f, "Hsap")
    assert has_soi is True


def test_analyze_fasta_detects_soi_absent(tmp_path):
    f = tmp_path / "hg1.fasta"
    write_fasta(f, [("Mmus_g1", "ACGT"), ("Drer_g1", "TTTT")])
    _, has_soi, _ = select_hgs.analyze_fasta(f, "Hsap")
    assert has_soi is False


def test_analyze_fasta_no_soi_filter_always_true(tmp_path):
    f = tmp_path / "hg1.fasta"
    write_fasta(f, [("Mmus_g1", "ACGT")])
    _, has_soi, _ = select_hgs.analyze_fasta(f, None)
    assert has_soi is True


def test_analyze_fasta_empty_file(tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    nseq, _, n_sps = select_hgs.analyze_fasta(f, None)
    assert nseq == 0
    assert n_sps == 0


def test_analyze_fasta_header_without_underscore_not_counted_as_species(tmp_path):
    """Headers without '_' contribute no species name."""
    f = tmp_path / "hg1.fasta"
    write_fasta(f, [("GENE1", "ACGT"), ("GENE2", "TTTT")])
    _, _, n_sps = select_hgs.analyze_fasta(f, None)
    assert n_sps == 0


# ---------------------------------------------------------------------------
# iter_fasta_files
# ---------------------------------------------------------------------------

def test_iter_fasta_files_returns_correct_extensions(tmp_path):
    for name in ["a.fasta", "b.fa", "c.faa", "d.fna", "e.txt", "f.csv"]:
        (tmp_path / name).write_text(">seq\nACGT\n")

    found = {p.name for p in select_hgs.iter_fasta_files(tmp_path)}
    assert found == {"a.fasta", "b.fa", "c.faa", "d.fna"}
    assert "e.txt" not in found
    assert "f.csv" not in found


def test_iter_fasta_files_empty_dir(tmp_path):
    results = list(select_hgs.iter_fasta_files(tmp_path))
    assert results == []


def test_iter_fasta_files_case_insensitive_extension(tmp_path):
    (tmp_path / "X.FASTA").write_text(">seq\nACGT\n")
    found = {p.name for p in select_hgs.iter_fasta_files(tmp_path)}
    assert "X.FASTA" in found
