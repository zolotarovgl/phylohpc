"""Tests for workflow/get_seqstat.py"""

import os
import sys
import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import get_seqstat


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_fasta(path, sequences):
    """sequences: list of (header, seq) tuples; seq may be a list of lines."""
    with open(path, "w") as fh:
        for hdr, seq in sequences:
            fh.write(f">{hdr}\n")
            if isinstance(seq, list):
                for part in seq:
                    fh.write(part + "\n")
            else:
                fh.write(seq + "\n")


# ---------------------------------------------------------------------------
# read_fasta_lengths
# ---------------------------------------------------------------------------

def test_read_fasta_lengths_single_sequence(tmp_path):
    f = tmp_path / "a.fasta"
    write_fasta(f, [("s1", "ACGTACGT")])
    assert get_seqstat.read_fasta_lengths(str(f)) == [8]


def test_read_fasta_lengths_multiple_sequences(tmp_path):
    f = tmp_path / "a.fasta"
    write_fasta(f, [("s1", "ACGT"), ("s2", "TTTTTT"), ("s3", "AA")])
    result = get_seqstat.read_fasta_lengths(str(f))
    assert sorted(result) == [2, 4, 6]


def test_read_fasta_lengths_multiline_sequence(tmp_path):
    f = tmp_path / "ml.fasta"
    # Two-line sequence: ACGT + ACGT = 8 chars
    write_fasta(f, [("s1", ["ACGT", "ACGT"])])
    assert get_seqstat.read_fasta_lengths(str(f)) == [8]


def test_read_fasta_lengths_empty_file(tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    assert get_seqstat.read_fasta_lengths(str(f)) == []


def test_read_fasta_lengths_ignores_blank_lines(tmp_path):
    f = tmp_path / "blank.fasta"
    f.write_text(">s1\nACGT\n\nACGT\n")
    # blank line inside sequence body is skipped; total = 8
    result = get_seqstat.read_fasta_lengths(str(f))
    assert result == [8]


# ---------------------------------------------------------------------------
# process_file
# ---------------------------------------------------------------------------

def test_process_file_outputs_correct_stats(tmp_path, capsys):
    f = tmp_path / "mycluster.fasta"
    # lengths: 4, 6, 8 → median = 6
    write_fasta(f, [("s1", "ACGT"), ("s2", "ACGTAC"), ("s3", "ACGTACGT")])
    get_seqstat.process_file(str(f))
    out = capsys.readouterr().out.strip()
    parts = out.split("\t")
    assert parts[0] == "mycluster"
    assert parts[1] == "3"   # n sequences
    assert parts[2] == "6"   # median length


def test_process_file_empty_fasta(tmp_path, capsys):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    get_seqstat.process_file(str(f))
    out = capsys.readouterr().out.strip()
    parts = out.split("\t")
    assert parts[1] == "0"
    assert parts[2] == "0"


def test_process_file_uses_basename_without_extension(tmp_path, capsys):
    f = tmp_path / "myfile.fasta"
    write_fasta(f, [("s1", "ACGT")])
    get_seqstat.process_file(str(f))
    out = capsys.readouterr().out.strip()
    assert out.startswith("myfile\t")
