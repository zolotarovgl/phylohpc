"""Tests for workflow/extract_clade.py"""

import os
import sys
import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))

# extract_clade.py requires ete3; skip the whole module when it is absent
# (e.g. when running under the system Python rather than the phylo conda env).
pytest.importorskip("ete3")

import subprocess

PYTHON = sys.executable

# Minimal species tree with named internal nodes
TREE_NEWICK = "((Hsap:1,Mmus:1)Mammalia:1,(Dmel:1,Cele:1)Ecdysozoa:1)Bilateria;"


def run_extract(tmp_path, node_name, tree_content=TREE_NEWICK):
    tree_file = tmp_path / "test.tree"
    tree_file.write_text(tree_content)
    script = Path(__file__).parent.parent / "workflow" / "extract_clade.py"
    result = subprocess.run(
        [PYTHON, str(script),
         "--tree",       str(tree_file),
         "--node",       node_name,
         "--out_prefix", str(tmp_path / node_name)],
        capture_output=True, text=True, cwd=str(tmp_path),
    )
    return result


# ── Happy-path ────────────────────────────────────────────────────────────────

def test_bilateria_in_species(tmp_path):
    r = run_extract(tmp_path, "Bilateria")
    assert r.returncode == 0
    in_sp = (tmp_path / "Bilateria.in_species.txt").read_text().splitlines()
    assert sorted(in_sp) == ["Cele", "Dmel", "Hsap", "Mmus"]


def test_bilateria_ignore_species_empty(tmp_path):
    """Root clade — nothing to ignore."""
    r = run_extract(tmp_path, "Bilateria")
    assert r.returncode == 0
    ign = (tmp_path / "Bilateria.ignore_species.txt").read_text().strip()
    assert ign == ""


def test_mammalia_in_species(tmp_path):
    r = run_extract(tmp_path, "Mammalia")
    assert r.returncode == 0
    in_sp = (tmp_path / "Mammalia.in_species.txt").read_text().splitlines()
    assert sorted(in_sp) == ["Hsap", "Mmus"]


def test_mammalia_ignore_species(tmp_path):
    r = run_extract(tmp_path, "Mammalia")
    assert r.returncode == 0
    ign = sorted((tmp_path / "Mammalia.ignore_species.txt").read_text().splitlines())
    assert ign == ["Cele", "Dmel"]


def test_pruned_tree_written(tmp_path):
    r = run_extract(tmp_path, "Mammalia")
    assert r.returncode == 0
    pruned = tmp_path / "Mammalia.pruned.tree"
    assert pruned.exists()
    content = pruned.read_text()
    assert "Hsap" in content
    assert "Mmus" in content
    # Out-of-clade species must not appear in the pruned tree
    assert "Dmel" not in content
    assert "Cele" not in content


def test_ecdysozoa_correct_split(tmp_path):
    r = run_extract(tmp_path, "Ecdysozoa")
    assert r.returncode == 0
    in_sp  = sorted((tmp_path / "Ecdysozoa.in_species.txt").read_text().splitlines())
    ign_sp = sorted((tmp_path / "Ecdysozoa.ignore_species.txt").read_text().splitlines())
    assert in_sp  == ["Cele", "Dmel"]
    assert ign_sp == ["Hsap", "Mmus"]


# ── Error handling ─────────────────────────────────────────────────────────────

def test_unknown_node_exits_nonzero(tmp_path):
    r = run_extract(tmp_path, "Nonexistent")
    assert r.returncode != 0
    assert "not found" in r.stderr.lower() or "error" in r.stderr.lower()


def test_single_species_clade_exits_nonzero(tmp_path):
    """A clade with only one leaf cannot be meaningfully analysed."""
    tree = "((A:1)SingleChild:1,B:1)Root;"
    r = run_extract(tmp_path, "SingleChild", tree_content=tree)
    assert r.returncode != 0
