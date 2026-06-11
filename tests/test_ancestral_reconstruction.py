"""Tests for workflow/ancestral_reconstruction.py

PastML is an external process, so these tests cover:
  - argument validation and the empty-PAM fast-exit path
  - output parsing helpers (parse_marginal_file, get_root_prob, classify_support)
  - the tree-info helper (get_tree_info)

Full end-to-end PastML execution is not exercised here because it requires
the pastml binary and a non-trivial phylogeny; that belongs in an integration
test run against real data.
"""

import os
import sys
import math
import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
pytest.importorskip("ete3")
import ancestral_reconstruction as ar


SIMPLE_TREE = "((Hsap:1,Mmus:1)Mammalia:1,(Dmel:1,Cele:1)Ecdysozoa:1)Bilateria;"


# ── classify_support ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("p, expected", [
    (1.00, "present"),
    (0.95, "present"),
    (0.90, "present"),
    (0.89, "likely_present"),
    (0.50, "likely_present"),
    (0.49, "likely_absent"),
    (0.10, "likely_absent"),
    (0.09, "absent"),
    (0.00, "absent"),
])
def test_classify_support_thresholds(p, expected):
    assert ar.classify_support(p) == expected


def test_classify_support_nan():
    assert ar.classify_support(float("nan")) == "uncertain"


# ── get_tree_info ─────────────────────────────────────────────────────────────

def test_get_tree_info_root_name(tmp_path):
    t = tmp_path / "tree.nwk"
    t.write_text(SIMPLE_TREE)
    root_name, internal, leaves = ar.get_tree_info(str(t))
    assert root_name == "Bilateria"


def test_get_tree_info_internal_nodes(tmp_path):
    t = tmp_path / "tree.nwk"
    t.write_text(SIMPLE_TREE)
    _, internal, _ = ar.get_tree_info(str(t))
    assert "Mammalia"  in internal
    assert "Ecdysozoa" in internal
    assert "Bilateria" in internal


def test_get_tree_info_leaves(tmp_path):
    t = tmp_path / "tree.nwk"
    t.write_text(SIMPLE_TREE)
    _, _, leaves = ar.get_tree_info(str(t))
    assert leaves == {"Hsap", "Mmus", "Dmel", "Cele"}


def test_get_tree_info_unnamed_root(tmp_path):
    t = tmp_path / "tree.nwk"
    t.write_text("(A:1,B:1);")   # no root name
    root_name, _, _ = ar.get_tree_info(str(t))
    assert root_name is None


# ── parse_marginal_file ───────────────────────────────────────────────────────

def test_parse_marginal_file_basic(tmp_path):
    f = tmp_path / "marginal_probabilities.og0.tab"
    f.write_text("node\t0\t1\nBilateria\t0.08\t0.92\nMammalia\t0.01\t0.99\n")
    s = ar.parse_marginal_file(f)
    assert s["Bilateria"] == pytest.approx(0.92)
    assert s["Mammalia"]  == pytest.approx(0.99)


def test_parse_marginal_file_missing(tmp_path):
    s = ar.parse_marginal_file(tmp_path / "nonexistent.tab")
    assert s.empty


def test_parse_marginal_file_empty(tmp_path):
    f = tmp_path / "p.tab"
    f.write_text("")
    s = ar.parse_marginal_file(f)
    assert s.empty


def test_parse_marginal_file_no_state1_column(tmp_path):
    f = tmp_path / "p.tab"
    f.write_text("node\t0\nBilateria\t1.0\n")
    s = ar.parse_marginal_file(f)
    assert s.empty


# ── get_root_prob ─────────────────────────────────────────────────────────────

def test_get_root_prob_by_name():
    probs = pd.Series({"Bilateria": 0.92, "Mammalia": 0.99})
    p = ar.get_root_prob(probs, "Bilateria", {"Bilateria", "Mammalia"},
                         {"Hsap", "Mmus"}, "OG1")
    assert p == pytest.approx(0.92)


def test_get_root_prob_fallback_to_non_leaf():
    """Root name not in index — fall back to non-leaf nodes."""
    probs = pd.Series({"SomeInternal": 0.75, "Hsap": 1.0, "Mmus": 0.0})
    p = ar.get_root_prob(probs, "MissingRoot", {"SomeInternal"},
                         {"Hsap", "Mmus"}, "OG1")
    assert p == pytest.approx(0.75)


def test_get_root_prob_empty_series():
    p = ar.get_root_prob(pd.Series(dtype=float), "Root", set(), set(), "OG1")
    assert math.isnan(p)


def test_get_root_prob_all_leaves():
    """If every entry in probs is a leaf, return nan."""
    probs = pd.Series({"Hsap": 1.0, "Mmus": 0.0})
    p = ar.get_root_prob(probs, None, set(), {"Hsap", "Mmus"}, "OG1")
    assert math.isnan(p)


# ── main() empty-PAM fast-exit ────────────────────────────────────────────────

def test_main_empty_pam_writes_empty_outputs(tmp_path):
    tree_f = tmp_path / "tree.nwk"
    tree_f.write_text(SIMPLE_TREE)

    pam_f = tmp_path / "pam.tsv"
    pam_f.write_text("species\n")    # header only, no OG columns

    out_f   = tmp_path / "states.tsv"
    nodes_f = tmp_path / "nodes.tsv"

    with patch("sys.argv", [
        "ar",
        "--pam",        str(pam_f),
        "--tree",       str(tree_f),
        "--output",     str(out_f),
        "--node_probs", str(nodes_f),
        "--node",       "Bilateria",
    ]):
        ar.main()

    states = pd.read_csv(out_f,   sep="\t")
    nodes  = pd.read_csv(nodes_f, sep="\t")
    assert states.empty
    assert nodes.empty
