"""Tests for workflow/visualize_ancestry.py"""

import json
import os
import sys
import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
pytest.importorskip("ete3")
import visualize_ancestry as va


# ── tree_to_dict ──────────────────────────────────────────────────────────────

def make_tree(newick):
    from ete3 import Tree
    return Tree(newick, format=1)


def test_tree_to_dict_leaf():
    t = make_tree("(Hsap:1.0,Mmus:1.0)Mammalia;")
    d = va.tree_to_dict(t)
    assert d["name"] == "Mammalia"
    assert len(d["children"]) == 2
    names = {c["name"] for c in d["children"]}
    assert names == {"Hsap", "Mmus"}


def test_tree_to_dict_leaf_flag():
    t = make_tree("(Hsap:1,Mmus:1)Root;")
    d = va.tree_to_dict(t)
    for child in d["children"]:
        assert child.get("leaf") is True


def test_tree_to_dict_internal_no_leaf_flag():
    t = make_tree("((A:1,B:1)AB:1,C:1)Root;")
    d = va.tree_to_dict(t)
    assert "leaf" not in d
    # AB is an internal node
    ab = next(c for c in d["children"] if c["name"] == "AB")
    assert "leaf" not in ab


def test_tree_to_dict_branch_lengths():
    t = make_tree("(Hsap:2.5,Mmus:3.0)Root;")
    d = va.tree_to_dict(t)
    leaves = {c["name"]: c["dist"] for c in d["children"]}
    assert leaves["Hsap"] == pytest.approx(2.5)
    assert leaves["Mmus"] == pytest.approx(3.0)


def test_tree_to_dict_nested():
    t = make_tree("((Hsap:1,Mmus:1)Mammalia:1,(Dmel:1,Cele:1)Ecdysozoa:1)Bilateria;")
    d = va.tree_to_dict(t)
    assert d["name"] == "Bilateria"
    child_names = {c["name"] for c in d["children"]}
    assert child_names == {"Mammalia", "Ecdysozoa"}


# ── parse_og_name ─────────────────────────────────────────────────────────────

def test_parse_og_name_full():
    r = va.parse_og_name("tfs.HG001.OG1")
    assert r == {"prefix": "tfs", "hg": "HG001", "og_num": "OG1"}


def test_parse_og_name_two_parts():
    r = va.parse_og_name("tfs.HG001")
    assert r["prefix"] == "tfs"
    assert r["hg"] == "HG001"
    assert r["og_num"] == ""


def test_parse_og_name_single():
    r = va.parse_og_name("tfs")
    assert r["prefix"] == "tfs"
    assert r["hg"] == ""
    assert r["og_num"] == ""


# ── HTML generation (integration) ────────────────────────────────────────────

SIMPLE_TREE = "((Hsap:1,Mmus:1)Mammalia:1,(Dmel:1,Cele:1)Ecdysozoa:1)Bilateria;"

def write_test_inputs(tmp_path):
    # Tree
    tree_f = tmp_path / "test.pruned.tree"
    tree_f.write_text(SIMPLE_TREE)

    # node_probs.tsv
    node_probs = pd.DataFrame([
        {"og": "tfs.HG001.OG1", "node": "Bilateria",  "P_present": 0.92},
        {"og": "tfs.HG001.OG1", "node": "Mammalia",   "P_present": 0.99},
        {"og": "tfs.HG001.OG1", "node": "Ecdysozoa",  "P_present": 0.85},
        {"og": "tfs.HG001.OG2", "node": "Bilateria",  "P_present": 0.10},
        {"og": "tfs.HG001.OG2", "node": "Mammalia",   "P_present": 0.05},
        {"og": "tfs.HG001.OG2", "node": "Ecdysozoa",  "P_present": 0.15},
    ])
    np_f = tmp_path / "test.node_probs.tsv"
    node_probs.to_csv(np_f, sep="\t", index=False)

    # ancestral_states.tsv
    states = pd.DataFrame([
        {"og": "tfs.HG001.OG1", "n_present": 4, "n_total": 4, "P_at_root": 0.92, "support": "present"},
        {"og": "tfs.HG001.OG2", "n_present": 1, "n_total": 4, "P_at_root": 0.10, "support": "likely_absent"},
    ])
    st_f = tmp_path / "test.ancestral_states.tsv"
    states.to_csv(st_f, sep="\t", index=False)

    # pam.tsv
    pam = pd.DataFrame(
        {"tfs.HG001.OG1": [1, 1, 1, 1], "tfs.HG001.OG2": [1, 0, 0, 0]},
        index=["Hsap", "Mmus", "Dmel", "Cele"]
    )
    pam.index.name = "species"
    pam_f = tmp_path / "test.pam.tsv"
    pam.to_csv(pam_f, sep="\t")

    return tree_f, np_f, st_f, pam_f


def test_html_generated(tmp_path):
    tree_f, np_f, st_f, pam_f = write_test_inputs(tmp_path)
    out_f = tmp_path / "out.html"
    args = ["visualize_ancestry",
            "--tree",       str(tree_f),
            "--node_probs", str(np_f),
            "--states",     str(st_f),
            "--pam",        str(pam_f),
            "--output",     str(out_f),
            "--node",       "Bilateria"]
    with patch("sys.argv", args):
        va.main()
    assert out_f.exists()
    assert out_f.stat().st_size > 1000


def test_html_contains_node_name(tmp_path):
    tree_f, np_f, st_f, pam_f = write_test_inputs(tmp_path)
    out_f = tmp_path / "out.html"
    with patch("sys.argv", ["v", "--tree", str(tree_f), "--node_probs", str(np_f),
                             "--states", str(st_f), "--pam", str(pam_f),
                             "--output", str(out_f), "--node", "Bilateria"]):
        va.main()
    html = out_f.read_text()
    assert "Bilateria" in html


def test_html_embeds_og_data(tmp_path):
    tree_f, np_f, st_f, pam_f = write_test_inputs(tmp_path)
    out_f = tmp_path / "out.html"
    with patch("sys.argv", ["v", "--tree", str(tree_f), "--node_probs", str(np_f),
                             "--states", str(st_f), "--pam", str(pam_f),
                             "--output", str(out_f), "--node", "Bilateria"]):
        va.main()
    html = out_f.read_text()
    assert "tfs.HG001.OG1" in html
    assert "tfs.HG001.OG2" in html


def test_html_embeds_valid_json(tmp_path):
    """NODE_PROBS and OG_META placeholders must be replaced with valid JSON."""
    tree_f, np_f, st_f, pam_f = write_test_inputs(tmp_path)
    out_f = tmp_path / "out.html"
    with patch("sys.argv", ["v", "--tree", str(tree_f), "--node_probs", str(np_f),
                             "--states", str(st_f), "--pam", str(pam_f),
                             "--output", str(out_f), "--node", "Bilateria"]):
        va.main()
    html = out_f.read_text()
    # None of the placeholder tokens should remain
    assert "%%" not in html


def test_html_no_pam_still_generates(tmp_path):
    """--pam is optional; HTML should still be produced."""
    tree_f, np_f, st_f, _ = write_test_inputs(tmp_path)
    out_f = tmp_path / "out.html"
    with patch("sys.argv", ["v", "--tree", str(tree_f), "--node_probs", str(np_f),
                             "--states", str(st_f),
                             "--output", str(out_f), "--node", "Bilateria"]):
        va.main()
    assert out_f.exists()
    assert out_f.stat().st_size > 1000
