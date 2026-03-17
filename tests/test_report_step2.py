"""Tests for workflow/report_step2.py"""

import json
import os
import sys
import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import report_step2 as r2


# ── get_species_prefix ────────────────────────────────────────────────────────

def test_prefix_underscore():
    assert r2.get_species_prefix("Mmus_ENSMUSP001") == "Mmus"

def test_prefix_dot():
    assert r2.get_species_prefix("Nvec.NVE12345") == "Nvec"

def test_prefix_no_sep():
    assert r2.get_species_prefix("Mmus") == "Mmus"


# ── load_og_csv ───────────────────────────────────────────────────────────────

def test_load_og_csv_basic(tmp_path):
    csv = tmp_path / "hg1.ortholog_groups.csv"
    csv.write_text(
        "Mmus_g1\tOG_0001\tOG_0001\tRef_gene\n"
        "Hsap_g1\tOG_0001\tOG_0001\tRef_gene\n"
        "Nvec_g2\tOG_0002\tOG_0002\tRef_gene2\n"
    )
    ogs = r2.load_og_csv(csv)
    assert set(ogs.keys()) == {"OG_0001", "OG_0002"}
    assert sorted(ogs["OG_0001"]) == ["Hsap_g1", "Mmus_g1"]
    assert ogs["OG_0002"] == ["Nvec_g2"]

def test_load_og_csv_missing(tmp_path):
    assert r2.load_og_csv(tmp_path / "nonexistent.csv") == {}

def test_load_og_csv_skips_comments(tmp_path):
    csv = tmp_path / "hg.csv"
    csv.write_text("# comment\nMmus_g1\tOG_0001\tOG_0001\t-\n")
    ogs = r2.load_og_csv(csv)
    assert "OG_0001" in ogs


# ── gene_tree_to_dict ─────────────────────────────────────────────────────────

def test_gene_tree_to_dict_leaf():
    try:
        from ete3 import Tree
    except ImportError:
        pytest.skip("ete3 not installed")
    t = Tree("(Mmus_g1:0.1);")
    d = r2.gene_tree_to_dict(t.get_leaves()[0])
    assert d["leaf"] is True
    assert d["name"] == "Mmus_g1"
    assert d["species"] == "Mmus"

def test_gene_tree_to_dict_internal_og():
    try:
        from ete3 import Tree
    except ImportError:
        pytest.skip("ete3 not installed")
    t = Tree("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);")
    d = r2.gene_tree_to_dict(t)
    # root has one child (the OG node wrapping two leaves)
    assert "children" in d
    og_node = d["children"][0]
    assert og_node["name"] == "OG_0001"
    assert len(og_node["children"]) == 2
    leaves = og_node["children"]
    assert all(c["leaf"] for c in leaves)
    species = {c["species"] for c in leaves}
    assert species == {"Mmus", "Hsap"}


# ── load_possvm_trees ─────────────────────────────────────────────────────────

def test_load_possvm_trees_basic(tmp_path):
    try:
        from ete3 import Tree
    except ImportError:
        pytest.skip("ete3 not installed")

    nwk = tmp_path / "Pref.TF.HG00001.treefile.ortholog_groups.newick"
    nwk.write_text("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2,Nvec_g2:0.3);\n")

    csv = tmp_path / "Pref.TF.HG00001.treefile.ortholog_groups.csv"
    csv.write_text("Mmus_g1\tOG_0001\tOG_0001\t-\nHsap_g1\tOG_0001\tOG_0001\t-\n")

    trees, all_species = r2.load_possvm_trees(tmp_path)
    assert len(trees) == 1
    rec = trees[0]
    assert rec["hg"] == "HG00001"
    assert rec["family"] == "TF"
    assert rec["n_leaves"] == 3
    assert sorted(rec["species"]) == ["Hsap", "Mmus", "Nvec"]
    assert "OG_0001" in rec["ogs"]
    assert sorted(all_species) == ["Hsap", "Mmus", "Nvec"]

def test_load_possvm_trees_empty_dir(tmp_path):
    trees, species = r2.load_possvm_trees(tmp_path)
    assert trees == []
    assert species == []


# ── HTML generation ───────────────────────────────────────────────────────────

def test_html_output(tmp_path):
    try:
        from ete3 import Tree
    except ImportError:
        pytest.skip("ete3 not installed")

    nwk = tmp_path / "Pre.FAM.HG00001.treefile.ortholog_groups.newick"
    nwk.write_text("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);")

    out = tmp_path / "report.html"
    r2.main(["--possvm_dir", str(tmp_path), "--output", str(out)])

    html = out.read_text()
    assert "Gene Tree Explorer" in html
    assert "Mmus_g1" in html
    assert "OG_0001" in html
    # Valid JSON embedded
    start = html.index("const TREES")
    chunk = html[start: start + 2000]
    json_start = chunk.index("[")
    # should not raise
    json.loads(chunk[json_start: chunk.index(";", json_start)])
