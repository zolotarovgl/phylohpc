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
    pytest.importorskip("ete3")
    from ete3 import Tree
    t = Tree("(Mmus_g1:0.1);")
    d = r2.gene_tree_to_dict(t.get_leaves()[0])
    assert d["leaf"] is True
    assert d["name"] == "Mmus_g1"
    assert d["species"] == "Mmus"

def test_gene_tree_to_dict_internal_og():
    pytest.importorskip("ete3")
    from ete3 import Tree
    t = Tree("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);")
    d = r2.gene_tree_to_dict(t)
    og_node = d["children"][0]
    assert og_node["name"] == "OG_0001"
    assert len(og_node["children"]) == 2
    assert all(c["leaf"] for c in og_node["children"])
    assert {c["species"] for c in og_node["children"]} == {"Mmus", "Hsap"}


# ── load_possvm_trees ─────────────────────────────────────────────────────────

def test_load_possvm_trees_basic(tmp_path):
    pytest.importorskip("ete3")
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
    assert rec["n_ogs"] == 1
    assert "OG_0001" in rec["og_names"]
    # Heavy fields present
    assert "tree_dict" in rec
    assert "ogs" in rec
    assert sorted(all_species) == ["Hsap", "Mmus", "Nvec"]

def test_load_possvm_trees_empty_dir(tmp_path):
    trees, species = r2.load_possvm_trees(tmp_path)
    assert trees == []
    assert species == []


# ── HTML generation ───────────────────────────────────────────────────────────

def test_html_output_structure(tmp_path):
    """TREE_INDEX (not TREES) and per-HG lazy script tags must be present."""
    pytest.importorskip("ete3")
    nwk = tmp_path / "Pre.TF.HG00001.treefile.ortholog_groups.newick"
    nwk.write_text("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);")

    out = tmp_path / "report.html"
    r2.main(["--possvm_dir", str(tmp_path), "--output", str(out)])
    html = out.read_text()

    assert "Gene Tree Explorer" in html
    # New lazy-loading structure
    assert "TREE_INDEX" in html
    assert "TREES" not in html or "TREE_INDEX" in html  # TREE_INDEX must be there
    assert 'id="treedata-' in html
    # Leaf names appear only in lazy tag, not in index
    assert "Mmus_g1" in html
    assert "OG_0001" in html

    # TREE_INDEX is valid JSON
    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert len(parsed) == 1
    assert parsed[0]["hg"] == "HG00001"
    assert "tree_dict" not in parsed[0]   # heavy data must NOT be in the index
    assert "ogs" not in parsed[0]


def test_html_family_grouping(tmp_path):
    """Two HGs from different families both appear in TREE_INDEX."""
    pytest.importorskip("ete3")
    (tmp_path / "Pre.TF.HG00001.treefile.ortholog_groups.newick").write_text(
        "((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);"
    )
    (tmp_path / "Pre.RBP.HG00002.treefile.ortholog_groups.newick").write_text(
        "((Mmus_g2:0.1,Hsap_g2:0.1)OG_0002:0.2);"
    )

    out = tmp_path / "report.html"
    r2.main(["--possvm_dir", str(tmp_path), "--output", str(out)])
    html = out.read_text()

    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 5000]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])

    families = {r["family"] for r in parsed}
    assert families == {"TF", "RBP"}


def test_lazy_script_tags(tmp_path):
    """Each HG must have exactly one treedata- script tag in the HTML."""
    pytest.importorskip("ete3")
    for hg, fam in [("HG00001", "TF"), ("HG00002", "RBP")]:
        (tmp_path / f"Pre.{fam}.{hg}.treefile.ortholog_groups.newick").write_text(
            "((Mmus_g1:0.1,Hsap_g1:0.1)OG_X:0.2);"
        )

    out = tmp_path / "report.html"
    r2.main(["--possvm_dir", str(tmp_path), "--output", str(out)])
    html = out.read_text()

    assert html.count('id="treedata-') == 2
    # Content of each tag should be valid JSON with "tree" and "ogs" keys
    import re
    for m in re.finditer(r'<script type="application/json" id="treedata-[^"]+">([^<]+)</script>', html):
        detail = json.loads(m.group(1))
        assert "tree" in detail
        assert "ogs" in detail
