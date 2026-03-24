"""Tests for workflow/report_step2.py"""

import json
import os
import re
import sys
import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import report_step2 as r2


def write_family_info(tmp_path, *families):
    path = tmp_path / "family_info.tsv"
    lines = [
        f"{family}\tPF00001\t-\t-\t-\tcategory\t{family}_class\n"
        for family in families
    ]
    path.write_text("".join(lines))
    return path


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
    t = Tree("(Mmus_g1:0.1);", format=1)
    d = r2.gene_tree_to_dict(t.get_leaves()[0])
    assert d["leaf"] is True
    assert d["name"] == "Mmus_g1"
    assert d["species"] == "Mmus"

def test_gene_tree_to_dict_internal_og():
    pytest.importorskip("ete3")
    from ete3 import Tree
    t = Tree("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);", format=1)
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
    assert "tree_dict" in rec
    assert "ogs" in rec
    assert sorted(all_species) == ["Hsap", "Mmus", "Nvec"]

def test_load_possvm_trees_empty_dir(tmp_path):
    trees, species = r2.load_possvm_trees(tmp_path)
    assert trees == []
    assert species == []


# ── Heatmap data loaders ─────────────────────────────────────────────────────

def test_load_family_info(tmp_path):
    csv = tmp_path / "info.csv"
    csv.write_text("MyFam\tc1\tc2\tc3\tc4\tc5\tMyClass\n")
    info = r2.load_family_info(str(csv))
    assert info == {"MyFam": "MyClass", "MyClass": "MyClass"}

def test_load_family_info_accepts_genefam_format(tmp_path):
    csv = tmp_path / "genefam.csv"
    csv.write_text("Forkhead\tForkhead\t1.1\t3\tGA\tTF\ttfs\n")
    info = r2.load_family_info(str(csv))
    assert info["Forkhead"] == "tfs"
    assert info["tfs"] == "tfs"

def test_load_family_details_accepts_genefam_format(tmp_path):
    csv = tmp_path / "genefam.csv"
    csv.write_text("Forkhead\tForkhead\t1.1\t3\tGA\tTF\ttfs\n")
    details = r2.load_family_details(str(csv))
    assert details["Forkhead"]["pfam"] == ["Forkhead"]
    assert details["Forkhead"]["category"] == "TF"
    assert details["Forkhead"]["cls"] == "tfs"

def test_build_family_records(tmp_path):
    gl = tmp_path / "Pre.Wnt.genes.list"
    gl.write_text("Mmus_g1\nHsap_g2\nMmus_g3\n")
    recs = r2.build_family_records(tmp_path, {"Wnt": "sig"})
    assert len(recs) == 1
    assert recs[0]["family"] == "Wnt"
    assert recs[0]["class"] == "sig"
    assert recs[0]["total"] == 3
    assert recs[0]["species_counts"]["Mmus"] == 2

def test_build_hg_records(tmp_path):
    fa = tmp_path / "Pre.Wnt.HG001.fasta"
    fa.write_text(">Mmus_g1\nACGT\n>Hsap_g2\nACGT\n")
    recs = r2.build_hg_records(tmp_path, {"Wnt": "sig"})
    assert len(recs) == 1
    assert recs[0]["hg"] == "HG001"
    assert recs[0]["total"] == 2

def test_build_family_records_missing_dir(tmp_path):
    recs = r2.build_family_records(tmp_path / "nope", {})
    assert recs == []


# ── Clade groupings ──────────────────────────────────────────────────────────

def test_parse_clade_groupings():
    pytest.importorskip("ete3")
    # Create a minimal species tree with named internal nodes
    from pathlib import Path
    import tempfile
    nwk = "((A:1,B:1)Clade1:1,(C:1,D:1,E:1)Clade2:1)Root;"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".nwk", delete=False) as f:
        f.write(nwk)
        f.flush()
        groups = r2.parse_clade_groupings(Path(f.name))
    os.unlink(f.name)
    # Root has 2 named children: Clade1 and Clade2, with 5 species total
    assert len(groups) >= 1
    root_entry = [g for g in groups if g["name"] == "Root"]
    assert len(root_entry) == 1
    rg = root_entry[0]["groups"]
    assert rg["A"] == "Clade1"
    assert rg["B"] == "Clade1"
    assert rg["C"] == "Clade2"

def test_parse_clade_groupings_no_file():
    groups = r2.parse_clade_groupings(Path("/nonexistent.nwk"))
    assert groups == []


# ── HTML generation ───────────────────────────────────────────────────────────

def test_html_output_structure(tmp_path):
    """TREE_INDEX and per-HG lazy script tags must be present."""
    pytest.importorskip("ete3")
    nwk = tmp_path / "Pre.TF.HG00001.treefile.ortholog_groups.newick"
    nwk.write_text("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);")
    family_info = write_family_info(tmp_path, "TF")

    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(tmp_path),
        "--family_info", str(family_info),
        "--output", str(out),
    ])
    html = out.read_text()

    assert "Step 2 Report" in html
    assert "TREE_INDEX" in html
    assert 'id="treedata-' in html
    assert "Mmus_g1" in html
    assert "OG_0001" in html

    # TREE_INDEX is valid JSON
    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert len(parsed) == 1
    assert parsed[0]["hg"] == "HG00001"
    assert "tree_dict" not in parsed[0]
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
    family_info = write_family_info(tmp_path, "TF", "RBP")

    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(tmp_path),
        "--family_info", str(family_info),
        "--output", str(out),
    ])
    html = out.read_text()

    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 5000]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])
    families = {r["family"] for r in parsed}
    assert families == {"TF", "RBP"}


def test_lazy_script_tags(tmp_path):
    """Each HG must have exactly one treedata- script tag."""
    pytest.importorskip("ete3")
    for hg, fam in [("HG00001", "TF"), ("HG00002", "RBP")]:
        (tmp_path / f"Pre.{fam}.{hg}.treefile.ortholog_groups.newick").write_text(
            "((Mmus_g1:0.1,Hsap_g1:0.1)OG_X:0.2);"
        )
    family_info = write_family_info(tmp_path, "TF", "RBP")

    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(tmp_path),
        "--family_info", str(family_info),
        "--output", str(out),
    ])
    html = out.read_text()

    assert html.count('id="treedata-') == 2
    for m in re.finditer(
        r'<script type="application/json" id="treedata-[^"]+">([^<]+)</script>', html
    ):
        detail = json.loads(m.group(1))
        assert "tree" in detail
        assert "ogs" in detail


def test_heatmap_data_in_html(tmp_path):
    """FAMILY_DATA and HG_DATA are present when search/cluster dirs provided."""
    pytest.importorskip("ete3")
    # Create search data
    search = tmp_path / "search"
    search.mkdir()
    (search / "Pre.Wnt.genes.list").write_text("Mmus_g1\nHsap_g2\n")
    # Create cluster data
    cluster = tmp_path / "clusters"
    cluster.mkdir()
    (cluster / "Pre.Wnt.HG001.fasta").write_text(">Mmus_g1\nACGT\n>Hsap_g2\nACGT\n")
    # Create possvm data
    possvm = tmp_path / "possvm"
    possvm.mkdir()
    (possvm / "Pre.Wnt.HG001.treefile.ortholog_groups.newick").write_text(
        "((Mmus_g1:0.1,Hsap_g2:0.1)OG_X:0.2);"
    )
    family_info = write_family_info(tmp_path, "Wnt")

    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(possvm),
        "--search_dir", str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_info),
        "--output", str(out),
    ])
    html = out.read_text()

    assert "FAMILY_DATA" in html
    assert "HG_DATA" in html
    # Parse FAMILY_DATA
    start = html.index("const FAMILY_DATA")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    fam = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert len(fam) == 1
    assert fam[0]["family"] == "Wnt"


def test_clade_data_empty_when_no_tree(tmp_path):
    """CLADE_DATA is [] when --species_tree not given."""
    pytest.importorskip("ete3")
    nwk = tmp_path / "Pre.TF.HG00001.treefile.ortholog_groups.newick"
    nwk.write_text("((Mmus_g1:0.1,Hsap_g1:0.1)OG_0001:0.2);")
    family_info = write_family_info(tmp_path, "TF")

    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(tmp_path),
        "--family_info", str(family_info),
        "--output", str(out),
    ])
    html = out.read_text()

    start = html.index("const CLADE_DATA")
    chunk = html[start: start + 200]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert parsed == []


# ── prev-tree (GeneRax toggle) ─────────────────────────────────────────────────

def test_prev_tree_in_lazy_data(tmp_path):
    """When --possvm_prev_dir given, lazy data includes prev_tree/prev_ogs."""
    pytest.importorskip("ete3")
    possvm = tmp_path / "possvm"; possvm.mkdir()
    prev   = tmp_path / "prev";   prev.mkdir()
    # Same HG in both dirs
    for d in [possvm, prev]:
        (d / "Pre.TF.HG001.treefile.ortholog_groups.newick").write_text(
            "((Mmus_g1:0.1,Hsap_g1:0.1)OG_X:0.2);"
        )
    family_info = write_family_info(tmp_path, "TF")
    out = tmp_path / "report.html"
    r2.main(["--possvm_dir", str(possvm), "--possvm_prev_dir", str(prev),
             "--family_info", str(family_info), "--output", str(out)])
    html = out.read_text()

    # Parse the single treedata- script tag
    m = re.search(r'<script type="application/json" id="treedata-[^"]+">([^<]+)</script>', html)
    assert m, "No treedata- script tag found"
    detail = json.loads(m.group(1))
    assert "prev_tree" in detail
    assert "prev_ogs" in detail


def test_has_prev_in_index(tmp_path):
    """TREE_INDEX records have has_prev=True when prev dir has matching HG."""
    pytest.importorskip("ete3")
    possvm = tmp_path / "possvm"; possvm.mkdir()
    prev   = tmp_path / "prev";   prev.mkdir()
    for d in [possvm, prev]:
        (d / "Pre.TF.HG001.treefile.ortholog_groups.newick").write_text(
            "((Mmus_g1:0.1,Hsap_g1:0.1)OG_X:0.2);"
        )
    family_info = write_family_info(tmp_path, "TF")
    out = tmp_path / "report.html"
    r2.main(["--possvm_dir", str(possvm), "--possvm_prev_dir", str(prev),
             "--family_info", str(family_info), "--output", str(out)])
    html = out.read_text()

    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert parsed[0]["has_prev"] is True


def test_has_prev_false_without_prev_dir(tmp_path):
    """TREE_INDEX records have has_prev=False when --possvm_prev_dir not given."""
    pytest.importorskip("ete3")
    (tmp_path / "Pre.TF.HG001.treefile.ortholog_groups.newick").write_text(
        "((Mmus_g1:0.1,Hsap_g1:0.1)OG_X:0.2);"
    )
    family_info = write_family_info(tmp_path, "TF")
    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(tmp_path),
        "--family_info", str(family_info),
        "--output", str(out),
    ])
    html = out.read_text()

    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    parsed = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert parsed[0]["has_prev"] is False


def test_prev_only_report_marks_no_generax(tmp_path):
    """Prev-only HGs should not enable GeneRax-specific UI state or counts."""
    pytest.importorskip("ete3")
    prev = tmp_path / "prev"
    prev.mkdir()
    (prev / "Pre.TF.HG001.treefile.ortholog_groups.newick").write_text(
        "((Mmus_g1:0.1,Hsap_g1:0.1)OG_X:0.2);"
    )
    family_info_path = write_family_info(tmp_path, "TF")

    out = tmp_path / "report.html"
    r2.main([
        "--possvm_dir", str(tmp_path / "missing-possvm"),
        "--possvm_prev_dir", str(prev),
        "--family_info", str(family_info_path),
        "--output", str(out),
    ])
    html = out.read_text()

    start = html.index("const HAVE_GENERAX")
    chunk = html[start: start + 100]
    assert "false" in chunk.lower()

    start = html.index("const TREE_INDEX")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    tree_index = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert tree_index[0]["source"] == "prev"

    start = html.index("const FAMILY_INFO")
    chunk = html[start: start + 3000]
    json_start = chunk.index("[")
    family_info = json.loads(chunk[json_start: chunk.index(";", json_start)])
    assert family_info[0]["n_trees"] == 1
    assert family_info[0]["n_generax"] == 0
