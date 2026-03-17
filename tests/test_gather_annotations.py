"""Tests for workflow/gather_annotations.py"""

import os
import sys
import pytest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import gather_annotations as ga


# ---------------------------------------------------------------------------
# load_species_ids
# ---------------------------------------------------------------------------

def test_load_species_ids_from_file(tmp_path):
    f = tmp_path / "ids.txt"
    f.write_text("Hsap\nMmus\nDrer\n")
    result = ga.load_species_ids(str(f))
    assert result == ["Hsap", "Mmus", "Drer"]


def test_load_species_ids_from_file_strips_blank_lines(tmp_path):
    f = tmp_path / "ids.txt"
    f.write_text("Hsap\n\nMmus\n")
    result = ga.load_species_ids(str(f))
    assert result == ["Hsap", "Mmus"]


def test_load_species_ids_from_string():
    result = ga.load_species_ids("Hsap")
    assert result == ["Hsap"]


# ---------------------------------------------------------------------------
# collect_tmp_anno
# ---------------------------------------------------------------------------

def test_collect_tmp_anno_basic(tmp_path):
    tree_dir = tmp_path / "tree_out"
    tree_dir.mkdir()
    csv = tree_dir / "pfam.groups.csv"
    csv.write_text(
        "Hsap_gene1\tOG1\thomolog\t1.0\n"
        "Mmus_gene2\tOG2\tortho\t0.9\n"
    )
    result = ga.collect_tmp_anno([str(tree_dir)], None, "Hsap")
    assert "Hsap_gene1" in result
    assert result["Hsap_gene1"] == ["OG1", "homolog", "1.0"]
    assert "Mmus_gene2" not in result


def test_collect_tmp_anno_first_dir_wins(tmp_path):
    """First tree_dir in the list has priority (first-wins semantics)."""
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()
    dir2.mkdir()

    (dir1 / "pfam.groups.csv").write_text("Hsap_g1\tOG1\thomolog\t1.0\n")
    (dir2 / "pfam.groups.csv").write_text("Hsap_g1\tOG2\tortho\t0.8\n")

    result = ga.collect_tmp_anno([str(dir1), str(dir2)], None, "Hsap")
    assert result["Hsap_g1"][0] == "OG1"


def test_collect_tmp_anno_prefix_filter(tmp_path):
    tree_dir = tmp_path / "tree_out"
    tree_dir.mkdir()
    (tree_dir / "fam1.groups.csv").write_text("Hsap_g1\tOG1\th\t1.0\n")
    (tree_dir / "fam2.groups.csv").write_text("Hsap_g2\tOG2\th\t1.0\n")

    result = ga.collect_tmp_anno([str(tree_dir)], "fam1", "Hsap")
    assert "Hsap_g1" in result
    assert "Hsap_g2" not in result


def test_collect_tmp_anno_missing_dir_ignored(tmp_path):
    result = ga.collect_tmp_anno([str(tmp_path / "nonexistent")], None, "Hsap")
    assert result == {}


# ---------------------------------------------------------------------------
# collect_ids_todo
# ---------------------------------------------------------------------------

def test_collect_ids_todo_basic(tmp_path):
    search_dir = tmp_path / "search"
    search_dir.mkdir()
    (search_dir / "pfam.genes.list").write_text(
        "Hsap_gene1 extra\nMmus_gene2 extra\nHsap_gene3 extra\n"
    )
    result = ga.collect_ids_todo(str(search_dir), None, "Hsap")
    assert sorted(result) == ["Hsap_gene1", "Hsap_gene3"]


def test_collect_ids_todo_prefix_filter(tmp_path):
    search_dir = tmp_path / "search"
    search_dir.mkdir()
    (search_dir / "fam1.genes.list").write_text("Hsap_g1 x\n")
    (search_dir / "fam2.genes.list").write_text("Hsap_g2 x\n")
    result = ga.collect_ids_todo(str(search_dir), "fam1", "Hsap")
    assert result == ["Hsap_g1"]


def test_collect_ids_todo_empty_dir(tmp_path):
    search_dir = tmp_path / "empty_search"
    search_dir.mkdir()
    result = ga.collect_ids_todo(str(search_dir), None, "Hsap")
    assert result == []


# ---------------------------------------------------------------------------
# collect_pep2hg
# ---------------------------------------------------------------------------

def test_collect_pep2hg_basic(tmp_path):
    cluster_dir = tmp_path / "clusters"
    cluster_dir.mkdir()
    (cluster_dir / "pfam_cluster.tsv").write_text(
        "HG001\tHsap_gene1\n"
        "HG001\tMmus_gene2\n"
        "HG002\tHsap_gene3\n"
    )
    result = ga.collect_pep2hg(str(cluster_dir), "Hsap")
    assert result == {
        "Hsap_gene1": "pfam.HG001",
        "Hsap_gene3": "pfam.HG002",
    }
    assert "Mmus_gene2" not in result


def test_collect_pep2hg_missing_dir_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ga.collect_pep2hg(str(tmp_path / "nonexistent"), "Hsap")


def test_collect_pep2hg_short_lines_skipped(tmp_path):
    cluster_dir = tmp_path / "clusters"
    cluster_dir.mkdir()
    (cluster_dir / "pfam_cluster.tsv").write_text("only_one_field\nHG1\tHsap_g1\n")
    result = ga.collect_pep2hg(str(cluster_dir), "Hsap")
    assert "Hsap_g1" in result


# ---------------------------------------------------------------------------
# build_result
# ---------------------------------------------------------------------------

def test_build_result_annotated_gene():
    tmp_anno = {"Hsap_g1": ["OG1", "homolog", "1.0"]}
    ids_todo = ["Hsap_g1"]
    gene2class = {}
    pep2hg = {}
    rows = ga.build_result(tmp_anno, ids_todo, gene2class, pep2hg)
    assert rows == [["Hsap_g1", "OG1", "homolog", "1.0"]]


def test_build_result_unclassified_with_hg():
    tmp_anno = {}
    ids_todo = ["Hsap_g2"]
    gene2class = {}
    pep2hg = {"Hsap_g2": "pfam.HG001"}
    rows = ga.build_result(tmp_anno, ids_todo, gene2class, pep2hg)
    assert rows[0][1] == "pfam.HG001:Unclassified"


def test_build_result_unclassified_with_gene_class():
    tmp_anno = {}
    ids_todo = ["Hsap_g3"]
    gene2class = {"Hsap_g3": "fam1"}
    pep2hg = {}
    rows = ga.build_result(tmp_anno, ids_todo, gene2class, pep2hg)
    assert rows[0][1] == "fam1:Unclassified"


def test_build_result_fully_unclassified():
    rows = ga.build_result({}, ["Hsap_g4"], {}, {})
    assert rows == [["Hsap_g4", "Unclassified"]]


def test_build_result_hg_takes_priority_over_gene_class():
    """When a gene is unclassified but has both hg and gene_class, hg wins."""
    tmp_anno = {}
    ids_todo = ["Hsap_g5"]
    gene2class = {"Hsap_g5": "fam1"}
    pep2hg = {"Hsap_g5": "pfam.HG002"}
    rows = ga.build_result(tmp_anno, ids_todo, gene2class, pep2hg)
    assert rows[0][1] == "pfam.HG002:Unclassified"


# ---------------------------------------------------------------------------
# write_split_outputs
# ---------------------------------------------------------------------------

def test_write_split_outputs_creates_files(tmp_path):
    rows = [
        ["Hsap_g1", "pfam.HG1:Classified", "x"],
        ["Hsap_g2", "pfam.HG2:Classified", "y"],
        ["Hsap_g3", "Unclassified"],
    ]
    ga.write_split_outputs(rows, "Hsap", tmp_path)

    assert (tmp_path / "Hsap.pfam.tsv").exists()
    assert (tmp_path / "Hsap.Unclassified.tsv").exists()


def test_write_split_outputs_correct_content(tmp_path):
    rows = [
        ["Hsap_g1", "pfam.HG1:OG1", "homolog"],
        ["Hsap_g2", "pfam.HG2:OG2", "ortho"],
    ]
    ga.write_split_outputs(rows, "Hsap", tmp_path)
    content = (tmp_path / "Hsap.pfam.tsv").read_text()
    assert "Hsap_g1" in content
    assert "Hsap_g2" in content
