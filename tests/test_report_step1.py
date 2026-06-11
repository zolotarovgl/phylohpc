"""Tests for workflow/report_step1.py"""

import json
import os
import sys
import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import report_step1 as r1


# ── get_species_prefix ────────────────────────────────────────────────────────

def test_prefix_underscore():
    assert r1.get_species_prefix("Mmus_ENSMUSP001") == "Mmus"


def test_prefix_dot():
    assert r1.get_species_prefix("Nvec.NVE12345") == "Nvec"


def test_prefix_no_sep():
    assert r1.get_species_prefix("Mmus") == "Mmus"


def test_prefix_multiple_underscores():
    assert r1.get_species_prefix("Hsap_ENSP00000_extra") == "Hsap"


# ── parse_genes_list ──────────────────────────────────────────────────────────

def test_parse_genes_list_basic(tmp_path):
    f = tmp_path / "tfs.bZIP.genes.list"
    f.write_text("Mmus_g1\nHsap_g2\nNvec_g3\n")
    assert r1.parse_genes_list(f) == ["Mmus_g1", "Hsap_g2", "Nvec_g3"]


def test_parse_genes_list_skips_comments(tmp_path):
    f = tmp_path / "x.genes.list"
    f.write_text("# comment\nMmus_g1\n")
    assert r1.parse_genes_list(f) == ["Mmus_g1"]


def test_parse_genes_list_missing():
    assert r1.parse_genes_list("/nonexistent/path.txt") == []


def test_parse_genes_list_empty(tmp_path):
    f = tmp_path / "empty.list"
    f.write_text("")
    assert r1.parse_genes_list(f) == []


# ── parse_fasta_species ───────────────────────────────────────────────────────

def test_parse_fasta_species_basic(tmp_path):
    f = tmp_path / "hg.fasta"
    f.write_text(">Mmus_g1 desc\nACGT\n>Hsap_g2\nACGT\n>Mmus_g3\nACGT\n")
    counts = r1.parse_fasta_species(f)
    assert counts["Mmus"] == 2
    assert counts["Hsap"] == 1


def test_parse_fasta_species_missing():
    assert r1.parse_fasta_species("/nonexistent/file.fasta") == {}


def test_parse_fasta_species_empty(tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    assert r1.parse_fasta_species(f) == {}


# ── load_family_info ──────────────────────────────────────────────────────────

def test_load_family_info_basic(tmp_path):
    csv = tmp_path / "info.csv"
    csv.write_text("bZIP\tprofile\t1.1\t3\tGA\tTF\ttfs\n"
                   "Chromo\tprofile\t1.1\t3\tGA\tchromatin\tchr\n")
    info = r1.load_family_info(csv)
    assert info["bZIP"] == "tfs"
    assert info["Chromo"] == "chr"


def test_load_family_info_missing():
    assert r1.load_family_info("/nonexistent.csv") == {}


def test_load_family_info_short_lines(tmp_path):
    csv = tmp_path / "info.csv"
    csv.write_text("bZIP\tprofile\n")  # only 2 columns — should be skipped
    assert r1.load_family_info(csv) == {}


# ── build_family_records ──────────────────────────────────────────────────────

def test_build_family_records_basic(tmp_path):
    search = tmp_path / "search"
    search.mkdir()
    (search / "tfs.bZIP.genes.list").write_text(
        "Mmus_g1\nMmus_g2\nHsap_g1\n"
    )
    family_info = {"bZIP": "tfs"}
    records = r1.build_family_records(search, family_info)
    assert len(records) == 1
    rec = records[0]
    assert rec["id"] == "tfs.bZIP"
    assert rec["class"] == "tfs"
    assert rec["species_counts"]["Mmus"] == 2
    assert rec["species_counts"]["Hsap"] == 1
    assert rec["total"] == 3


def test_build_family_records_empty_file_skipped(tmp_path):
    search = tmp_path / "search"
    search.mkdir()
    (search / "tfs.bZIP.genes.list").write_text("")
    records = r1.build_family_records(search, {})
    assert records == []


def test_build_family_records_no_dir(tmp_path):
    records = r1.build_family_records(tmp_path / "nonexistent", {})
    assert records == []


def test_build_family_records_multiple(tmp_path):
    search = tmp_path / "search"
    search.mkdir()
    (search / "tfs.bZIP.genes.list").write_text("Mmus_g1\n")
    (search / "chr.Chromo.genes.list").write_text("Hsap_g1\nNvec_g1\n")
    records = r1.build_family_records(search, {"bZIP": "tfs", "Chromo": "chr"})
    ids = {r["id"] for r in records}
    assert ids == {"tfs.bZIP", "chr.Chromo"}


# ── build_hg_records ──────────────────────────────────────────────────────────

def test_build_hg_records_basic(tmp_path):
    cluster = tmp_path / "clusters"
    cluster.mkdir()
    (cluster / "tfs.bZIP.1.fasta").write_text(
        ">Mmus_g1\nACGT\n>Hsap_g1\nACGT\n"
    )
    records = r1.build_hg_records(cluster, {"bZIP": "tfs"})
    assert len(records) == 1
    rec = records[0]
    assert rec["id"] == "tfs.bZIP.1"
    assert rec["class"] == "tfs"
    assert rec["species_counts"]["Mmus"] == 1
    assert rec["species_counts"]["Hsap"] == 1


def test_build_hg_records_empty_fasta_skipped(tmp_path):
    cluster = tmp_path / "clusters"
    cluster.mkdir()
    (cluster / "tfs.bZIP.1.fasta").write_text("")
    records = r1.build_hg_records(cluster, {})
    assert records == []


def test_build_hg_records_no_dir(tmp_path):
    assert r1.build_hg_records(tmp_path / "missing", {}) == []


# ── main() / HTML generation ──────────────────────────────────────────────────

def write_test_inputs(tmp_path):
    search = tmp_path / "search"
    cluster = tmp_path / "clusters"
    search.mkdir()
    cluster.mkdir()

    (search / "tfs.bZIP.genes.list").write_text(
        "Mmus_g1\nMmus_g2\nHsap_g1\nDmel_g1\n"
    )
    (search / "chr.Chromo.genes.list").write_text(
        "Mmus_g3\nHsap_g2\n"
    )
    (cluster / "tfs.bZIP.HG001.fasta").write_text(
        ">Mmus_g1\nACGT\n>Hsap_g1\nACGT\n"
    )

    family_csv = tmp_path / "info.csv"
    family_csv.write_text(
        "bZIP\tprofile\t1.1\t3\tGA\tTF\ttfs\n"
        "Chromo\tprofile\t1.1\t3\tGA\tchromatin\tchr\n"
    )
    return search, cluster, family_csv


def test_html_generated(tmp_path):
    search, cluster, family_csv = write_test_inputs(tmp_path)
    out = tmp_path / "report.html"
    with patch("sys.argv", [
        "report_step1",
        "--search_dir",  str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_csv),
        "--tree",        str(tmp_path / "nonexistent.nwk"),   # no tree → empty order
        "--output",      str(out),
    ]):
        r1.main()
    assert out.exists()
    assert out.stat().st_size > 2000


def test_html_no_placeholders(tmp_path):
    """All %%PLACEHOLDER%% tokens must have been substituted."""
    search, cluster, family_csv = write_test_inputs(tmp_path)
    out = tmp_path / "report.html"
    with patch("sys.argv", [
        "report_step1",
        "--search_dir",  str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_csv),
        "--tree",        str(tmp_path / "nonexistent.nwk"),
        "--output",      str(out),
    ]):
        r1.main()
    html = out.read_text()
    assert "%%" not in html


def test_html_contains_family_name(tmp_path):
    search, cluster, family_csv = write_test_inputs(tmp_path)
    out = tmp_path / "report.html"
    with patch("sys.argv", [
        "report_step1",
        "--search_dir",  str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_csv),
        "--tree",        str(tmp_path / "nonexistent.nwk"),
        "--output",      str(out),
    ]):
        r1.main()
    html = out.read_text()
    assert "bZIP" in html
    assert "Chromo" in html


def test_html_embeds_valid_json(tmp_path):
    """FAMILY_DATA and HG_DATA must be valid JSON arrays."""
    search, cluster, family_csv = write_test_inputs(tmp_path)
    out = tmp_path / "report.html"
    with patch("sys.argv", [
        "report_step1",
        "--search_dir",  str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_csv),
        "--tree",        str(tmp_path / "nonexistent.nwk"),
        "--output",      str(out),
    ]):
        r1.main()
    html = out.read_text()
    # Extract FAMILY_DATA JSON and parse it
    import re
    m = re.search(r"const FAMILY_DATA\s*=\s*(\[.*?\]);", html, re.DOTALL)
    assert m, "FAMILY_DATA not found in HTML"
    data = json.loads(m.group(1))
    assert isinstance(data, list)
    assert any(d["id"] == "tfs.bZIP" for d in data)


def test_html_empty_inputs(tmp_path):
    """Empty search and cluster dirs should still produce valid HTML."""
    search = tmp_path / "search"; search.mkdir()
    cluster = tmp_path / "clusters"; cluster.mkdir()
    family_csv = tmp_path / "info.csv"; family_csv.write_text("")
    out = tmp_path / "report.html"
    with patch("sys.argv", [
        "report_step1",
        "--search_dir",  str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_csv),
        "--tree",        str(tmp_path / "nonexistent.nwk"),
        "--output",      str(out),
    ]):
        r1.main()
    assert out.exists()
    html = out.read_text()
    assert "%%" not in html


def test_html_with_tree(tmp_path):
    """When a valid newick tree is supplied, SPECIES_ORDER should be non-empty."""
    pytest.importorskip("ete3")
    search, cluster, family_csv = write_test_inputs(tmp_path)
    tree_f = tmp_path / "tree.nwk"
    tree_f.write_text("((Mmus:1,Hsap:1)Mammalia:1,Dmel:1)Root;")
    out = tmp_path / "report.html"
    with patch("sys.argv", [
        "report_step1",
        "--search_dir",  str(search),
        "--cluster_dir", str(cluster),
        "--family_info", str(family_csv),
        "--tree",        str(tree_f),
        "--output",      str(out),
    ]):
        r1.main()
    html = out.read_text()
    import re
    m = re.search(r"const SPECIES_ORDER\s*=\s*(\[.*?\]);", html, re.DOTALL)
    assert m
    order = json.loads(m.group(1))
    assert "Mmus" in order
    assert "Hsap" in order
