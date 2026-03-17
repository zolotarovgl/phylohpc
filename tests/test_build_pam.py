"""Tests for workflow/build_pam.py"""

import os
import sys
import pytest
import pandas as pd
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import build_pam as bp


# ── get_species_prefix ────────────────────────────────────────────────────────

def test_prefix_underscore():
    assert bp.get_species_prefix("Mmus_ENSMUSP001") == "Mmus"

def test_prefix_dot():
    assert bp.get_species_prefix("Nvec.NVE12345") == "Nvec"

def test_prefix_no_separator():
    assert bp.get_species_prefix("Mmus") == "Mmus"

def test_prefix_multiple_underscores():
    assert bp.get_species_prefix("Hsap_ENSP00000_extra") == "Hsap"


# ── load_gene_og_table ────────────────────────────────────────────────────────

def test_load_basic(tmp_path):
    csv = tmp_path / "hg1.ortholog_groups.csv"
    csv.write_text(
        "Mmus_gene1\tOG1\tref_OG1\tGene1\n"
        "Hsap_gene2\tOG1\tref_OG1\tGene1\n"
        "Nvec_gene3\tOG2\tref_OG2\tGene2\n"
    )
    df = bp.load_gene_og_table([str(csv)], {"Mmus", "Hsap", "Nvec"})
    assert len(df) == 3
    assert set(df["og"]) == {"OG1", "OG2"}


def test_load_filters_out_of_clade(tmp_path):
    csv = tmp_path / "hg1.ortholog_groups.csv"
    csv.write_text(
        "Mmus_gene1\tOG1\tref\tGene1\n"
        "Outgroup_gene\tOG1\tref\tGene1\n"
    )
    df = bp.load_gene_og_table([str(csv)], {"Mmus"})
    assert len(df) == 1
    assert df.iloc[0]["species"] == "Mmus"


def test_load_empty_file(tmp_path):
    csv = tmp_path / "empty.csv"
    csv.write_text("")
    df = bp.load_gene_og_table([str(csv)], {"Mmus"})
    assert df.empty


def test_load_missing_file(tmp_path):
    df = bp.load_gene_og_table([str(tmp_path / "nonexistent.csv")], {"Mmus"})
    assert df.empty


def test_load_skips_comment_lines(tmp_path):
    csv = tmp_path / "hg.csv"
    csv.write_text(
        "# header comment\n"
        "Mmus_g1\tOG1\tref\tGene1\n"
    )
    df = bp.load_gene_og_table([str(csv)], {"Mmus"})
    assert len(df) == 1


def test_load_multiple_files(tmp_path):
    csv1 = tmp_path / "hg1.csv"
    csv2 = tmp_path / "hg2.csv"
    csv1.write_text("Mmus_g1\tOG1\tref\tGene1\n")
    csv2.write_text("Mmus_g2\tOG2\tref\tGene2\n")
    df = bp.load_gene_og_table([str(csv1), str(csv2)], {"Mmus"})
    assert len(df) == 2
    assert set(df["og"]) == {"OG1", "OG2"}


# ── PAM construction (via main()) ─────────────────────────────────────────────

def run_build_pam(tmp_path, csv_contents, species, min_presence=1):
    """Helper: write CSV files and species list, run main(), return PAM df."""
    csv_paths = []
    for i, content in enumerate(csv_contents):
        p = tmp_path / f"hg{i}.csv"
        p.write_text(content)
        csv_paths.append(str(p))

    sp_file = tmp_path / "species.txt"
    sp_file.write_text("\n".join(species) + "\n")

    out_file = tmp_path / "pam.tsv"

    import sys
    from unittest.mock import patch
    args = ["build_pam",
            "--csvs", *csv_paths,
            "--species", str(sp_file),
            "--min_presence", str(min_presence),
            "--output", str(out_file)]
    with patch("sys.argv", args):
        bp.main()

    return pd.read_csv(out_file, sep="\t", index_col=0)


def test_pam_binary_values(tmp_path):
    content = (
        "Mmus_g1\tOG1\tref\tGene1\n"
        "Hsap_g1\tOG1\tref\tGene1\n"
        "Nvec_g1\tOG2\tref\tGene2\n"
    )
    pam = run_build_pam(tmp_path, [content], ["Mmus", "Hsap", "Nvec"])
    assert set(pam.values.flatten()).issubset({0, 1})


def test_pam_correct_presence(tmp_path):
    content = (
        "Mmus_g1\tOG1\tref\tGene1\n"
        "Hsap_g1\tOG1\tref\tGene1\n"
    )
    pam = run_build_pam(tmp_path, [content], ["Mmus", "Hsap", "Nvec"])
    assert pam.loc["Mmus", "OG1"] == 1
    assert pam.loc["Hsap", "OG1"] == 1
    assert pam.loc["Nvec", "OG1"] == 0


def test_pam_all_species_present_in_rows(tmp_path):
    content = "Mmus_g1\tOG1\tref\tGene1\n"
    pam = run_build_pam(tmp_path, [content], ["Mmus", "Hsap", "Nvec"])
    assert set(pam.index) == {"Mmus", "Hsap", "Nvec"}


def test_pam_min_presence_filter(tmp_path):
    # OG1 present in 2 species, OG2 only in 1
    content = (
        "Mmus_g1\tOG1\tref\tGene1\n"
        "Hsap_g1\tOG1\tref\tGene1\n"
        "Nvec_g1\tOG2\tref\tGene2\n"
    )
    pam = run_build_pam(tmp_path, [content], ["Mmus", "Hsap", "Nvec"], min_presence=2)
    assert "OG1" in pam.columns
    assert "OG2" not in pam.columns


def test_pam_multiple_genes_same_og_counts_once(tmp_path):
    """Two genes from the same species in the same OG → presence = 1, not 2."""
    content = (
        "Mmus_g1\tOG1\tref\tGene1\n"
        "Mmus_g2\tOG1\tref\tGene1\n"
    )
    pam = run_build_pam(tmp_path, [content], ["Mmus", "Hsap"])
    assert pam.loc["Mmus", "OG1"] == 1


def test_pam_empty_csvs_writes_empty_matrix(tmp_path):
    pam = run_build_pam(tmp_path, [""], ["Mmus", "Hsap"])
    assert pam.shape[1] == 0          # no OG columns
    assert set(pam.index) == {"Mmus", "Hsap"}
