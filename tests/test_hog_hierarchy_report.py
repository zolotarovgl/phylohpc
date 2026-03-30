"""Tests for the hOG hierarchy report builder."""

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))

from build_hog_report import process_hg
from visualize_hog_hierarchy import build_data, og_short_py


def write_possvm_csv(path: Path, rows: list[tuple[str, str, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "gene\torthogroup\torthogroup_support\treference_ortholog\treference_support\n"
        + "\n".join("\t".join(row) for row in rows)
        + "\n"
    )


def write_in_species(path: Path, species: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(species) + "\n")


def test_og_short_py_preserves_like_label():
    og_id = "Bilateria.tfs.Forkhead.HG6.0:like:Foxj1"
    assert og_short_py(og_id, "tfs.Forkhead.HG6", "Bilateria") == "like:Foxj1"


def test_build_hog_report_preserves_like_orthogroups(tmp_path):
    ancestry_dir = tmp_path / "ancestry"
    levels = ["Metazoa", "Bilateria"]
    hg = "tfs.Forkhead.HG6"

    write_in_species(ancestry_dir / "Metazoa" / "Metazoa.in_species.txt", ["A", "B", "C"])
    write_in_species(ancestry_dir / "Bilateria" / "Bilateria.in_species.txt", ["A", "B"])

    write_possvm_csv(
        ancestry_dir / "Metazoa" / "possvm" / f"{hg}.ortholog_groups.csv",
        [
            ("A_g1", "Metazoa.tfs.Forkhead.HG6.0:like:Foxj1", "1.0", "NA", "NA"),
            ("B_g1", "Metazoa.tfs.Forkhead.HG6.0:like:Foxj1", "1.0", "NA", "NA"),
            ("C_g1", "Metazoa.tfs.Forkhead.HG6.1:Foxj1", "1.0", "Foxj1", "1.0"),
        ],
    )
    write_possvm_csv(
        ancestry_dir / "Bilateria" / "possvm" / f"{hg}.ortholog_groups.csv",
        [
            ("A_g1", "Bilateria.tfs.Forkhead.HG6.0:like:Foxj1", "1.0", "NA", "NA"),
            ("B_g1", "Bilateria.tfs.Forkhead.HG6.0:like:Foxj1", "1.0", "NA", "NA"),
            ("A_g2", "Bilateria.tfs.Forkhead.HG6.1:Foxj1", "1.0", "Foxj1", "1.0"),
        ],
    )

    stats, links = process_hg(hg, levels, ancestry_dir)
    stats_df = pd.DataFrame(stats)
    links_df = pd.DataFrame(links)

    like_ogs = {row["og"] for row in stats if ":like:" in row["og"]}
    assert like_ogs == {
        "Metazoa.tfs.Forkhead.HG6.0:like:Foxj1",
        "Bilateria.tfs.Forkhead.HG6.0:like:Foxj1",
    }
    totals = {(row["level"], row["og"]): row["n_total_species"] for row in stats}
    assert totals[("Metazoa", "Metazoa.tfs.Forkhead.HG6.0:like:Foxj1")] == 3
    assert totals[("Bilateria", "Bilateria.tfs.Forkhead.HG6.0:like:Foxj1")] == 2

    data = build_data(stats_df, links_df, levels, trees={})
    meta = data["hgs"][hg]

    node_map = {node["id"]: node for node in meta["nodes"]}
    assert "Metazoa.tfs.Forkhead.HG6.0:like:Foxj1" in node_map
    assert "Bilateria.tfs.Forkhead.HG6.0:like:Foxj1" in node_map
    assert node_map["Metazoa.tfs.Forkhead.HG6.0:like:Foxj1"]["n_total_species"] == 3
    assert node_map["Bilateria.tfs.Forkhead.HG6.0:like:Foxj1"]["n_total_species"] == 2
    assert "like:foxj1" in meta["search_text"].lower()
