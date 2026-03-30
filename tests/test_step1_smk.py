"""
Tests for step1.smk — a simple Snakemake port of step1.nf.

Three layers:
  1. genefam parsing       — pure Python, no Snakemake required
  2. target/path shaping   — pure Python
  3. Snakemake dry-runs    — skipped when snakemake is not on PATH
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

SMK = Path(__file__).parent.parent / "step1.smk"

HAS_SNAKEMAKE = shutil.which("snakemake") is not None
snakemake_only = pytest.mark.skipif(
    not HAS_SNAKEMAKE, reason="snakemake not on PATH"
)


def parse_genefam(content: str):
    rows = []
    prefs = {}
    for raw in content.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = [c.strip() for c in raw.split("\t")]
        cols = [c for c in cols if c]
        family = cols[0]
        pref = cols[-1]
        rows.append(family)
        prefs[family] = pref
    return rows, prefs


def cluster_targets(families: list, prefs: dict, cluster_dir: str) -> list:
    return [f"{cluster_dir}/{prefs[f]}.{f}_cluster.tsv" for f in families]


def test_parse_genefam_basic():
    families, prefs = parse_genefam("Kinase\tPF00069\tKIN\n")
    assert families == ["Kinase"]
    assert prefs == {"Kinase": "KIN"}


def test_parse_genefam_skips_blank_comment_and_trailing_empty():
    content = "# c\n\nKinase\tPF00069\tKIN\t\nTF\tPF00010\tTRX\n"
    families, prefs = parse_genefam(content)
    assert families == ["Kinase", "TF"]
    assert prefs == {"Kinase": "KIN", "TF": "TRX"}


def test_cluster_targets_shape():
    out = cluster_targets(["Kinase", "TF"], {"Kinase": "KIN", "TF": "TRX"}, "results_step1_smk/clusters")
    assert out == [
        "results_step1_smk/clusters/KIN.Kinase_cluster.tsv",
        "results_step1_smk/clusters/TRX.TF_cluster.tsv",
    ]


def _make_workspace(tmp_path: Path) -> Path:
    genefam = tmp_path / "genefam.tsv"
    genefam.write_text("Kinase\tPF00069\tKIN\nTF\tPF00010\tTRX\n")
    infasta = tmp_path / "input.fasta"
    infasta.write_text(">KIN_gene1\nACGT\n>TRX_gene1\nACGT\n")
    species_list = tmp_path / "species_list"
    species_list.write_text("KIN\nTRX\n")
    cfg = {
        "genefam_info": str(genefam),
        "species_list": str(species_list),
        "infasta": str(infasta),
        "use_prepare": "0",
        "prepare_db_dir": "/tmp/db",
        "search_dir": str(tmp_path / "step1_out" / "search"),
        "cluster_dir": str(tmp_path / "step1_out" / "clusters"),
        "pfam_db": "/tmp/Pfam-A.hmm",
        "domain_expand": "30",
        "s1_ncpu": "2",
        "s2_ncpu": "2",
        "s2_inflation": "1.1",
        "max_n": "3000",
    }
    cfg_file = tmp_path / "test_cfg.yaml"
    cfg_file.write_text("\n".join(f"{k}: {v}" for k, v in cfg.items()) + "\n")
    return cfg_file


def _dry_run(tmp_path: Path) -> subprocess.CompletedProcess:
    cfg_file = _make_workspace(tmp_path)
    env = dict(os.environ)
    env["XDG_CACHE_HOME"] = "/tmp/.cache"
    return subprocess.run(
        [
            "snakemake", "-s", str(SMK),
            "--dry-run", "--forceall",
            "--configfile", str(cfg_file),
        ],
        capture_output=True, text=True, cwd=str(tmp_path), env=env,
    )


@snakemake_only
def test_dryrun_step1_succeeds(tmp_path):
    r = _dry_run(tmp_path)
    assert r.returncode == 0, r.stderr


@snakemake_only
def test_dryrun_step1_mentions_search_and_cluster(tmp_path):
    r = _dry_run(tmp_path)
    assert r.returncode == 0, r.stderr
    stdout = r.stdout.lower()
    assert "rule search" in stdout
    assert "rule cluster" in stdout


@snakemake_only
def test_dryrun_step1_targets_cluster_tsvs(tmp_path):
    r = _dry_run(tmp_path)
    assert r.returncode == 0, r.stderr
    assert "_cluster.tsv" in r.stdout
