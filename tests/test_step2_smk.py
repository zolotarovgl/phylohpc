"""
Tests for step2.smk — the Snakemake port of step2.nf.

Three layers:
  1. Config/flag coercion  — pure Python, no Snakemake required
  2. _report_inputs logic  — pure Python reimplementation
  3. IDs file parsing      — pure Python
  4. Snakemake dry-runs    — skipped when snakemake is not on PATH
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

# Path to the workflow file under test
SMK = Path(__file__).parent.parent / "step2.smk"

HAS_SNAKEMAKE = shutil.which("snakemake") is not None
snakemake_only = pytest.mark.skipif(
    not HAS_SNAKEMAKE, reason="snakemake not on PATH"
)


# ── helpers that mirror step2.smk logic ───────────────────────────────────────

def coerce_run_generax(value) -> bool:
    """Mirror: bool(int(config.get('run_generax', 0)))"""
    return bool(int(value))


def pvm_infix(run_generax: bool) -> str:
    """Mirror: PVM_INFIX = 'generax' if RUN_GENERAX else 'treefile'"""
    return "generax" if run_generax else "treefile"


def parse_ids(content: str) -> list:
    """Mirror the IDs-file parsing in step2.smk."""
    return [line.strip() for line in content.splitlines()
            if line.strip() and not line.startswith("#")]


def report_inputs(ids: list, outdir: str, run_generax: bool) -> list:
    """
    Pure reimplementation of _report_inputs(wildcards) from step2.smk.
    Returns the list of expected input paths for the report rule.
    """
    infix = pvm_infix(run_generax)
    pvm = [f"{outdir}/possvm/{i}.{infix}.ortholog_groups.newick" for i in ids]
    pvm_prev = (
        [f"{outdir}/possvm_prev/{i}.treefile.ortholog_groups.newick" for i in ids]
        if run_generax else []
    )
    return pvm + pvm_prev


# ── 1. config / flag coercion ─────────────────────────────────────────────────

def test_run_generax_zero_is_false():
    assert coerce_run_generax("0") is False


def test_run_generax_one_is_true():
    assert coerce_run_generax("1") is True


def test_run_generax_int_zero_is_false():
    assert coerce_run_generax(0) is False


def test_run_generax_int_one_is_true():
    assert coerce_run_generax(1) is True


def test_pvm_infix_without_generax():
    assert pvm_infix(False) == "treefile"


def test_pvm_infix_with_generax():
    assert pvm_infix(True) == "generax"


# ── 2. IDs file parsing ───────────────────────────────────────────────────────

def test_ids_basic():
    assert parse_ids("HG001\nHG002\n") == ["HG001", "HG002"]


def test_ids_skips_blank_lines():
    assert parse_ids("HG001\n\nHG002\n") == ["HG001", "HG002"]


def test_ids_skips_comment_lines():
    assert parse_ids("# comment\nHG001\n# another\nHG002\n") == ["HG001", "HG002"]


def test_ids_strips_whitespace():
    assert parse_ids("  HG001  \n\tHG002\t\n") == ["HG001", "HG002"]


def test_ids_empty_file():
    assert parse_ids("") == []


def test_ids_only_comments_and_blanks():
    assert parse_ids("# comment\n\n# another\n") == []


# ── 3. _report_inputs logic ───────────────────────────────────────────────────

def test_report_inputs_no_generax_two_ids():
    out = report_inputs(["HG1", "HG2"], "results", False)
    assert len(out) == 2
    assert all("possvm/" in p and "treefile" in p for p in out)
    assert not any("possvm_prev" in p for p in out)


def test_report_inputs_with_generax_two_ids():
    out = report_inputs(["HG1", "HG2"], "results", True)
    assert len(out) == 4
    pvm_paths = [p for p in out if "possvm_prev" not in p]
    prev_paths = [p for p in out if "possvm_prev" in p]
    assert len(pvm_paths) == 2
    assert len(prev_paths) == 2
    assert all("generax" in p for p in pvm_paths)
    assert all("treefile" in p for p in prev_paths)


def test_report_inputs_generax_prev_paths_use_treefile_infix():
    """possvm_prev always uses the treefile infix, even when GeneRax is on."""
    out = report_inputs(["HG1"], "results", True)
    prev = [p for p in out if "possvm_prev" in p]
    assert len(prev) == 1
    assert "treefile" in prev[0]


def test_report_inputs_empty_ids_no_generax():
    assert report_inputs([], "results", False) == []


def test_report_inputs_empty_ids_with_generax():
    assert report_inputs([], "results", True) == []


def test_report_inputs_correct_outdir_prefix():
    out = report_inputs(["HG1"], "/data/out", False)
    assert out[0].startswith("/data/out/possvm/")


# ── 4. Snakemake dry-run helpers ──────────────────────────────────────────────

def _make_workspace(tmp_path: Path, ids: list, run_generax: bool = False) -> Path:
    """Create a minimal directory + config for a Snakemake dry-run."""
    (tmp_path / "results" / "clusters").mkdir(parents=True)
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text("\n".join(ids) + "\n")
    for hg in ids:
        (tmp_path / "results" / "clusters" / f"{hg}.fasta").write_text(
            f">{hg}_seq\nACGT\n"
        )
    cfg = {
        "ids":           str(ids_file),
        "outdir":        str(tmp_path / "results"),
        "refnames":      str(tmp_path / "refnames.csv"),
        "species_tree":  str(tmp_path / "sptree.nwk"),
        "refspecies":    "Mmus",
        "tree_method":   "fasttree",
        "iqtree2_model": "LG",
        "mafft_opt":     '""',
        "run_generax":   "1" if run_generax else "0",
        "ncpu_generax":  "4",
        "max_spr":       "2",
        "subs_model":    "LG",
    }
    cfg_file = tmp_path / "test_cfg.yaml"
    cfg_file.write_text(
        "\n".join(f"{k}: {v}" for k, v in cfg.items()) + "\n"
    )
    return cfg_file


def _dry_run(tmp_path: Path, ids: list, run_generax: bool = False) -> subprocess.CompletedProcess:
    cfg_file = _make_workspace(tmp_path, ids, run_generax)
    return subprocess.run(
        [
            "snakemake", "-s", str(SMK),
            "--dry-run", "--forceall",
            "--configfile", str(cfg_file),
        ],
        capture_output=True, text=True, cwd=str(tmp_path),
    )


# ── 5. Snakemake dry-run tests ────────────────────────────────────────────────

@snakemake_only
def test_dryrun_no_generax_succeeds(tmp_path):
    r = _dry_run(tmp_path, ["HG001"], run_generax=False)
    assert r.returncode == 0, r.stderr


@snakemake_only
def test_dryrun_no_generax_omits_generax_rule(tmp_path):
    r = _dry_run(tmp_path, ["HG001"], run_generax=False)
    assert r.returncode == 0, r.stderr
    # The word "generax" should not appear as a job name
    assert "rule generax" not in r.stdout.lower()


@snakemake_only
def test_dryrun_with_generax_succeeds(tmp_path):
    r = _dry_run(tmp_path, ["HG001"], run_generax=True)
    assert r.returncode == 0, r.stderr


@snakemake_only
def test_dryrun_with_generax_includes_generax_rule(tmp_path):
    r = _dry_run(tmp_path, ["HG001"], run_generax=True)
    assert r.returncode == 0, r.stderr
    assert "generax" in r.stdout.lower()


@snakemake_only
def test_dryrun_empty_ids_succeeds(tmp_path):
    r = _dry_run(tmp_path, [], run_generax=False)
    assert r.returncode == 0, r.stderr
