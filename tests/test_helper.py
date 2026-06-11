"""Tests for workflow/helper.py"""

import json
import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import helper


# ---------------------------------------------------------------------------
# parse_bash_config
# ---------------------------------------------------------------------------

def test_parse_bash_config_basic(tmp_path):
    cfg_file = tmp_path / "config.sh"
    cfg_file.write_text("FOO=bar\nBAZ=123\n")
    result = helper.parse_bash_config(str(cfg_file))
    assert result == {"FOO": "bar", "BAZ": "123"}


def test_parse_bash_config_ignores_comments_and_blank_lines(tmp_path):
    cfg_file = tmp_path / "config.sh"
    cfg_file.write_text("# comment\n\nKEY=value\n")
    result = helper.parse_bash_config(str(cfg_file))
    assert result == {"KEY": "value"}


def test_parse_bash_config_value_with_equals(tmp_path):
    """Values containing '=' should be kept intact (split on first '=' only)."""
    cfg_file = tmp_path / "config.sh"
    cfg_file.write_text("URL=http://example.com?a=1\n")
    result = helper.parse_bash_config(str(cfg_file))
    assert result == {"URL": "http://example.com?a=1"}


def test_parse_bash_config_empty_file(tmp_path):
    cfg_file = tmp_path / "config.sh"
    cfg_file.write_text("")
    result = helper.parse_bash_config(str(cfg_file))
    assert result == {}


# ---------------------------------------------------------------------------
# parse_genefam
# ---------------------------------------------------------------------------

GENEFAM_LINE = "MYFAM\tHMM1,HMM2\t2.5\t4\tga\tgroupA\tpfx"


def test_parse_genefam_basic(tmp_path):
    gf_file = tmp_path / "gf.tsv"
    gf_file.write_text(GENEFAM_LINE + "\n")
    result = helper.parse_genefam(str(gf_file))
    assert "pfx.MYFAM" in result
    rec = result["pfx.MYFAM"]
    assert rec["family"] == "MYFAM"
    assert rec["hmms"] == ["HMM1", "HMM2"]
    assert rec["min_seq"] == 4
    assert rec["threshold"] == "ga"
    assert rec["group"] == "groupA"
    assert rec["prefix"] == "pfx"


def test_parse_genefam_no_prefix_append(tmp_path):
    gf_file = tmp_path / "gf.tsv"
    gf_file.write_text(GENEFAM_LINE + "\n")
    result = helper.parse_genefam(str(gf_file), append_prefix=False)
    assert "MYFAM" in result
    assert "pfx.MYFAM" not in result


def test_parse_genefam_multiple_families(tmp_path):
    lines = "\n".join([
        "FAM1\tHMM1\t1.4\t3\tga\tgA\tpA",
        "FAM2\tHMM2\t2.0\t5\ttc\tgB\tpB",
    ])
    gf_file = tmp_path / "gf.tsv"
    gf_file.write_text(lines + "\n")
    result = helper.parse_genefam(str(gf_file))
    assert len(result) == 2
    assert "pA.FAM1" in result
    assert "pB.FAM2" in result


def test_parse_genefam_bad_line_skipped(tmp_path, capsys):
    """Lines with wrong number of fields should trigger a warning but not crash."""
    gf_file = tmp_path / "gf.tsv"
    gf_file.write_text("BADLINE\tonly_two_fields\n" + GENEFAM_LINE + "\n")
    result = helper.parse_genefam(str(gf_file))
    captured = capsys.readouterr()
    assert "Unknown genefam format" in captured.out
    assert len(result) == 1  # only the valid line was parsed


# ---------------------------------------------------------------------------
# read_json
# ---------------------------------------------------------------------------

def test_read_json_valid(tmp_path):
    jf = tmp_path / "data.json"
    jf.write_text(json.dumps({"key": "val"}))
    assert helper.read_json(str(jf)) == {"key": "val"}


def test_read_json_missing_file(tmp_path, capsys):
    result = helper.read_json(str(tmp_path / "nonexistent.json"))
    assert result == {}
    assert "WARNING" in capsys.readouterr().out


def test_read_json_empty_file(tmp_path, capsys):
    jf = tmp_path / "empty.json"
    jf.write_text("")
    result = helper.read_json(str(jf))
    assert result == {}
    assert "WARNING" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# parse_mem  (inner function of check_job – tested via the module directly)
# ---------------------------------------------------------------------------
# The parse_mem logic is duplicated in check_job.py as a standalone function.
# We expose the same logic via helper.check_job's inner scope, but it is
# easier to test the identical standalone version in check_job.py.
# Here we at least verify the helper.check_job returns the UNKNOWN sentinel
# when called with a non-existent job ID (i.e., when squeue/sacct are absent
# or return nothing). This guards against regressions in the output structure.

def test_check_job_returns_dict_with_expected_keys(monkeypatch):
    import subprocess

    def fake_run(cmd, **kwargs):
        result = subprocess.CompletedProcess(cmd, 0)
        result.stdout = ""
        result.stderr = ""
        return result

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = helper.check_job("99999")
    assert set(result.keys()) == {
        "jobid", "state", "exitcode", "elapsed", "timelimit", "reqmem", "maxrss"
    }
    assert result["state"] == "UNKNOWN"
