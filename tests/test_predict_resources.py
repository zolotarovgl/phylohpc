"""Tests for workflow/predict_resources.py"""

import os
import sys
import math
import textwrap
import pytest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import predict_resources as pr


# ---------------------------------------------------------------------------
# Helper: write a small FASTA to a tmp file
# ---------------------------------------------------------------------------

def write_fasta(path, sequences):
    """sequences: list of (header, seq) tuples"""
    with open(path, "w") as fh:
        for hdr, seq in sequences:
            fh.write(f">{hdr}\n{seq}\n")


# ---------------------------------------------------------------------------
# get_nseq
# ---------------------------------------------------------------------------

def test_get_nseq_basic(tmp_path):
    f = tmp_path / "a.fasta"
    write_fasta(str(f), [("s1", "ACGT"), ("s2", "TTTT"), ("s3", "AAAA")])
    assert pr.get_nseq(str(f)) == 3


def test_get_nseq_single(tmp_path):
    f = tmp_path / "a.fasta"
    write_fasta(str(f), [("only", "ACGT")])
    assert pr.get_nseq(str(f)) == 1


def test_get_nseq_empty_file(tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    assert pr.get_nseq(str(f)) == 0


# ---------------------------------------------------------------------------
# get_mlen
# ---------------------------------------------------------------------------

def test_get_mlen_returns_median(tmp_path):
    f = tmp_path / "a.fasta"
    # lengths: 2, 4, 6  → median = 4
    write_fasta(str(f), [("s1", "AC"), ("s2", "ACGT"), ("s3", "ACGTAC")])
    assert pr.get_mlen(str(f)) == 4


def test_get_mlen_multiline_sequence(tmp_path):
    f = tmp_path / "ml.fasta"
    f.write_text(">seq1\nACGT\nACGT\n")  # total length = 8
    assert pr.get_mlen(str(f)) == 8


def test_get_mlen_empty_file(tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    assert pr.get_mlen(str(f)) == 0


def test_get_mlen_single_sequence(tmp_path):
    f = tmp_path / "single.fasta"
    write_fasta(str(f), [("s1", "AAAAAAAAAA")])  # length 10
    assert pr.get_mlen(str(f)) == 10


# ---------------------------------------------------------------------------
# predict_model
# ---------------------------------------------------------------------------

def test_predict_model_matches_manual_calculation():
    coefs = {
        "(Intercept)": 1.0,
        "log(nseq)": 0.5,
        "log(mlen)": 0.3,
        "log(nseq):log(mlen)": 0.1,
    }
    nseq, mlen = 10, 100
    ln_n = math.log(nseq)
    ln_m = math.log(mlen)
    expected = math.exp(1.0 + 0.5 * ln_n + 0.3 * ln_m + 0.1 * ln_n * ln_m)
    result = pr.predict_model(coefs, nseq, mlen)
    assert abs(result - expected) < 1e-9


def test_predict_model_positive_output():
    coefs = {
        "(Intercept)": 2.0,
        "log(nseq)": 1.0,
        "log(mlen)": 1.0,
        "log(nseq):log(mlen)": 0.5,
    }
    assert pr.predict_model(coefs, 50, 300) > 0


# ---------------------------------------------------------------------------
# round_base
# ---------------------------------------------------------------------------

def test_round_base_above_base():
    # 1100 → rounded up to next multiple of 1024 = 2048
    result = pr.round_base(np.array([1100.0]), 1024)
    assert result[0] == 2048.0


def test_round_base_below_base_unchanged():
    # 500 < 1024 → stays 500
    result = pr.round_base(np.array([500.0]), 1024)
    assert result[0] == 500.0


def test_round_base_exact_multiple():
    result = pr.round_base(np.array([2048.0]), 1024)
    assert result[0] == 2048.0


def test_round_base_vectorised():
    arr = np.array([70.0, 130.0])
    result = pr.round_base(arr, 60)
    # 70 ≥ 60 → ceil(70/60)*60 = 2*60 = 120
    # 130 ≥ 60 → ceil(130/60)*60 = 3*60 = 180
    assert result[0] == 120.0
    assert result[1] == 180.0


# ---------------------------------------------------------------------------
# convert_mem
# ---------------------------------------------------------------------------

def test_convert_mem_rounds_up_to_100():
    result = pr.convert_mem(np.array([150.0]))
    assert result[0] == "200.MB"


def test_convert_mem_exact_100():
    result = pr.convert_mem(np.array([100.0]))
    assert result[0] == "100.MB"


def test_convert_mem_suffix():
    result = pr.convert_mem(np.array([350.0]))
    assert result[0].endswith(".MB")


# ---------------------------------------------------------------------------
# convert_time
# ---------------------------------------------------------------------------

def test_convert_time_respects_minimum():
    # value < min_time=5 → should be bumped to 5
    result = pr.convert_time(np.array([2.0]))
    assert result[0] == "5.min"


def test_convert_time_normal_value():
    result = pr.convert_time(np.array([7.3]))
    assert result[0] == "8.min"


def test_convert_time_suffix():
    result = pr.convert_time(np.array([10.0]))
    assert result[0].endswith(".min")


def test_convert_time_custom_min():
    result = pr.convert_time(np.array([3.0]), min_time=10)
    assert result[0] == "10.min"


# ---------------------------------------------------------------------------
# predict_resources
# ---------------------------------------------------------------------------

DEFAULTS = {
    "ALIGN": {
        "mem": 500,
        "time": 30,
        "large_nseq_threshold": 1000,
        "large_mem": 8000,
        "large_time": 300,
    }
}

COEFS = {
    "(Intercept)": 5.0,
    "log(nseq)": 0.8,
    "log(mlen)": 0.6,
    "log(nseq):log(mlen)": 0.1,
}

MODELS_WITH_BOTH = {"ALIGN": {"mem": COEFS, "time": COEFS}}
MODELS_EMPTY = {}


def test_predict_resources_uses_model_when_available():
    mem, time = pr.predict_resources("ALIGN", MODELS_WITH_BOTH, DEFAULTS, 50, 300)
    expected_mem = pr.predict_model(COEFS, 50, 300)
    expected_time = pr.predict_model(COEFS, 50, 300)
    assert mem == max(expected_mem, DEFAULTS["ALIGN"]["mem"])
    assert time == max(expected_time, DEFAULTS["ALIGN"]["time"])


def test_predict_resources_falls_back_to_default_when_no_model():
    mem, time = pr.predict_resources("ALIGN", MODELS_EMPTY, DEFAULTS, 5, 50)
    assert mem == DEFAULTS["ALIGN"]["mem"]
    assert time == DEFAULTS["ALIGN"]["time"]


def test_predict_resources_large_family_override():
    # nseq = 2000, which is >= threshold of 1000
    mem, time = pr.predict_resources("ALIGN", MODELS_EMPTY, DEFAULTS, 2000, 300)
    assert mem == DEFAULTS["ALIGN"]["large_mem"]
    assert time == DEFAULTS["ALIGN"]["large_time"]


def test_predict_resources_default_is_minimum():
    # Force model to predict very small values; defaults must win.
    tiny_coefs = {
        "(Intercept)": -10.0,
        "log(nseq)": 0.0,
        "log(mlen)": 0.0,
        "log(nseq):log(mlen)": 0.0,
    }
    models = {"ALIGN": {"mem": tiny_coefs, "time": tiny_coefs}}
    mem, time = pr.predict_resources("ALIGN", models, DEFAULTS, 5, 50)
    assert mem == DEFAULTS["ALIGN"]["mem"]
    assert time == DEFAULTS["ALIGN"]["time"]
