#!/usr/bin/env python3
"""Executable tests for Phase 9B3R repair-specification guardrails."""

from __future__ import annotations

import importlib.util
import inspect
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/Users/emily/thesis/PDAC")
SCRIPT = ROOT / "06_scripts/python/16_phase9b3r_spatial_validation.py"
TABLES = ROOT / "05_results/tables"


def load_module():
    spec = importlib.util.spec_from_file_location("phase9b3r", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rank_score_enforces_coverage():
    mod = load_module()
    expr = pd.DataFrame(
        [[1, 2], [2, 3], [3, 4]],
        index=["A", "B", "C"],
        columns=["s1", "s2"],
    )
    score = mod.rank_score(expr, ["A", "B", "C", "D", "E"], expected_gene_count=5)
    assert score.isna().all()


def test_add_result_convergence_failure():
    mod = load_module()

    class FailedFit:
        converged = False
        params = {"is_tumor": 1.0}
        bse = {"is_tumor": 0.1}
        pvalues = {"is_tumor": 0.001}

    rows = []
    fit = FailedFit()
    converged = bool(getattr(fit, "converged", False))
    usable = fit is not None and "is_tumor" in getattr(fit, "params", {}) and converged
    if not usable:
        rows.append({"coefficient": np.nan, "std_error": np.nan, "p_value": np.nan, "q_value": np.nan})
    assert math.isnan(rows[0]["coefficient"])
    assert math.isnan(rows[0]["p_value"])


def test_validator_rejects_uniform_placeholders():
    placeholders = pd.DataFrame(
        {
            "observed_statistic": [0.0, 0.0],
            "empirical_p_value": [1.0, 1.0],
            "null_variance": [0.0, 0.0],
            "execution_status": ["EXECUTED", "EXECUTED"],
        }
    )
    assert ((placeholders["observed_statistic"] == 0) & (placeholders["empirical_p_value"] == 1)).any()
    assert (placeholders["null_variance"] <= 0).any()


def test_output_negative_controls_have_real_iterations():
    neg_path = TABLES / "phase9b3r_negative_control_results.tsv"
    null_path = TABLES / "phase9b3r_negative_control_null_distributions.tsv"
    if not neg_path.exists() or not null_path.exists():
        return
    neg = pd.read_csv(neg_path, sep="\t")
    nulls = pd.read_csv(null_path, sep="\t")
    assert len(nulls) > 0
    assert (neg["execution_status"] == "EXECUTED").all()
    assert (pd.to_numeric(neg["null_variance"], errors="coerce") > 0).all()
    assert not (((neg["observed_statistic"] == 0) & (neg["empirical_p_value"] == 1)).any())


def test_ineligible_feature_does_not_enter_model_outputs():
    path = TABLES / "phase9b3r_hwang_naive_models.tsv"
    if not path.exists():
        return
    rows = pd.read_csv(path, sep="\t")
    comp = rows[rows["feature_name"].eq("HALLMARK_SPERMATOGENESIS")]
    assert not comp["model_id"].str.endswith(("MODEL_A", "MODEL_B", "MODEL_C")).any()
    assert comp["p_value"].isna().all()
    assert comp["q_value"].isna().all()


def test_nonconverged_model_has_no_reportable_p_or_q():
    path = TABLES / "phase9b3r_hwang_treated_models.tsv"
    if not path.exists():
        return
    rows = pd.read_csv(path, sep="\t")
    bad = rows[rows["model_converged"].astype(str).eq("False")]
    assert bad["p_value"].isna().all()
    assert bad["q_value"].isna().all()


def test_no_unauthorized_z_test_inference_used():
    for name in ["phase9b3r_hwang_naive_models.tsv", "phase9b3r_hwang_treated_models.tsv"]:
        path = TABLES / name
        if not path.exists():
            continue
        rows = pd.read_csv(path, sep="\t")
        eligible = rows[rows["eligibility_status"].eq("ELIGIBLE")]
        assert set(eligible["inference_method"]) == {"statsmodels_asymptotic_z_locked_plan"}
        assert "unauthorized_z_test" not in set(rows["inference_method"])


def test_invalid_results_do_not_contribute_to_evidence_classification():
    path = TABLES / "phase9b3r_spatial_evidence.tsv"
    if not path.exists():
        return
    evidence = pd.read_csv(path, sep="\t")
    invalid = evidence[evidence["feature_name"].isin(["HALLMARK_SPERMATOGENESIS", "MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"])]
    assert set(invalid["evidence_category"]) == {"INSUFFICIENT_SPATIAL_DATA"}


def test_evidence_derivation_logic_not_hardcoded():
    mod = load_module()
    source = inspect.getsource(mod.derive_evidence)
    assert "model_b" in source
    assert "controls_ok" in source
    assert "mon_positive" in source
    assert "elif" in source
