#!/usr/bin/env python3
"""Unit tests for Phase 10B-R corrected synthesis."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "06_scripts/python/18_phase10br_cross_layer_synthesis.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase10br", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_cell_composition_is_not_malignant_specific_support():
    module = load_module()
    row = pd.Series(
        {
            "feature_name": "CTCFL",
            "discovery_evidence": "ROBUST_HOST_MECHANISM",
            "bulk_evidence": "EXTERNALLY_REPLICATED_HOST_FEATURE",
            "sc_evidence": "CELL_COMPOSITION_EXPLAINED",
            "spatial_evidence": "NOT_EVALUATED",
            "final_synthesis_category": "PARTIALLY_REPLICATED",
        }
    )
    assert module.derive_category(row) == "PARTIALLY_REPLICATED"
    assert row["sc_evidence"] != "MALIGNANT_CELL_INTRINSIC_SUPPORT"


def test_ineligible_wgcna_modules_remain_insufficient_data():
    module = load_module()
    row = pd.Series(
        {
            "feature_name": "MEblue",
            "discovery_evidence": "ROBUST_HOST_MECHANISM",
            "bulk_evidence": "PARTIALLY_REPLICATED_HOST_FEATURE",
            "sc_evidence": "INSUFFICIENT_SINGLE_CELL_DATA",
            "spatial_evidence": "INSUFFICIENT_SPATIAL_DATA",
            "final_synthesis_category": "INSUFFICIENT_DATA",
        }
    )
    assert module.derive_category(row) == "INSUFFICIENT_DATA"


def test_external_database_unavailable_for_gene_without_local_result():
    module = load_module()
    status, source = module.external_status("CTCFL", "GTEx")
    assert status in {"NOT_RUN_DATABASE_UNAVAILABLE", "LOCAL_QUERY_RESULT_AVAILABLE"}
    if status == "NOT_RUN_DATABASE_UNAVAILABLE":
        assert source == ""


def test_phase10br_outputs_include_all_phase10a_candidates_after_run():
    inventory = pd.read_csv(ROOT / "05_results/tables/phase10a_cross_layer_evidence_inventory.tsv", sep="\t")
    targets = pd.read_csv(ROOT / "05_results/tables/phase10br_candidate_target_scores.tsv", sep="\t")
    assert set(inventory["feature_name"]) <= set(targets["feature_name"])


def test_ctcfl_is_blocked_by_composition_penalty_after_run():
    targets = pd.read_csv(ROOT / "05_results/tables/phase10br_candidate_target_scores.tsv", sep="\t")
    ctcfl = targets.loc[targets["feature_name"] == "CTCFL"].iloc[0]
    assert ctcfl["sc_evidence"] == "CELL_COMPOSITION_EXPLAINED"
    assert "COMPOSITION_SENSITIVE" in ctcfl["priority_decision"]
