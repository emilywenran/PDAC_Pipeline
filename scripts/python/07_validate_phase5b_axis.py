#!/usr/bin/env python3
"""Validate Phase 5B continuous axis outputs."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "05_results" / "tables"
EXPECTED_ANALYSES = {
    "AXIS_MOFFITT50_PRIMARY",
    "AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY",
    "AXIS_SECONDARY",
    "AXIS_OUTLIER_EXCL",
    "AXIS_RAW_COUNTS",
    "AXIS_UNSUP_CENTROID",
    "AXIS_LEAVE_ONE_OUT",
}


def fail(msg: str) -> None:
    raise SystemExit(f"VALIDATION FAILED: {msg}")


def main() -> None:
    scores = pd.read_csv(TABLE / "phase5b_sample_continuous_scores.tsv", sep="\t")
    cent = pd.read_csv(TABLE / "phase5b_centroid_distance_scores.tsv", sep="\t")
    coverage = pd.read_csv(TABLE / "phase5b_signature_coverage.tsv", sep="\t")
    decision = pd.read_csv(TABLE / "phase5b_overall_decision.tsv", sep="\t")

    if scores["patient_id"].nunique() != 62:
        fail(f"expected 62 unique patients, observed {scores['patient_id'].nunique()}")

    dup = scores.duplicated(["patient_id", "expression_sample_id", "analysis_id", "scoring_method"])
    if dup.any():
        fail(f"duplicated sample-method-analysis records: {int(dup.sum())}")

    analyses = set(scores["analysis_id"].dropna().unique())
    if analyses != EXPECTED_ANALYSES:
        fail(f"analysis IDs mismatch: {sorted(analyses)}")

    if len(analyses) != 7:
        fail("expected exactly seven analysis IDs")

    if coverage["missing_genes"].sum() != 0:
        fail("signature coverage is incomplete without documented zero missing genes")

    numeric_cols = scores.select_dtypes(include=[np.number]).columns
    if np.isinf(scores[numeric_cols].to_numpy()).any():
        fail("infinite score values detected")

    loo = cent[cent["analysis_id"].eq("AXIS_LEAVE_ONE_OUT")]
    if loo.empty:
        fail("leave-one-out centroid rows missing")
    if loo["sample_included_in_centroid"].astype(str).str.lower().isin(["true", "1"]).any():
        fail("centroid leakage detected in leave-one-out analysis")

    primary = scores[scores["analysis_id"].eq("AXIS_MOFFITT50_PRIMARY")]
    med = primary.groupby("public_subtype")["basal_classical_contrast"].median()
    if not (med["Basal"] > med["Hybrid"] > med["Classical"]):
        fail("primary contrast score direction is not Basal > Hybrid > Classical")

    valid_decisions = {
        "TWO_POLES_WITH_INTERMEDIATE_CONTINUUM",
        "TWO_POLES_WITH_COACTIVATED_HYBRID",
        "HETEROGENEOUS_HYBRID_STATES",
        "NO_CLEAR_CONTINUOUS_AXIS",
        "INCONCLUSIVE",
    }
    if decision.loc[0, "overall_decision"] not in valid_decisions:
        fail("overall decision is not one of the locked categories")

    prohibited = [c for c in scores.columns if "optimized" in c.lower() or "cutoff" in c.lower()]
    if prohibited:
        fail(f"public-label-driven threshold optimization columns detected: {prohibited}")

    required_files = [
        "phase5b_centroid_distance_scores.tsv",
        "phase5b_score_method_concordance.tsv",
        "phase5b_method_sensitive_samples.tsv",
        "phase5b_public_group_score_comparison.tsv",
        "phase5b_ordered_trend_tests.tsv",
        "phase5b_hybrid_state_assessment.tsv",
        "phase5b_axis_distribution_tests.tsv",
        "phase5b_axis_stability_relationships.tsv",
        "phase5b_sensitivity_summary.tsv",
        "phase5b_category_transition_summary.tsv",
    ]
    missing = [f for f in required_files if not (TABLE / f).exists()]
    if missing:
        fail(f"missing required output files: {missing}")

    print("Phase 5B validation passed.")


if __name__ == "__main__":
    main()
