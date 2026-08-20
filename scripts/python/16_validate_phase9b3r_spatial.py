#!/usr/bin/env python3
"""Validate Phase 9B3R corrected spatial validation outputs."""

from __future__ import annotations

import ast
import csv
import math
import sys
from pathlib import Path


ROOT = Path("/Users/emily/thesis/PDAC")
TABLES = ROOT / "05_results/tables"
FIGS = ROOT / "05_results/figures"
ANALYSIS = ROOT / "04_analysis/09_external_validation"
SCRIPT = ROOT / "06_scripts/python/16_phase9b3r_spatial_validation.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def evidence_is_hardcoded() -> bool:
    source = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(source)
    derive = next((n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "derive_evidence"), None)
    if derive is None:
        return True
    names = {n.id for n in ast.walk(derive) if isinstance(n, ast.Name)}
    required_inputs = {"coverage", "naive", "treated", "mon", "neg", "model_a", "model_b", "model_c", "controls_ok"}
    has_branching = any(isinstance(n, ast.If) for n in ast.walk(derive))
    return not (required_inputs <= names and has_branching)


def main() -> int:
    errors: list[str] = []
    required = [
        "phase9b3r_runtime_validation.tsv",
        "phase9b3r_dataset_qc.tsv",
        "phase9b3r_feature_coverage.tsv",
        "phase9b3r_hwang_naive_models.tsv",
        "phase9b3r_hwang_treated_models.tsv",
        "phase9b3r_moncada_exploratory_results.tsv",
        "phase9b3r_negative_control_results.tsv",
        "phase9b3r_negative_control_null_distributions.tsv",
        "phase9b3r_spatial_evidence.tsv",
        "phase9b3r_cross_cohort_synthesis.tsv",
    ]
    for name in required:
        path = TABLES / name
        if not path.exists() or path.stat().st_size == 0:
            fail(errors, f"Missing or empty table: {name}")
    for path in [
        ANALYSIS / "PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS.md",
        ANALYSIS / "PHASE9B3R_CORRECTION_LOG.md",
        FIGS / "phase9b3r_hwang_primary_models.pdf",
    ]:
        if not path.exists() or path.stat().st_size == 0:
            fail(errors, f"Missing or empty output: {path.relative_to(ROOT)}")
    if errors:
        print("Phase 9B3R validation FAILED")
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1

    coverage = read_tsv(TABLES / "phase9b3r_feature_coverage.tsv")
    ineligible = {
        (r["dataset_id"], r["feature_name"])
        for r in coverage
        if r.get("formal_inference_status") == "INSUFFICIENT_SPATIAL_DATA"
    }
    if ("HWANG_GSE202051_NAIVE", "HALLMARK_SPERMATOGENESIS") not in ineligible:
        fail(errors, "HALLMARK_SPERMATOGENESIS is not marked insufficient in naive cohort")
    for row in coverage:
        if row.get("feature_layer") == "WGCNA_module" and number(row.get("coverage_fraction", "")) < 0.80:
            if row.get("formal_inference_status") != "INSUFFICIENT_SPATIAL_DATA":
                fail(errors, f"Ineligible WGCNA module not gated: {row.get('feature_name')}")

    model_rows = []
    for table in ["phase9b3r_hwang_naive_models.tsv", "phase9b3r_hwang_treated_models.tsv"]:
        rows = read_tsv(TABLES / table)
        model_rows.extend(rows)
        for row in rows:
            key = (row.get("cohort_id"), row.get("feature_name"))
            if key in ineligible:
                if row.get("model_id", "").endswith(("MODEL_A", "MODEL_B", "MODEL_C")):
                    fail(errors, f"Ineligible feature entered formal model in {table}: {key}")
                if math.isfinite(number(row.get("p_value", ""))) or math.isfinite(number(row.get("q_value", ""))):
                    fail(errors, f"Ineligible feature retained p/q in {table}: {key}")
            if row.get("model_converged") == "False":
                for col in ["coefficient", "std_error", "ci_low", "ci_high", "p_value", "q_value"]:
                    if math.isfinite(number(row.get(col, ""))):
                        fail(errors, f"Nonconverged model retained {col}: {row.get('model_id')}")
            if row.get("inference_method") == "asymptotic_z" or row.get("inference_method") == "unauthorized_z_test":
                fail(errors, f"Unauthorized inference method label in {table}: {row.get('model_id')}")
            if row.get("eligibility_status") == "ELIGIBLE" and row.get("inference_method") != "statsmodels_asymptotic_z_locked_plan":
                fail(errors, f"Missing locked inference method label in {table}: {row.get('model_id')}")

    neg = read_tsv(TABLES / "phase9b3r_negative_control_results.tsv")
    nulls = read_tsv(TABLES / "phase9b3r_negative_control_null_distributions.tsv")
    required_controls = {
        "coordinate permutation",
        "size-matched random gene set",
        "expression-matched random gene set",
        "unrelated Hallmark pathway",
        "label permutation",
        "leakage control",
    }
    observed_controls = {r.get("control_type") for r in neg if r.get("execution_status") == "EXECUTED"}
    missing = required_controls - observed_controls
    if missing:
        fail(errors, f"Missing executed negative controls: {sorted(missing)}")
    if not nulls:
        fail(errors, "No iteration-level null-distribution rows")
    null_counts: dict[tuple[str, str, str], int] = {}
    for row in nulls:
        key = (row.get("dataset_id", ""), row.get("control_type", ""), row.get("control_id", ""))
        null_counts[key] = null_counts.get(key, 0) + 1
    for row in neg:
        if row.get("execution_status") != "EXECUTED":
            fail(errors, f"Negative control not executed: {row}")
            continue
        iterations = int(number(row.get("iterations", "")))
        key = (row.get("dataset_id", ""), row.get("control_type", ""), row.get("control_id", ""))
        if null_counts.get(key, 0) < iterations:
            fail(errors, f"Negative control lacks iteration rows: {key}")
        if number(row.get("null_variance", "")) <= 0:
            fail(errors, f"Negative control has zero/null variance: {key}")
        if number(row.get("observed_statistic", "")) == 0 and number(row.get("empirical_p_value", "")) == 1:
            fail(errors, f"Placeholder-like negative control accepted: {key}")

    evidence = read_tsv(TABLES / "phase9b3r_spatial_evidence.tsv")
    ev_by_feature = {r["feature_name"]: r for r in evidence}
    if ev_by_feature.get("HALLMARK_SPERMATOGENESIS", {}).get("evidence_category") != "INSUFFICIENT_SPATIAL_DATA":
        fail(errors, "Ineligible comparator contributes as biological evidence")
    for module in ["MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"]:
        if ev_by_feature.get(module, {}).get("evidence_category") != "INSUFFICIENT_SPATIAL_DATA":
            fail(errors, f"Ineligible WGCNA module contributes to evidence: {module}")
    if evidence_is_hardcoded():
        fail(errors, "Evidence category appears hardcoded instead of derived")

    text = (ANALYSIS / "PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS.md").read_text(encoding="utf-8")
    if "READY_FOR_PHASE9B3C2_COMPLETE_INDEPENDENT_REVIEW" not in text:
        fail(errors, "Final readiness decision missing from corrected report")
    if "SUPERSEDED_BY_PHASE9B3R" not in (ANALYSIS / "PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md").read_text(encoding="utf-8"):
        fail(errors, "Original Phase 9B3B report is not marked superseded")

    if errors:
        print("Phase 9B3R validation FAILED")
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1
    print("Phase 9B3R validation PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
