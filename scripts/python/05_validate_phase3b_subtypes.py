#!/usr/bin/env python3
"""Validate Phase 3B subtype reproduction outputs."""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results" / "tables"
SCRIPT = ROOT / "06_scripts/python/05_phase3b_reproduce_subtypes.py"

PRIMARY = TABLE_DIR / "phase3b_primary_subtype_assignments.tsv"
ALL_METHODS = TABLE_DIR / "phase3b_all_method_assignments.tsv"
RUNTIME = TABLE_DIR / "phase3b_signature_runtime_validation.tsv"
METRICS = TABLE_DIR / "phase3b_method_agreement_metrics.tsv"
CONFUSIONS = TABLE_DIR / "phase3b_confusion_matrices.tsv"
DISCORDANT = TABLE_DIR / "phase3b_discordant_samples.tsv"
SENSITIVITY = TABLE_DIR / "phase3b_sensitivity_summary.tsv"
REPORT = ROOT / "04_analysis/05_subtype_reproduction/PHASE3B_SUBTYPE_REPRODUCTION.md"

ALLOWED_PUBLIC = {"Basal", "Hybrid", "Classical"}
ALLOWED_BY_METHOD = {
    "GSE172356_original": {"Basal", "Hybrid", "Classical"},
    "Moffitt": {"Basal", "Classical", "Others"},
    "PurIST": {"Basal-like", "Classical"},
}
EXPECTED_METHODS = set(ALLOWED_BY_METHOD)
EXPECTED_GENES_USED = {"GSE172356_original": 94, "Moffitt": 49, "PurIST": 16}
EXPECTED_PATIENTS = 62


def fail(errors: list[str], message: str):
    errors.append(message)


def require_file(errors: list[str], path: Path):
    if not path.exists():
        fail(errors, f"Missing required file: {path.relative_to(ROOT)}")


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def validate_locked_constants(errors: list[str]):
    text = SCRIPT.read_text()
    expected_snippets = [
        "PRIMARY_SLICE_SIZES = [17, 23, 22]",
        'PRIMARY_SLICE_LABELS = ["Basal", "Hybrid", "Classical"]',
        "MOFFITT_SLICE_SIZES = [27, 17, 18]",
        'MOFFITT_SLICE_LABELS = ["Classical", "Basal", "Others"]',
        "PURIST_INTERCEPT = -6.815",
        "PURIST_CUTOFF = 0.5",
        '("GPR87", "REG4", 1.994)',
        '("KRT6A", "ANXA10", 2.031)',
        '("BCAR3", "GATA6", 1.618)',
        '("PTGES", "CLDN18", 0.922)',
        '("ITGA3", "LGALS4", 1.059)',
        '("C16orf74", "DDC", 0.929)',
        '("S100A2", "SLC40A1", 2.505)',
        '("KRT5", "CLRN3", 0.485)',
    ]
    for snippet in expected_snippets:
        if snippet not in text:
            fail(errors, f"Locked parameter missing or changed in script: {snippet}")

    banned = ["GridSearchCV", "RandomizedSearchCV", "fit(", "train_test_split", "SMOTE", "SelectKBest", "feature_selection"]
    for token in banned:
        if token in text:
            fail(errors, f"Potential accidental training/optimization token found in Phase 3B script: {token}")


def main():
    errors: list[str] = []
    for path in [PRIMARY, ALL_METHODS, RUNTIME, METRICS, CONFUSIONS, DISCORDANT, SENSITIVITY, REPORT, SCRIPT]:
        require_file(errors, path)
    if errors:
        raise SystemExit("\n".join(errors))

    validate_locked_constants(errors)

    primary = read(PRIMARY)
    all_methods = read(ALL_METHODS)
    runtime = read(RUNTIME)
    metrics = read(METRICS)
    sensitivity = read(SENSITIVITY)
    report_text = REPORT.read_text()

    if primary["patient_id"].nunique() != EXPECTED_PATIENTS:
        fail(errors, f"Primary output has {primary['patient_id'].nunique()} unique patients, expected 62.")
    if len(primary) != EXPECTED_PATIENTS:
        fail(errors, f"Primary output has {len(primary)} rows, expected 62.")
    if primary["patient_id"].duplicated().any():
        fail(errors, "Primary output contains duplicated patient_id values.")
    if primary["expression_sample_id"].duplicated().any():
        fail(errors, "Primary output contains duplicated expression_sample_id values.")
    if set(primary["original_public_subtype"]) - ALLOWED_PUBLIC:
        fail(errors, "Primary output contains disallowed public subtype labels.")
    if set(primary["reproduced_subtype"]) - ALLOWED_PUBLIC:
        fail(errors, "Primary output contains disallowed reproduced subtype labels.")

    counts = primary["reproduced_subtype"].value_counts().to_dict()
    if counts.get("Basal") != 17 or counts.get("Hybrid") != 23 or counts.get("Classical") != 22:
        fail(errors, f"Primary reproduced class counts are not 17/23/22: {counts}")

    if all_methods["patient_id"].nunique() != EXPECTED_PATIENTS:
        fail(errors, "All-method output does not contain 62 unique patients.")
    observed_methods = set(all_methods["method_name"])
    if observed_methods != EXPECTED_METHODS:
        fail(errors, f"All-method methods are {observed_methods}, expected {EXPECTED_METHODS}.")

    per_method_counts = all_methods.groupby("method_name")["patient_id"].nunique().to_dict()
    for method in EXPECTED_METHODS:
        if per_method_counts.get(method) != EXPECTED_PATIENTS:
            fail(errors, f"{method} has {per_method_counts.get(method)} patients, expected 62.")
        duplicated = all_methods[all_methods["method_name"].eq(method)]["patient_id"].duplicated().any()
        if duplicated:
            fail(errors, f"{method} contains duplicated patient IDs.")

    for method, allowed in ALLOWED_BY_METHOD.items():
        observed = set(all_methods.loc[all_methods["method_name"].eq(method), "predicted_subtype"])
        if observed - allowed:
            fail(errors, f"{method} has disallowed labels: {observed - allowed}")
        genes_used = set(all_methods.loc[all_methods["method_name"].eq(method), "genes_used"].astype(int))
        if genes_used != {EXPECTED_GENES_USED[method]}:
            fail(errors, f"{method} genes_used values {genes_used}, expected {EXPECTED_GENES_USED[method]}.")

    if runtime["status"].ne("PASS").any():
        fail(errors, "Runtime validation table contains non-PASS rows.")
    if not runtime["validation_item"].eq("required_expression_scale").any():
        fail(errors, "Runtime validation does not document required expression scales.")
    if not runtime["validation_item"].eq("pair_direction_coefficients").any():
        fail(errors, "Runtime validation does not document PurIST gene-pair directions and coefficients.")

    exact = metrics[
        metrics["method_name"].eq("GSE172356_original")
        & metrics["analysis_set"].eq("full_62")
        & metrics["metric"].eq("exact_agreement")
    ]["value"]
    if exact.empty or abs(float(exact.iloc[0]) - 1.0) > 1e-12:
        fail(errors, "Primary full_62 exact agreement is not 1.0.")

    required_sensitivity = {
        "full_62",
        "exclude_phase2b_outlier_candidates",
        "log2_median_centering_stress_test",
        "alternative_missingness_gene_median_imputed",
        "alternative_missingness_zero_filled",
    }
    observed_sensitivity = set(sensitivity["analysis_set"])
    if not required_sensitivity.issubset(observed_sensitivity):
        fail(errors, f"Sensitivity summary missing analyses: {required_sensitivity - observed_sensitivity}")

    if re.search(r"\bGridSearchCV\b|\bRandomizedSearchCV\b|\bSMOTE\b|\bSelectKBest\b", report_text, flags=re.IGNORECASE):
        fail(errors, "Report contains language suggesting model tuning or feature selection tooling.")
    if "public labels were not used as model-training inputs" not in report_text:
        fail(errors, "Report does not explicitly state that public labels were not used as model-training inputs.")
    if "score direction" not in report_text.lower():
        fail(errors, "Report does not explicitly document score directions.")

    if errors:
        raise SystemExit("Phase 3B validation failed:\n" + "\n".join(f"- {e}" for e in errors))

    print("Phase 3B validation passed.")
    print("Validated 62 unique patients, 3 applicable methods, locked gene counts, locked parameters, and documented score directions.")


if __name__ == "__main__":
    main()
