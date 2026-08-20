#!/usr/bin/env python3
"""Validate corrected Phase 9B2R single-cell result substance."""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results/tables"
FIG_DIR = ROOT / "05_results/figures"
ANALYSIS_DIR = ROOT / "04_analysis/09_external_validation"

MODULES = {"MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"}
PLACEHOLDER_STATUSES = {"PASS_DESCRIPTIVE_CONTROL_GENERATED", "PLACEHOLDER", "TO_VERIFY", "TO_VERIFY_FOR_SOME_CONTROLS"}
REQUIRED_TABLES = {
    "phase9b2r_runtime_validation.tsv",
    "phase9b2r_module_transfer_coverage.tsv",
    "phase9b2r_negative_control_results.tsv",
    "phase9b2r_tf_evidence_classification.tsv",
    "phase9b2r_malignant_feature_axis_associations.tsv",
    "phase9b2r_cellular_source_evidence.tsv",
    "phase9b2r_cellular_source_models.tsv",
    "phase9b2r_cell_composition_sensitivity.tsv",
    "phase9b2r_tumor_control_descriptive.tsv",
    "phase9b2r_core_analysis_verification.tsv",
}
REQUIRED_FIGURES = {
    "phase9b2r_cohort_cell_counts.pdf",
    "phase9b2r_cell_annotation_markers.pdf",
    "phase9b2r_malignant_cell_audit.pdf",
    "phase9b2r_moffitt_axis_by_cell_type.pdf",
    "phase9b2r_malignant_axis_by_patient.pdf",
    "phase9b2r_hallmark_cellular_source.pdf",
    "phase9b2r_tf_activity_cellular_source.pdf",
    "phase9b2r_malignant_feature_axis_heatmap.pdf",
    "phase9b2r_cell_composition_sensitivity.pdf",
    "phase9b2r_tumor_control_descriptive.pdf",
    "phase9b2r_negative_control_summary.pdf",
    "phase9b2r_cellular_source_evidence_summary.pdf",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_missing(value: str | None) -> bool:
    return value is None or value == "" or value.upper() in {"NA", "NAN", "NULL"}


def as_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def main() -> int:
    errors: list[str] = []

    for name in REQUIRED_TABLES:
        path = TABLE_DIR / name
        if not path.exists() or path.stat().st_size == 0:
            errors.append(f"Missing or empty required table: {path.relative_to(ROOT)}")
    for name in REQUIRED_FIGURES:
        path = FIG_DIR / name
        if not path.exists() or path.stat().st_size == 0:
            errors.append(f"Missing or empty required figure: {path.relative_to(ROOT)}")

    old_report = ANALYSIS_DIR / "PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md"
    if old_report.exists() and "SUPERSEDED_BY_PHASE9B2R" not in old_report.read_text(encoding="utf-8", errors="ignore")[:500]:
        errors.append("Superseded Phase 9B2 report is not marked SUPERSEDED_BY_PHASE9B2R near the top.")

    report = ANALYSIS_DIR / "PHASE9B2R_CORRECTED_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md"
    if not report.exists() or report.stat().st_size == 0:
        errors.append("Missing corrected Phase 9B2R report.")
    else:
        text = report.read_text(encoding="utf-8", errors="ignore")
        for required in ["Ochrobactrum was not tested", "PENG_CRA001160", "INSUFFICIENT_SINGLE_CELL_DATA", "ready for complete independent review"]:
            if required not in text:
                errors.append(f"Corrected report missing required statement: {required}")

    if errors:
        print("Phase 9B2R validation FAILED:", file=sys.stderr)
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1

    cov = read_tsv(TABLE_DIR / "phase9b2r_module_transfer_coverage.tsv")
    for row in cov:
        module = row["module_name"]
        coverage = as_float(row["coverage_fraction"])
        if module in MODULES and coverage < 0.80 and row.get("eligibility") == "ELIGIBLE":
            errors.append(f"Low-coverage module {module} is marked eligible.")
        for col in ("total_discovery_genes", "mapped_genes", "detected_genes", "duplicate_mappings"):
            if is_missing(row.get(col)):
                errors.append(f"Module coverage missing {col} for {module}.")

    assoc = read_tsv(TABLE_DIR / "phase9b2r_malignant_feature_axis_associations.tsv")
    assoc_features = {r["feature_name"] for r in assoc}
    bad_modules = MODULES & assoc_features
    if bad_modules:
        errors.append("Ineligible modules appear in corrected malignant-axis associations: " + ", ".join(sorted(bad_modules)))

    models = read_tsv(TABLE_DIR / "phase9b2r_cellular_source_models.tsv")
    model_features = {r["feature_name"] for r in models}
    bad_model_modules = MODULES & model_features
    if bad_model_modules:
        errors.append("Ineligible modules appear in corrected cellular-source models: " + ", ".join(sorted(bad_model_modules)))

    neg = read_tsv(TABLE_DIR / "phase9b2r_negative_control_results.tsv")
    if not neg:
        errors.append("Negative-control table is empty.")
    required_controls = {
        "size-matched randomized gene sets",
        "expression-matched randomized gene sets",
        "unrelated Hallmark pathway controls",
        "nonselected TF regulon controls",
        "patient-label permutation",
        "cell-type-label permutation",
    }
    observed_controls = {r["control_type"] for r in neg}
    missing_controls = required_controls - observed_controls
    if missing_controls:
        errors.append("Missing required negative-control types: " + ", ".join(sorted(missing_controls)))
    for row in neg:
        status = row.get("execution_status", "")
        if status in PLACEHOLDER_STATUSES:
            errors.append(f"Placeholder negative-control status for {row.get('target_feature')}: {status}")
        if status not in {"EXECUTED", "TECHNICALLY_INAPPLICABLE", "FAILED_WITH_REASON"}:
            errors.append(f"Invalid negative-control execution status: {status}")
        if row["control_type"] in {"patient-label permutation", "cell-type-label permutation"}:
            if row.get("iteration_count") != "1000":
                errors.append(f"Permutation iteration count is not 1000 for {row.get('target_feature')}.")
        if row["control_type"] in {"size-matched randomized gene sets", "expression-matched randomized gene sets"}:
            if row.get("iteration_count") != "100":
                errors.append(f"Module random-control iteration count is not 100 for {row.get('target_feature')}.")
        if status == "EXECUTED":
            for col in ("empirical_p_value", "candidate_statistic", "control_statistic"):
                if is_missing(row.get(col)):
                    errors.append(f"Executed control missing {col}: {row.get('control_type')} {row.get('target_feature')}")
        if status == "TECHNICALLY_INAPPLICABLE" and not row.get("failure_reason"):
            errors.append(f"Technically inapplicable control lacks reason: {row.get('target_feature')}")

    evidence = read_tsv(TABLE_DIR / "phase9b2r_cellular_source_evidence.tsv")
    ev_by_feature = {r["feature_name"]: r for r in evidence}
    for module in MODULES:
        row = ev_by_feature.get(module)
        if not row:
            errors.append(f"Module missing from evidence table: {module}")
        elif row.get("final_category") != "INSUFFICIENT_SINGLE_CELL_DATA":
            errors.append(f"Low-coverage module {module} has invalid category {row.get('final_category')}.")
    bad_supported = [r["feature_name"] for r in evidence if r["feature_name"] in MODULES and "SUPPORT" in r.get("final_category", "")]
    if bad_supported:
        errors.append("Low-coverage modules classified as supported: " + ", ".join(sorted(bad_supported)))

    tf = read_tsv(TABLE_DIR / "phase9b2r_tf_evidence_classification.tsv")
    categories = {r.get("evidence_category", "") for r in tf}
    if len(categories) <= 1:
        errors.append("All TFs have one blanket evidence category.")
    for row in tf:
        if row.get("activity_calculation_status") == "EXECUTED" and row.get("evidence_category") == "TO_VERIFY":
            errors.append(f"Executed TF remains TO_VERIFY without rule-based category: {row.get('feature_name')}")
        method_text = " ".join(row.values()).lower()
        if "expression proxy" in method_text or "tf_symbol_expression" in method_text:
            errors.append(f"TF expression proxy detected for {row.get('feature_name')}")

    for fig in ("phase9b2r_malignant_feature_axis_heatmap.pdf",):
        # Substance is checked through the source table used for this figure.
        if bad_modules:
            errors.append(f"Ineligible feature would appear in corrected formal result figure {fig}.")

    if errors:
        print("Phase 9B2R validation FAILED:", file=sys.stderr)
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1

    print("Phase 9B2R validation passed.")
    print("canonical_dataset_id=PENG_CRA001160")
    print("result=READY_FOR_COMPLETE_INDEPENDENT_REVIEW")
    return 0


if __name__ == "__main__":
    sys.exit(main())
