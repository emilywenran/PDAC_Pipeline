#!/usr/bin/env python3
"""Validate Phase 9B3B spatial validation outputs and guardrails."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
TABLES = ROOT / "05_results/tables"
REPORT = ROOT / "04_analysis/09_external_validation/PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fail(msg: str, errors: list[str]) -> None:
    errors.append(msg)


def main() -> int:
    errors: list[str] = []
    required = [
        "phase9b3b_runtime_validation.tsv",
        "phase9b3b_dataset_qc.tsv",
        "phase9b3b_feature_coverage.tsv",
        "phase9b3b_hwang_naive_models.tsv",
        "phase9b3b_hwang_treated_models.tsv",
        "phase9b3b_moncada_exploratory_results.tsv",
        "phase9b3b_negative_control_results.tsv",
        "phase9b3b_spatial_evidence.tsv",
        "phase9b3b_cross_cohort_synthesis.tsv",
    ]
    for name in required:
        path = TABLES / name
        if not path.exists() or path.stat().st_size == 0:
            fail(f"Missing or empty required table: {name}", errors)
    if not REPORT.exists():
        fail("Missing PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md", errors)
    if errors:
        print("Phase 9B3B validation FAILED")
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1

    runtime = read_tsv(TABLES / "phase9b3b_runtime_validation.tsv")
    flags = {r["check_id"]: r["status"] for r in runtime}
    for check in [
        "roi_pairing_preserved",
        "patient_replicate_preserved",
        "hwang_cohorts_not_pooled",
        "moncada_exploratory_only",
        "platform_matrices_not_merged",
        "target_genes_not_used_for_compartment_assignment",
        "no_microbiome_or_causality_claims",
    ]:
        if flags.get(check) != "PASS":
            fail(f"Runtime guardrail failed or missing: {check}", errors)

    for table in ["phase9b3b_hwang_naive_models.tsv", "phase9b3b_hwang_treated_models.tsv"]:
        rows = read_tsv(TABLES / table)
        if not rows:
            fail(f"No model rows in {table}", errors)
        for row in rows:
            if row.get("replicate_unit") != "patient":
                fail(f"{table} has non-patient replicate unit", errors)
            if row.get("model_id", "").endswith("MODEL_A") and "patient_id:ROI_id" not in row.get("random_effect_structure", ""):
                fail(f"{table} Model A ignores Hwang ROI pairing", errors)
            if row.get("cohort_id") == "POOLED":
                fail(f"{table} pools Hwang cohorts", errors)

    moncada = read_tsv(TABLES / "phase9b3b_moncada_exploratory_results.tsv")
    for row in moncada:
        if "formal" in row.get("evidence_claim", "").lower() or "replication" in row.get("evidence_claim", "").lower():
            fail("Moncada reported as formal replication", errors)
        if row.get("replicate_unit") != "patient":
            fail("Moncada replicate unit is not patient", errors)

    coverage = read_tsv(TABLES / "phase9b3b_feature_coverage.tsv")
    for row in coverage:
        if row.get("feature_layer") == "WGCNA_module":
            cov = float(row.get("coverage_fraction") or "0")
            if cov < 0.80 and row.get("formal_inference_status") != "INSUFFICIENT_SPATIAL_DATA":
                fail("WGCNA module below 80% coverage entered formal inference", errors)

    neg = read_tsv(TABLES / "phase9b3b_negative_control_results.tsv")
    required_controls = {
        "within-section coordinate permutation",
        "size-matched random gene set",
        "expression-matched random gene set",
        "unrelated Hallmark pathway",
        "label permutation",
        "leakage check",
    }
    observed = {r.get("control_type") for r in neg if r.get("execution_status") == "EXECUTED"}
    missing = required_controls - observed
    if missing:
        fail(f"Missing executed negative controls: {sorted(missing)}", errors)
    if any(r.get("execution_status") in {"TO_VERIFY", "PLACEHOLDER", "MISSING"} for r in neg):
        fail("Negative controls include placeholder/TO_VERIFY/missing rows", errors)

    text = REPORT.read_text(encoding="utf-8").lower()
    banned = ["microbial localization", "microbiome validation", "causal mediation", "causality"]
    for term in banned:
        if term in text:
            fail(f"Report contains banned claim term: {term}", errors)

    if errors:
        print("Phase 9B3B validation FAILED")
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1
    print("Phase 9B3B validation PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
