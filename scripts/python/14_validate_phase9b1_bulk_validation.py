#!/usr/bin/env python3
from pathlib import Path
import gzip
import sys

import pandas as pd

ROOT = Path("/Users/emily/thesis/PDAC")
TABLE = ROOT / "05_results/tables"
FIG = ROOT / "05_results/figures"
REPORT = ROOT / "04_analysis/09_external_validation/PHASE9B1_BULK_EXTERNAL_VALIDATION_RESULTS.md"

REQUIRED_TABLES = [
    "phase9b1_bulk_cohort_qc.tsv",
    "phase9b1_signature_coverage.tsv",
    "phase9b1_bulk_state_scores.tsv.gz",
    "phase9b1_bulk_host_feature_scores.tsv.gz",
    "phase9b1_module_transfer_coverage.tsv",
    "phase9b1_cohort_replication_results.tsv",
    "phase9b1_negative_control_results.tsv",
    "phase9b1_cross_cohort_synthesis.tsv",
    "phase9b1_host_feature_replication_evidence.tsv",
]

REQUIRED_FIGURES = [
    "phase9b1_bulk_cohort_qc.pdf",
    "phase9b1_axis_score_distributions.pdf",
    "phase9b1_pathway_replication_forest.pdf",
    "phase9b1_tf_replication_heatmap.pdf",
    "phase9b1_module_replication_forest.pdf",
    "phase9b1_cross_cohort_summary.pdf",
    "phase9b1_negative_control_summary.pdf",
]


def require(path):
    if not path.exists() or path.stat().st_size == 0:
        raise AssertionError(f"Missing or empty required file: {path}")


def read_table(name):
    path = TABLE / name
    require(path)
    return pd.read_csv(path, sep="\t")


def main():
    for name in REQUIRED_TABLES:
        require(TABLE / name)
    for name in REQUIRED_FIGURES:
        require(FIG / name)
    require(REPORT)

    qc = read_table("phase9b1_bulk_cohort_qc.tsv")
    expected = {"TCGA_PAAD", "GSE71729", "GSE62452"}
    observed = set(qc["dataset_id"])
    if observed != expected:
        raise AssertionError(f"Unexpected analyzed cohorts: {observed}; expected {expected}")
    if (qc["analyzed_samples"] < 30).any():
        raise AssertionError("At least one bulk cohort has fewer than 30 analyzed samples.")

    coverage = read_table("phase9b1_signature_coverage.tsv")
    key = coverage[coverage["signature_name"].isin(["Moffitt50_basal", "Moffitt50_classical"])]
    if key.empty or (key["coverage_fraction"] < 0.8).any():
        raise AssertionError("Moffitt50 signature coverage failed the locked 80% threshold.")

    repl = read_table("phase9b1_cohort_replication_results.tsv")
    families = set(repl["feature_family"])
    if not {"pathway", "tf", "module"}.issubset(families):
        raise AssertionError(f"Missing replication feature families: {families}")
    if {"feature_score ~ moffitt50_contrast"} != set(repl["model"]):
        raise AssertionError("Unexpected replication model specification.")

    evidence = read_table("phase9b1_host_feature_replication_evidence.tsv")
    allowed = {
        "EXTERNALLY_REPLICATED_HOST_FEATURE",
        "PARTIALLY_REPLICATED_HOST_FEATURE",
        "NOT_REPLICATED",
        "INSUFFICIENT_EXTERNAL_DATA",
        "TO_VERIFY",
    }
    if not set(evidence["phase9a_evidence_category"]).issubset(allowed):
        raise AssertionError("Evidence table contains non-Phase-9A category labels.")

    text = REPORT.read_text()
    banned = ["single-cell validation executed", "spatial validation executed", "microbiome validation executed"]
    if any(b in text.lower() for b in banned):
        raise AssertionError("Report appears to claim a prohibited validation layer was executed.")

    print("Phase 9B1 validation checks passed.")


if __name__ == "__main__":
    main()
