#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache-phase7a5")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def p(*parts: str) -> Path:
    return ROOT.joinpath(*parts)


def fail(message: str) -> None:
    raise SystemExit(f"VALIDATION_FAIL: {message}")


def require_file(path: Path) -> None:
    if not path.exists():
        fail(f"Missing required file: {path}")
    if path.stat().st_size == 0:
        fail(f"Required file is empty: {path}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    required = [
        p("01_metadata", "host_tme_covariates.tsv"),
        p("05_results", "tables", "phase7a5_host_covariate_qc.tsv"),
        p("05_results", "tables", "phase7a5_host_covariate_correlations.tsv"),
        p("05_results", "tables", "phase7a5_covariate_model_feasibility.tsv"),
        p("05_results", "figures", "phase7a5_host_covariate_distributions.pdf"),
        p("05_results", "figures", "phase7a5_host_covariate_correlation.pdf"),
    ]
    for path in required:
        require_file(path)

    cov = pd.read_csv(p("01_metadata", "host_tme_covariates.tsv"), sep="\t")
    crosswalk = pd.read_csv(p("01_metadata", "expression_sample_crosswalk.tsv"), sep="\t")
    manifest = pd.read_csv(p("01_metadata", "sample_manifest.tsv"), sep="\t")
    expr_header = pd.read_csv(
        p("03_processed", "expression", "GSE172356_expression_log2_analysis_ready.tsv.gz"),
        sep="\t",
        nrows=0,
    )

    expected_cols = [
        "patient_id",
        "expression_sample_id",
        "stromal_score",
        "immune_score",
        "estimate_score",
        "inferred_tumor_purity",
        "method",
        "method_version",
        "gene_coverage",
        "validation_status",
        "notes",
    ]
    missing = [col for col in expected_cols if col not in cov.columns]
    if missing:
        fail(f"host_tme_covariates.tsv missing columns: {missing}")
    if len(cov) != 62:
        fail(f"Expected 62 covariate rows, found {len(cov)}")
    if cov["patient_id"].nunique() != 62:
        fail("Expected 62 unique patient IDs")
    if cov["patient_id"].duplicated().any():
        fail("Duplicated patient IDs detected")
    if cov["expression_sample_id"].duplicated().any():
        fail("Duplicated expression sample IDs detected")
    if cov["patient_id"].tolist() != crosswalk["patient_id"].tolist():
        fail("Covariate patient order does not match expression crosswalk")
    if cov["expression_sample_id"].tolist() != crosswalk["expression_column"].tolist():
        fail("Covariate expression sample order does not match expression crosswalk")
    if list(expr_header.columns[1:]) != crosswalk["expression_column"].tolist():
        fail("Expression matrix sample order does not match project crosswalk")
    if set(cov["patient_id"]) != set(manifest["patient_id"]):
        fail("Covariate patients do not match sample manifest patients")

    score_cols = ["stromal_score", "immune_score", "estimate_score", "inferred_tumor_purity"]
    score_values = cov[score_cols].to_numpy(dtype=float)
    if np.isnan(score_values).any():
        fail("Missing ESTIMATE-derived scores detected")
    if not np.isfinite(score_values).all():
        fail("Non-finite ESTIMATE-derived scores detected")
    purity = cov["inferred_tumor_purity"].astype(float)
    if not ((purity >= 0).all() and (purity <= 1).all()):
        fail("Inferred tumor purity values must be bounded in [0, 1]")
    if (cov["validation_status"] != "PASS").any():
        fail("All host_tme_covariates validation_status values must be PASS")
    notes = " ".join(cov["notes"].astype(str).unique())
    if "subtype_and_microbiome_not_used_for_score_generation" not in notes:
        fail("Covariate notes must state subtype/microbiome information was not used for score generation")
    if not cov["method_version"].str.contains("estimate_1.0.13", regex=False).all():
        fail("method_version must record official estimate package version 1.0.13")

    qc = pd.read_csv(p("05_results", "tables", "phase7a5_host_covariate_qc.tsv"), sep="\t")
    if set(qc["metric"]) != set(score_cols):
        fail("QC table must summarize exactly the required score columns")
    if (qc["missing_count"].astype(int) != 0).any() or (qc["infinite_count"].astype(int) != 0).any():
        fail("QC table reports missing or infinite scores")

    corr = pd.read_csv(p("05_results", "tables", "phase7a5_host_covariate_correlations.tsv"), sep="\t")
    if len(corr) != 6:
        fail("Expected 6 pairwise host covariate correlations")
    if not corr["spearman_rho"].between(-1, 1).all():
        fail("Spearman rho outside [-1, 1]")

    feasibility = pd.read_csv(p("05_results", "tables", "phase7a5_covariate_model_feasibility.tsv"), sep="\t")
    expected_feas_cols = [
        "model_id",
        "covariates",
        "available_patients",
        "correlation_warning",
        "maximum_VIF",
        "condition_number",
        "model_permitted",
        "reason",
        "analysis_role",
        "notes",
    ]
    missing = [col for col in expected_feas_cols if col not in feasibility.columns]
    if missing:
        fail(f"Feasibility table missing columns: {missing}")
    for model_id in ["Model_3P", "Model_3I", "Model_3S"]:
        if model_id not in set(feasibility["model_id"]):
            fail(f"Missing prespecified sensitivity model {model_id}")
    if "Model_3ALL_TME" not in set(feasibility["model_id"]):
        fail("Missing combined TME screen row documenting whether all ESTIMATE covariates can be combined")
    if (feasibility["available_patients"].astype(int) < 62).any():
        fail("All Phase 7A.5 model feasibility rows should have 62 complete patients")
    if feasibility.loc[feasibility["model_id"] == "Model_0", "model_permitted"].iloc[0] != "YES":
        fail("Primary Model 0 feasibility must remain permitted")

    forbidden_outputs = [
        p("05_results", "tables", "phase7b_microbiome_associations.tsv"),
        p("05_results", "tables", "phase7a_genus_association_results.tsv"),
    ]
    for path in forbidden_outputs:
        if path.exists():
            fail(f"Forbidden microbiome association output exists during Phase 7A.5: {path}")

    print("VALIDATION_PASS: Phase 7A.5 host covariates validated")
    for path in required:
        print(f"{path.relative_to(ROOT)}\tsha256:{sha256(path)}")


if __name__ == "__main__":
    main()
