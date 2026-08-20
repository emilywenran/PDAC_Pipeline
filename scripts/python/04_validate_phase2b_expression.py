#!/usr/bin/env python3
"""Validate Phase 2B analysis-ready expression outputs."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"VALIDATION FAILED: {message}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest_sha(file_id: str) -> str:
    manifest = pd.read_csv(ROOT / "01_metadata/file_manifest.tsv", sep="\t", dtype=str)
    row = manifest.loc[manifest["file_id"] == file_id]
    require(row.shape[0] == 1, f"file manifest missing unique row for {file_id}")
    value = row.iloc[0]["md5"]
    require(value.startswith("sha256:"), f"manifest checksum for {file_id} is not sha256")
    return value.split("sha256:", 1)[1]


def main() -> None:
    audited_path = ROOT / "03_processed/expression/GSE172356_expression_audited.tsv.gz"
    norm_path = ROOT / "03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz"
    log2_path = ROOT / "03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz"
    xwalk_path = ROOT / "01_metadata/expression_sample_crosswalk.tsv"
    missing_sample_path = ROOT / "05_results/tables/phase2b_missingness_by_sample.tsv"
    missing_gene_path = ROOT / "05_results/tables/phase2b_missingness_by_gene.tsv"
    filtering_path = ROOT / "05_results/tables/phase2b_filtering_sensitivity.tsv"
    outlier_path = ROOT / "05_results/tables/phase2b_outlier_assessment.tsv"
    report_path = ROOT / "04_analysis/03_expression_qc/PHASE2B_ANALYSIS_READY_EXPRESSION.md"
    figures = [
        ROOT / "05_results/figures/phase2b_missingness_heatmap.pdf",
        ROOT / "05_results/figures/phase2b_transformed_pca.pdf",
        ROOT / "05_results/figures/phase2b_transformed_sample_correlation.pdf",
        ROOT / "05_results/figures/phase2b_sample_qc_summary.pdf",
    ]
    scripts = [
        ROOT / "06_scripts/python/04_phase2b_prepare_expression.py",
        ROOT / "06_scripts/python/04_validate_phase2b_expression.py",
    ]
    for path in [
        audited_path,
        norm_path,
        log2_path,
        xwalk_path,
        missing_sample_path,
        missing_gene_path,
        filtering_path,
        outlier_path,
        report_path,
        *figures,
        *scripts,
    ]:
        require(path.exists(), f"missing output {path}")
        require(path.stat().st_size > 0, f"empty output {path}")

    require(sha256(audited_path) == manifest_sha("GSE172356_expression_audited"), "Phase 2A audited matrix checksum changed")

    xwalk = pd.read_csv(xwalk_path, sep="\t", dtype=str)
    expected_samples = xwalk["expression_column"].tolist()
    require(len(expected_samples) == 62, "crosswalk does not contain 62 samples")
    require(not pd.Index(expected_samples).duplicated().any(), "crosswalk sample order contains duplicates")

    audited = pd.read_csv(audited_path, sep="\t", compression="gzip")
    require(audited.shape == (45140, 63), f"unexpected audited matrix shape {audited.shape}")
    require(audited.columns[1:].tolist() == expected_samples, "audited matrix sample order differs from expression crosswalk")
    require(int(audited.drop(columns=["gene"]).isna().sum().sum()) == 73202, "audited matrix missing count changed")

    norm = pd.read_csv(norm_path, sep="\t", compression="gzip")
    log2 = pd.read_csv(log2_path, sep="\t", compression="gzip")
    require(norm.shape == log2.shape, "normalized and log2 matrix shapes differ")
    require(norm.columns.tolist() == log2.columns.tolist(), "normalized and log2 matrix columns differ")
    require(norm.columns[1:].tolist() == expected_samples, "final sample order differs from expression crosswalk")
    require(norm.shape[1] == 63, f"expected 62 samples plus gene column, observed {norm.shape[1]}")
    require(norm.shape[0] > 0, "no genes retained")
    require(norm["gene"].is_unique, "final normalized matrix contains duplicated genes")
    require(log2["gene"].is_unique, "final log2 matrix contains duplicated genes")
    require(not pd.Index(norm.columns[1:]).duplicated().any(), "final normalized matrix contains duplicated samples")
    require(not pd.Index(log2.columns[1:]).duplicated().any(), "final log2 matrix contains duplicated samples")

    norm_expr = norm.drop(columns=["gene"])
    log2_expr = log2.drop(columns=["gene"])
    require(int(norm_expr.isna().sum().sum()) == 0, "final normalized matrix contains missing values")
    require(int(log2_expr.isna().sum().sum()) == 0, "final log2 matrix contains missing values")
    require(np.isfinite(norm_expr.to_numpy(dtype=float)).all(), "final normalized matrix contains non-finite values")
    require(np.isfinite(log2_expr.to_numpy(dtype=float)).all(), "final log2 matrix contains non-finite values")
    require((norm_expr.to_numpy(dtype=float) >= 0).all(), "final normalized matrix contains negative values")
    require(np.allclose(np.log2(norm_expr.to_numpy(dtype=float) + 1), log2_expr.to_numpy(dtype=float), rtol=1e-10, atol=1e-10), "log2 transform is not reproducible")

    missing_sample = pd.read_csv(missing_sample_path, sep="\t")
    missing_gene = pd.read_csv(missing_gene_path, sep="\t")
    filtering = pd.read_csv(filtering_path, sep="\t")
    outlier = pd.read_csv(outlier_path, sep="\t")
    require(missing_sample.shape[0] == 62, "missingness-by-sample table does not contain 62 rows")
    require(missing_gene.shape[0] == 45140, "missingness-by-gene table does not contain 45,140 rows")
    require(int(missing_sample["missing_count"].sum()) == 73202, "sample missing counts do not sum to 73,202")
    require(int(missing_gene["missing_count"].sum()) == 73202, "gene missing counts do not sum to 73,202")
    require(int(missing_sample["literal_NA"].sum()) == 73202, "literal NA source count is not 73,202")
    require(int(missing_sample["blank_field"].sum()) == 0, "blank source count is not zero")
    require(int(missing_sample["parse_failure"].sum()) == 0, "parse-failure source count is not zero")
    require(filtering["selected_for_primary_matrix"].sum() == 1, "filtering sensitivity table does not mark exactly one primary rule")
    selected_genes = int(filtering.loc[filtering["selected_for_primary_matrix"], "genes_retained"].iloc[0])
    require(selected_genes == norm.shape[0], "selected filtering gene count does not match final matrix")
    require(outlier.shape[0] == 62, "outlier assessment does not contain 62 rows")
    require(set(outlier["phase2b_sample_classification"]).issubset({"RETAIN", "RETAIN_WITH_SENSITIVITY_ANALYSIS", "EXCLUDE_RECOMMENDED", "TO_VERIFY"}), "invalid sample classification")
    require(int((outlier["phase2b_sample_classification"] == "EXCLUDE_RECOMMENDED").sum()) == 0, "Phase 2B recommended exclusion without approval")

    report = report_path.read_text()
    for text in [
        "literal `NA`",
        "complete-observation",
        "log2(normalized_count + 1)",
        "Phase 3 subtype reproduction may proceed",
        "TO_VERIFY",
    ]:
        require(text in report, f"report missing required text: {text}")

    print("Phase 2B expression validation passed.")


if __name__ == "__main__":
    main()
