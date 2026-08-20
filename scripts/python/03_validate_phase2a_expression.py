#!/usr/bin/env python3
"""Validate Phase 2A expression audit outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"VALIDATION FAILED: {message}")


def main() -> None:
    matrix_path = ROOT / "03_processed/expression/GSE172356_expression_audited.tsv.gz"
    annot_path = ROOT / "03_processed/expression/GSE172356_gene_annotation.tsv"
    xwalk_path = ROOT / "01_metadata/expression_sample_crosswalk.tsv"
    sample_qc_path = ROOT / "05_results/tables/phase2a_expression_sample_qc.tsv"
    gene_qc_path = ROOT / "05_results/tables/phase2a_expression_gene_qc.tsv"
    mapping_path = ROOT / "05_results/tables/phase2a_expression_mapping_summary.tsv"
    report_path = ROOT / "04_analysis/03_expression_qc/PHASE2A_EXPRESSION_AUDIT.md"
    figures = [
        ROOT / "05_results/figures/phase2a_expression_distribution.pdf",
        ROOT / "05_results/figures/phase2a_sample_correlation.pdf",
        ROOT / "05_results/figures/phase2a_expression_pca.pdf",
    ]

    for path in [matrix_path, annot_path, xwalk_path, sample_qc_path, gene_qc_path, mapping_path, report_path, *figures]:
        require(path.exists(), f"missing output {path}")
        require(path.stat().st_size > 0, f"empty output {path}")

    matrix = pd.read_csv(matrix_path, sep="\t", compression="gzip")
    expr = matrix.drop(columns=["gene"])
    require(matrix.shape == (45140, 63), f"unexpected matrix shape {matrix.shape}")
    require(expr.shape[1] == 62, f"unexpected sample count {expr.shape[1]}")
    require(matrix["gene"].is_unique, "gene identifiers are not unique")
    require(not pd.Index(expr.columns).duplicated().any(), "sample identifiers are duplicated")
    require(int(expr.isna().sum().sum()) == 73202, "unexpected missing-value count")
    require(int((~np.isfinite(expr.to_numpy(dtype=float)) & ~expr.isna().to_numpy()).sum()) == 0, "infinite values detected")
    require(int((expr < 0).sum().sum()) == 0, "negative values detected")

    xwalk = pd.read_csv(xwalk_path, sep="\t")
    require(xwalk.shape[0] == 62, f"unexpected crosswalk row count {xwalk.shape[0]}")
    require(set(xwalk["mapping_status"]) == {"MAPPED"}, "not all expression samples are mapped")
    require(xwalk["geo_sample_id"].nunique() == 62, "GEO IDs are not unique in expression crosswalk")
    require(xwalk["patient_id"].nunique() == 62, "patient IDs are not unique in expression crosswalk")

    sample_qc = pd.read_csv(sample_qc_path, sep="\t")
    gene_qc = pd.read_csv(gene_qc_path, sep="\t")
    require(sample_qc.shape[0] == 62, "unexpected sample QC row count")
    require(gene_qc.shape[0] == 45140, "unexpected gene QC row count")
    require(int(sample_qc["missing_count"].sum()) == 73202, "sample QC missing counts do not sum to matrix missingness")
    require(int(gene_qc["missing_count"].sum()) == 73202, "gene QC missing counts do not sum to matrix missingness")
    require(int(sample_qc["infinite_count"].sum()) == 0, "sample QC infinite counts are nonzero")
    require(int(gene_qc["infinite_count"].sum()) == 0, "gene QC infinite counts are nonzero")

    mapping = pd.read_csv(mapping_path, sep="\t")
    mapped = int(mapping.loc[mapping["mapping_status"] == "MAPPED", "sample_count"].sum())
    require(mapped == 62, f"unexpected mapped sample count {mapped}")

    report = report_path.read_text()
    for text in [
        "DESeq size-factor-normalized counts",
        "Matrix dimensions: 45140 gene rows x 62 expression samples",
        "Successfully mapped samples: 62 / 62",
        "TO_VERIFY",
    ]:
        require(text in report, f"report missing required text: {text}")

    print("Phase 2A expression validation passed.")


if __name__ == "__main__":
    main()
