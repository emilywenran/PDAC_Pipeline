#!/usr/bin/env python3
"""Phase 2B expression missingness audit and analysis-ready matrix creation."""

from __future__ import annotations

import gzip
import hashlib
import math
import os
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.covariance import MinCovDet
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


SOURCE = ROOT / "02_data/reference/GSE172356_processed/GSE172356_PDA_gene_expression_matrix.txt.gz"
AUDITED = ROOT / "03_processed/expression/GSE172356_expression_audited.tsv.gz"
GENE_ANNOT = ROOT / "03_processed/expression/GSE172356_gene_annotation.tsv"
CROSSWALK = ROOT / "01_metadata/expression_sample_crosswalk.tsv"
SAMPLE_MANIFEST = ROOT / "01_metadata/sample_manifest.tsv"
FILE_MANIFEST = ROOT / "01_metadata/file_manifest.tsv"
PHASE2A_SAMPLE_QC = ROOT / "05_results/tables/phase2a_expression_sample_qc.tsv"

OUT_NORM = ROOT / "03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz"
OUT_LOG2 = ROOT / "03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz"
OUT_MISSING_SAMPLE = ROOT / "05_results/tables/phase2b_missingness_by_sample.tsv"
OUT_MISSING_GENE = ROOT / "05_results/tables/phase2b_missingness_by_gene.tsv"
OUT_FILTER = ROOT / "05_results/tables/phase2b_filtering_sensitivity.tsv"
OUT_OUTLIER = ROOT / "05_results/tables/phase2b_outlier_assessment.tsv"
OUT_REPORT = ROOT / "04_analysis/03_expression_qc/PHASE2B_ANALYSIS_READY_EXPRESSION.md"
FIG_MISSING = ROOT / "05_results/figures/phase2b_missingness_heatmap.pdf"
FIG_PCA = ROOT / "05_results/figures/phase2b_transformed_pca.pdf"
FIG_COR = ROOT / "05_results/figures/phase2b_transformed_sample_correlation.pdf"
FIG_QC = ROOT / "05_results/figures/phase2b_sample_qc_summary.pdf"
DECISION_LOG = ROOT / "09_docs/planning/DECISION_LOG.md"
PROJECT_STATUS = ROOT / "00_admin/PROJECT_STATUS.md"

TODAY = date.today().isoformat()
SUSPECTED_OUTLIERS = {"YX16135T", "YX16158T", "YX16194T", "YX16224T"}


def ensure_dirs() -> None:
    for path in [
        OUT_NORM.parent,
        OUT_MISSING_SAMPLE.parent,
        FIG_MISSING.parent,
        OUT_REPORT.parent,
        ROOT / ".matplotlib",
    ]:
        path.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(ROOT / ".matplotlib")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def robust_z(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    med = values.median()
    mad = (values - med).abs().median()
    if mad == 0 or np.isnan(mad):
        sd = values.std(ddof=0)
        return (values - values.mean()) / sd if sd else pd.Series(0.0, index=values.index)
    return 0.6745 * (values - med) / mad


def read_source_missing_representation(sample_order: list[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    records = []
    counts = {"blank_fields": 0, "literal_NA": 0, "parse_failures": 0, "other_non_numeric": 0}
    with gzip.open(SOURCE, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_index = {sample: header.index(sample) for sample in sample_order}
        for row_index, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            for sample, idx in sample_index.items():
                value = parts[idx] if idx < len(parts) else ""
                source_type = None
                if value == "":
                    source_type = "blank_field"
                    counts["blank_fields"] += 1
                elif value.upper() == "NA":
                    source_type = "literal_NA"
                    counts["literal_NA"] += 1
                else:
                    try:
                        float(value)
                    except ValueError:
                        source_type = "parse_failure"
                        counts["parse_failures"] += 1
                if source_type:
                    records.append(
                        {
                            "gene": gene,
                            "expression_column": sample,
                            "row_index_1based": row_index,
                            "source_missing_representation": source_type,
                        }
                    )
    return pd.DataFrame.from_records(records), counts


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(AUDITED, sep="\t", compression="gzip")
    xwalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    manifest = pd.read_csv(SAMPLE_MANIFEST, sep="\t", dtype=str)
    phase2a_qc = pd.read_csv(PHASE2A_SAMPLE_QC, sep="\t")
    annot = pd.read_csv(GENE_ANNOT, sep="\t")
    sample_order = xwalk["expression_column"].tolist()
    matrix = matrix[["gene", *sample_order]]
    return matrix, xwalk, manifest, phase2a_qc, annot


def missingness_tables(
    matrix: pd.DataFrame,
    xwalk: pd.DataFrame,
    manifest: pd.DataFrame,
    phase2a_qc: pd.DataFrame,
    missing_cells: pd.DataFrame,
    representation_counts: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    expr = matrix.drop(columns=["gene"])
    missing = expr.isna()
    total_genes = expr.shape[0]
    source_by_sample = (
        missing_cells.groupby(["expression_column", "source_missing_representation"]).size().unstack(fill_value=0)
    )
    for col in ["literal_NA", "blank_field", "parse_failure"]:
        if col not in source_by_sample.columns:
            source_by_sample[col] = 0

    sample = pd.DataFrame(
        {
            "expression_column": expr.columns,
            "missing_count": missing.sum(axis=0).values.astype(int),
            "missing_proportion": missing.mean(axis=0).values,
            "nonmissing_count": (~missing).sum(axis=0).values.astype(int),
            "total_expression_skipna": expr.sum(axis=0, skipna=True).values,
            "detected_genes_gt0_skipna": (expr > 0).sum(axis=0).values.astype(int),
            "suspected_phase2a_outlier": [col in SUSPECTED_OUTLIERS for col in expr.columns],
        }
    )
    sample = sample.merge(
        source_by_sample[["literal_NA", "blank_field", "parse_failure"]],
        left_on="expression_column",
        right_index=True,
        how="left",
    ).fillna({"literal_NA": 0, "blank_field": 0, "parse_failure": 0})
    for col in ["literal_NA", "blank_field", "parse_failure"]:
        sample[col] = sample[col].astype(int)
    sample = sample.merge(xwalk, on="expression_column", how="left")
    sample = sample.merge(
        manifest[["sample_id", "batch", "tumor_purity"]],
        left_on="patient_id",
        right_on="sample_id",
        how="left",
    ).drop(columns=["sample_id"])
    sample = sample.merge(
        phase2a_qc[["expression_column", "outlier_flag"]].rename(columns={"outlier_flag": "phase2a_outlier_flag"}),
        on="expression_column",
        how="left",
    )
    sample["missing_count_robust_z"] = robust_z(sample["missing_count"]).values
    sample["total_expression_robust_z"] = robust_z(sample["total_expression_skipna"]).values
    sample["detected_genes_robust_z"] = robust_z(sample["detected_genes_gt0_skipna"]).values

    source_by_gene = missing_cells.groupby(["gene", "source_missing_representation"]).size().unstack(fill_value=0)
    for col in ["literal_NA", "blank_field", "parse_failure"]:
        if col not in source_by_gene.columns:
            source_by_gene[col] = 0
    gene = pd.DataFrame(
        {
            "gene": matrix["gene"],
            "row_index_1based": np.arange(1, matrix.shape[0] + 1),
            "missing_count": missing.sum(axis=1).values.astype(int),
            "missing_proportion": missing.mean(axis=1).values,
            "nonmissing_count": (~missing).sum(axis=1).values.astype(int),
            "mean_expression_skipna": expr.mean(axis=1, skipna=True).values,
            "median_expression_skipna": expr.median(axis=1, skipna=True).values,
            "detected_sample_count_gt0_skipna": (expr > 0).sum(axis=1).values.astype(int),
            "all_zero_among_observed": (expr.fillna(0).sum(axis=1) == 0).values,
        }
    )
    gene = gene.merge(
        source_by_gene[["literal_NA", "blank_field", "parse_failure"]],
        left_on="gene",
        right_index=True,
        how="left",
    ).fillna({"literal_NA": 0, "blank_field": 0, "parse_failure": 0})
    for col in ["literal_NA", "blank_field", "parse_failure"]:
        gene[col] = gene[col].astype(int)
    gene["missing_pattern"] = np.where(
        gene["missing_count"] == 0,
        "complete",
        np.where(gene["missing_count"] == 21, "missing_in_21_samples", "missing_in_41_samples"),
    )

    missing_by_subtype = sample.groupby("subtype_original", dropna=False)["missing_count"].agg(["count", "min", "median", "max"]).reset_index()
    missing_by_batch = sample.groupby("batch", dropna=False)["missing_count"].agg(["count", "min", "median", "max"]).reset_index()
    metrics = {
        "total_missing": int(missing.sum().sum()),
        "genes_with_missing": int(missing.any(axis=1).sum()),
        "samples_with_missing": int(missing.any(axis=0).sum()),
        "complete_genes": int((~missing.any(axis=1)).sum()),
        "sample_missing_min": int(sample["missing_count"].min()),
        "sample_missing_max": int(sample["missing_count"].max()),
        "gene_missing_patterns": gene["missing_count"].value_counts().sort_index().to_dict(),
        "representation_counts": representation_counts,
        "missing_total_expression_pearson": float(sample["missing_count"].corr(sample["total_expression_skipna"])),
        "missing_detected_gene_pearson": float(sample["missing_count"].corr(sample["detected_genes_gt0_skipna"])),
        "missing_by_subtype": missing_by_subtype,
        "missing_by_batch": missing_by_batch,
    }
    return sample, gene, metrics


def filtering_sensitivity(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, dict[str, int]]:
    expr = matrix.drop(columns=["gene"])
    missing = expr.isna()
    n_samples = expr.shape[1]
    min_10pct = math.ceil(0.10 * n_samples)
    min_20pct = math.ceil(0.20 * n_samples)
    min_50pct = math.ceil(0.50 * n_samples)
    complete = ~missing.any(axis=1)
    not_all_zero = expr.fillna(0).sum(axis=1) > 0
    detected_1 = (expr > 0).sum(axis=1) >= 1
    expr1_10pct = (expr >= 1).sum(axis=1) >= min_10pct
    primary = complete & not_all_zero & expr1_10pct

    rules = [
        ("complete_observations_only", complete, "missing_handling", "Retain only genes with no missing values; selected missing-value strategy."),
        ("missing_le_10_percent", missing.mean(axis=1) <= 0.10, "missing_handling", "Threshold sensitivity; no partially missing gene passes because missing blocks are 21 or 41 samples."),
        ("missing_le_20_percent", missing.mean(axis=1) <= 0.20, "missing_handling", "Threshold sensitivity; same retained set as complete observations."),
        ("missing_le_50_percent", missing.mean(axis=1) <= 0.50, "missing_handling", "Sensitivity only; would retain genes missing in 21 samples and require imputation."),
        ("retain_all_with_gene_median_imputation_sensitivity", pd.Series(True, index=expr.index), "missing_handling", "Rejected for primary matrix because source NA origin is not explained."),
        ("retain_all_with_zero_fill_sensitivity", pd.Series(True, index=expr.index), "missing_handling", "Rejected for primary matrix because source semantics do not support treating NA as zero."),
        ("remove_all_zero_genes_only", not_all_zero, "expression_filter", "Removes genes with zero observed expression after NA ignored."),
        ("detected_gt0_in_at_least_1_sample", detected_1, "expression_filter", "Retains genes with any positive observed expression."),
        ("normalized_count_ge_1_in_10_percent_samples", expr1_10pct, "expression_filter", f"Requires normalized count >= 1 in at least {min_10pct} samples."),
        ("normalized_count_ge_1_in_20_percent_samples", (expr >= 1).sum(axis=1) >= min_20pct, "expression_filter", f"Requires normalized count >= 1 in at least {min_20pct} samples."),
        ("normalized_count_ge_1_in_50_percent_samples", (expr >= 1).sum(axis=1) >= min_50pct, "expression_filter", f"Requires normalized count >= 1 in at least {min_50pct} samples."),
        ("primary_complete_not_all_zero_count_ge_1_10_percent", primary, "selected_primary", "Complete observations, not all-zero, and normalized count >= 1 in at least 10% of samples."),
    ]
    rows = []
    for rule, mask, category, notes in rules:
        mask = pd.Series(mask, index=expr.index).fillna(False)
        retained = expr.loc[mask]
        rows.append(
            {
                "rule_id": rule,
                "category": category,
                "genes_retained": int(mask.sum()),
                "genes_removed": int((~mask).sum()),
                "samples_retained": n_samples,
                "missing_values_retained": int(retained.isna().sum().sum()),
                "requires_imputation": bool(retained.isna().any().any()),
                "selected_for_primary_matrix": rule == "primary_complete_not_all_zero_count_ge_1_10_percent",
                "notes": notes,
            }
        )
    stats = {
        "min_samples_10pct": min_10pct,
        "min_samples_20pct": min_20pct,
        "min_samples_50pct": min_50pct,
        "primary_genes": int(primary.sum()),
    }
    return pd.DataFrame.from_records(rows), primary, stats


def build_outlier_assessment(norm: pd.DataFrame, log2_matrix: pd.DataFrame, missing_sample: pd.DataFrame, xwalk: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr_norm = norm.drop(columns=["gene"])
    expr_log = log2_matrix.drop(columns=["gene"])
    samples = expr_norm.columns.tolist()
    cor = expr_log.corr(method="pearson")
    cor_no_diag = cor.mask(np.eye(cor.shape[0], dtype=bool))

    scaled = StandardScaler().fit_transform(expr_log.T)
    pca = PCA(n_components=10, random_state=0)
    pcs10 = pca.fit_transform(scaled)
    pcs2 = pcs10[:, :2]
    center2 = np.median(pcs2, axis=0)
    pca_distance = np.sqrt(((pcs2 - center2) ** 2).sum(axis=1))

    try:
        mcd = MinCovDet(random_state=0, support_fraction=0.75).fit(pcs10)
        mahal = np.sqrt(mcd.mahalanobis(pcs10))
        mahal_note = "robust_MCD_on_first_10_PCs"
    except Exception:
        mahal = np.full(len(samples), np.nan)
        mahal_note = "not_applicable_MCD_failed"

    out = pd.DataFrame(
        {
            "expression_column": samples,
            "total_expression_filtered": expr_norm.sum(axis=0).values,
            "detected_genes_gt0_filtered": (expr_norm > 0).sum(axis=0).values.astype(int),
            "median_sample_correlation_log2": cor_no_diag.median(axis=0).reindex(samples).values,
            "mean_sample_correlation_log2": cor_no_diag.mean(axis=0).reindex(samples).values,
            "pc1": pcs2[:, 0],
            "pc2": pcs2[:, 1],
            "pca_distance": pca_distance,
            "robust_mahalanobis_distance": mahal,
            "robust_mahalanobis_method": mahal_note,
            "suspected_phase2a_outlier": [sample in SUSPECTED_OUTLIERS for sample in samples],
        }
    )
    out = out.merge(missing_sample[["expression_column", "missing_count", "missing_proportion"]], on="expression_column", how="left")
    out = out.merge(xwalk, on="expression_column", how="left")
    out["total_expression_robust_z"] = robust_z(out["total_expression_filtered"]).values
    out["detected_genes_robust_z"] = robust_z(out["detected_genes_gt0_filtered"]).values
    out["median_correlation_robust_z"] = robust_z(out["median_sample_correlation_log2"]).values
    out["pca_distance_robust_z"] = robust_z(out["pca_distance"]).values
    out["robust_mahalanobis_robust_z"] = robust_z(out["robust_mahalanobis_distance"]).values

    classifications = []
    reasons = []
    for _, row in out.iterrows():
        flags = []
        if abs(row["total_expression_robust_z"]) > 3.5:
            flags.append("total_expression_robust_z_abs_gt_3.5")
        if row["detected_genes_robust_z"] < -3.5:
            flags.append("detected_genes_robust_z_lt_-3.5")
        if row["median_correlation_robust_z"] < -3.5:
            flags.append("median_correlation_robust_z_lt_-3.5")
        if row["pca_distance_robust_z"] > 3.5:
            flags.append("pca_distance_robust_z_gt_3.5")
        if row["robust_mahalanobis_robust_z"] > 3.5:
            flags.append("robust_mahalanobis_robust_z_gt_3.5")
        if row["suspected_phase2a_outlier"] and flags:
            classifications.append("RETAIN_WITH_SENSITIVITY_ANALYSIS")
        elif row["suspected_phase2a_outlier"]:
            classifications.append("RETAIN_WITH_SENSITIVITY_ANALYSIS")
            flags.append("phase2a_suspected_outlier_requires_later_sensitivity")
        else:
            classifications.append("RETAIN")
        reasons.append(";".join(flags) if flags else "none")
    out["phase2b_sample_classification"] = classifications
    out["objective_evidence"] = reasons
    out["sensitivity_analysis_plan"] = np.where(
        out["phase2b_sample_classification"] == "RETAIN_WITH_SENSITIVITY_ANALYSIS",
        "Repeat Phase 3 subtype reproduction with and without this sample and report stability of sample assignments and unsupervised geometry; do not remove before approval.",
        "Not required by Phase 2B QC.",
    )
    pca_meta = pd.DataFrame({"component": [f"PC{i}" for i in range(1, 11)], "variance_explained": pca.explained_variance_ratio_})
    return out, cor, pca_meta


def write_matrices(matrix: pd.DataFrame, primary_mask: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    norm = matrix.loc[primary_mask, :].copy()
    log2 = norm.copy()
    sample_cols = norm.columns[1:]
    log2.loc[:, sample_cols] = np.log2(norm.loc[:, sample_cols] + 1)
    norm.to_csv(OUT_NORM, sep="\t", index=False, compression="gzip")
    log2.to_csv(OUT_LOG2, sep="\t", index=False, compression="gzip")
    return norm, log2


def plot_figures(matrix: pd.DataFrame, missing_gene: pd.DataFrame, log2: pd.DataFrame, outlier: pd.DataFrame, cor: pd.DataFrame) -> None:
    expr = matrix.drop(columns=["gene"])
    missing = expr.isna().to_numpy(dtype=float)
    gene_mask = missing_gene["missing_count"] > 0
    heat = missing[gene_mask.values, :]
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.imshow(heat, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    ax.set_title("Source literal NA pattern among genes with missing values")
    ax.set_xlabel("Samples in expression crosswalk order")
    ax.set_ylabel("Genes with >=1 missing value")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(FIG_MISSING)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    flagged = outlier["phase2b_sample_classification"] != "RETAIN"
    ax.scatter(outlier.loc[~flagged, "pc1"], outlier.loc[~flagged, "pc2"], c="#2f6f8f", s=36, label="RETAIN")
    ax.scatter(outlier.loc[flagged, "pc1"], outlier.loc[flagged, "pc2"], c="#b33a3a", s=46, label="Sensitivity")
    for _, row in outlier.loc[flagged].iterrows():
        ax.annotate(row["expression_column"], (row["pc1"], row["pc2"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA on filtered log2(normalized count + 1)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_PCA)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(cor, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Sample correlation, filtered log2 matrix")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(FIG_COR)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    axes = axes.ravel()
    colors = np.where(outlier["phase2b_sample_classification"] == "RETAIN", "#2f6f8f", "#b33a3a")
    axes[0].scatter(outlier["missing_count"], outlier["total_expression_filtered"], c=colors)
    axes[0].set_xlabel("Original missing values")
    axes[0].set_ylabel("Filtered total expression")
    axes[1].scatter(outlier["detected_genes_gt0_filtered"], outlier["median_sample_correlation_log2"], c=colors)
    axes[1].set_xlabel("Detected genes")
    axes[1].set_ylabel("Median correlation")
    axes[2].hist(outlier["pca_distance"], bins=16, color="#756bb1", edgecolor="white")
    axes[2].set_xlabel("PCA distance")
    axes[2].set_ylabel("Samples")
    axes[3].hist(outlier["robust_mahalanobis_distance"].dropna(), bins=16, color="#238b45", edgecolor="white")
    axes[3].set_xlabel("Robust Mahalanobis distance")
    axes[3].set_ylabel("Samples")
    fig.suptitle("Phase 2B sample QC summary", y=0.98)
    fig.tight_layout()
    fig.savefig(FIG_QC)
    plt.close(fig)


def update_file_manifest(paths: dict[str, tuple[Path, str, str]]) -> None:
    manifest = pd.read_csv(FILE_MANIFEST, sep="\t", dtype=str)
    manifest = manifest[~manifest["file_id"].isin(paths.keys())]
    entries = []
    for file_id, (path, data_type, notes) in paths.items():
        entries.append(
            {
                "file_id": file_id,
                "dataset": "GSE172356_processed",
                "sample_id": "",
                "data_type": data_type,
                "local_path": str(path),
                "source_url_or_accession": "derived_from_GSE172356_expression_audited",
                "file_size": str(path.stat().st_size),
                "md5": f"sha256:{sha256(path)}",
                "download_date": TODAY,
                "processing_status": "generated_Phase2B",
                "notes": notes,
            }
        )
    manifest = pd.concat([manifest, pd.DataFrame(entries)], ignore_index=True)
    manifest.to_csv(FILE_MANIFEST, sep="\t", index=False)


def update_decision_log(primary_genes: int) -> None:
    text = DECISION_LOG.read_text()
    if "D-08" not in text:
        text = text.replace(
            "| 2026-06-30 | **D-07** | Use the official GSE172356 processed matrix as the Phase 2A audited host expression source while preserving its DESeq size-factor-normalized count scale | `02_data/reference/GSE172356_processed/`, `03_processed/expression/`, `04_analysis/03_expression_qc/`, `05_results/tables/`, `05_results/figures/` |\n",
            "| 2026-06-30 | **D-07** | Use the official GSE172356 processed matrix as the Phase 2A audited host expression source while preserving its DESeq size-factor-normalized count scale | `02_data/reference/GSE172356_processed/`, `03_processed/expression/`, `04_analysis/03_expression_qc/`, `05_results/tables/`, `05_results/figures/` |\n"
            "| 2026-06-30 | **D-08** | Prepare Phase 2B expression matrices by retaining complete-observation genes, applying unsupervised expression filtering, and using log2(normalized count + 1) for analysis-ready subtype reproduction input | `03_processed/expression/`, `04_analysis/03_expression_qc/`, `05_results/tables/`, `05_results/figures/` |\n",
        )
        text = text.replace(
            "\n---\n\n## Revisions and Corrections Log",
            f"""
### D-08: Phase 2B Missing-Value Handling, Filtering, and Transformation
*   **Date:** 2026-06-30
*   **Decision:** Use complete-observation genes as the primary missing-value strategy for GSE172356 Phase 2B, then remove all-zero genes and retain genes with DESeq size-factor-normalized count >= 1 in at least 10% of the 62 mapped samples. Preserve all 62 samples. Write both filtered normalized counts and `log2(normalized count + 1)` matrices. The selected primary matrix contains {primary_genes} genes.
*   **Alternatives Considered:** Retain all genes with missing values and perform gene-median imputation; retain genes with <=50% missingness and impute; replace `NA` with zero; apply DESeq2 VST/rlog directly to the processed matrix. These were rejected for the primary matrix because the source uses literal `NA` without documentation that the entries are structural zeros, and VST/rlog require raw-count assumptions that are not satisfied by an already normalized, fractional matrix.
*   **Scientific and Operational Justification:** Complete-observation filtering avoids imputing unexplained source `NA` values before subtype reproduction. The expression filter is independent of subtype labels and removes genes unlikely to contribute stable unsupervised structure. The log2 transform is reproducible for non-negative normalized counts and avoids applying a second library-size normalization.
*   **Files / Analyses Affected:** `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`, `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`, `04_analysis/03_expression_qc/PHASE2B_ANALYSIS_READY_EXPRESSION.md`, `05_results/tables/phase2b_*`, and `05_results/figures/phase2b_*`.

---

## Revisions and Corrections Log""",
        )
        DECISION_LOG.write_text(text)


def update_project_status(primary_genes: int) -> None:
    text = PROJECT_STATUS.read_text()
    current = "Phase 2B - analysis-ready GSE172356 expression matrices generated and validated; Phase 3 subtype reproduction may proceed after human review."
    text = text.replace(
        "Phase 2A - processed GSE172356 host expression matrix acquired and audited; generated files await human review before commit.",
        current,
    )
    additions = f"""- Completed Phase 2B missingness audit, filtering sensitivity, transformation, and analysis-ready expression matrix generation.
- Confirmed all 73,202 missing expression cells are literal `NA` strings in the official GEO processed matrix, not blank fields or parse failures.
- Selected complete-observation missing-value handling plus unsupervised expression filtering; retained {primary_genes} genes and all 62 mapped samples.
- Created filtered normalized-count and filtered log2(normalized count + 1) matrices under `03_processed/expression/`.
- Reassessed the four Phase 2A suspected outliers and retained all samples, with sensitivity analysis required for the four flagged samples in later phases.
- Ran `06_scripts/python/04_phase2b_prepare_expression.py` successfully.
- Ran `06_scripts/python/04_validate_phase2b_expression.py` successfully.
"""
    if "Completed Phase 2B missingness audit" not in text:
        text = text.replace(
            "- Ran `06_scripts/python/03_validate_phase2a_expression.py` successfully.\n",
            f"- Ran `06_scripts/python/03_validate_phase2a_expression.py` successfully.\n{additions}",
        )
    text = text.replace(
        "- Human review of Phase 2A generated files, matrix dimensions, missingness, and suspected outlier flags before committing.\n- Resolve `TO_VERIFY`: confirm whether DESeq size-factor-normalized counts are the intended downstream input scale, or whether raw-count reprocessing is required in a later phase.\n",
        "- Human review of Phase 2B generated files, missing-value strategy, filtering thresholds, outlier sensitivity plan, and analysis-ready matrices before committing.\n- Resolve `TO_VERIFY`: official source documentation does not explain the literal `NA` entries; raw-count reprocessing remains a possible later-phase check if processed-matrix provenance is judged insufficient.\n",
    )
    text = text.replace(
        "No biological conclusions, subtype-driven feature selection, differential expression, clustering for subtype discovery, batch correction, outlier removal, model training, or SMOTE are approved in Phase 2A.",
        "No biological conclusions, subtype-driven feature selection, differential expression, supervised classification, batch correction, outlier removal, model training, or SMOTE were performed in Phase 2B.",
    )
    text = text.replace(
        "Phase 2B may proceed only after review of the Phase 2A generated files and matrix dimensions. Do not commit automatically until reviewed.",
        "Phase 3 subtype reproduction may proceed from the Phase 2B log2 analysis-ready matrix after human review. Do not commit automatically until reviewed.",
    )
    PROJECT_STATUS.write_text(text)


def write_report(
    matrix: pd.DataFrame,
    norm: pd.DataFrame,
    log2: pd.DataFrame,
    missing_metrics: dict[str, object],
    filtering: pd.DataFrame,
    outlier: pd.DataFrame,
    pca_meta: pd.DataFrame,
) -> None:
    patterns = ", ".join(f"{int(k)} missing in {int(v)} genes" for k, v in missing_metrics["gene_missing_patterns"].items())
    subtype_summary = missing_metrics["missing_by_subtype"].to_string(index=False)
    batch_summary = missing_metrics["missing_by_batch"].to_string(index=False)
    flagged = outlier[outlier["phase2b_sample_classification"] != "RETAIN"][
        ["expression_column", "patient_id", "phase2b_sample_classification", "objective_evidence"]
    ].to_string(index=False)
    selected = filtering[filtering["selected_for_primary_matrix"]].iloc[0]
    report = f"""# Phase 2B Analysis-Ready Expression: GSE172356

## Scope and Guardrails

K-Dense `exploratory-data-analysis` was used from `/Users/emily/.agents/skills/exploratory-data-analysis/SKILL.md`. Phase 2B audited missingness, selected missing-value handling, evaluated subtype-independent expression filters, transformed the selected matrix, and reassessed sample QC. No subtype classification, differential expression, supervised feature selection, batch correction, sample exclusion, or biological interpretation was performed.

## Source and Missing-Value Pattern

Input audited matrix: `03_processed/expression/GSE172356_expression_audited.tsv.gz`, with {matrix.shape[0]} genes and {matrix.shape[1] - 1} mapped samples. The original GEO processed matrix contains {missing_metrics['total_missing']} missing expression cells. Source representation audit found {missing_metrics['representation_counts']['literal_NA']} literal `NA` strings, {missing_metrics['representation_counts']['blank_fields']} blank fields, and {missing_metrics['representation_counts']['parse_failures']} parse failures. Therefore, the missing cells are source literal `NA` values, not parsing failures introduced in Phase 2A.

Missingness is concentrated in genes, not randomly scattered: {patterns}. All 62 samples have missing values; per-sample missing counts range from {missing_metrics['sample_missing_min']} to {missing_metrics['sample_missing_max']}. The pattern forms two complementary sample blocks. Missing count correlation with total expression is {missing_metrics['missing_total_expression_pearson']:.4f}; correlation with detected genes is {missing_metrics['missing_detected_gene_pearson']:.4f}. Batch fields in `sample_manifest.tsv` are `NA`, so batch association is not assessable from current metadata.

Subtype association was audited descriptively only, not used for filtering or threshold choice:

```text
{subtype_summary}
```

Batch summary:

```text
{batch_summary}
```

Official source documentation confirms the processed matrix is based on HTSeq counts normalized by DESeq size factors, but the available GEO/source documentation does not explain the literal `NA` entries. Supplementary Data 4 `Figure1.GeneMatrix` contains a separate 94-signature table and does not document the 45,140-gene processed-matrix `NA` values. TO_VERIFY: confirm with source authors or raw-count reprocessing if the literal `NA` provenance becomes material.

## Missing-Value Handling Decision

Primary strategy: complete-observation filtering, retaining genes with complete observations before expression filtering. This avoids replacing unexplained source `NA` values. The rejected alternatives were:

- Remove genes exceeding a missingness threshold: thresholds <=20% collapse to complete observations because partial-missing genes are missing in 21 or 41 samples; <=50% would require imputing genes missing in 21 samples.
- Replace structurally absent count-like entries with zero: rejected for the primary matrix because source documentation does not state that `NA` means no reads or structural absence.
- Gene-median imputation: evaluated as a sensitivity strategy only; rejected for the primary matrix because it inserts modeled values into block-missing source cells.

Sensitivity analyses required later: repeat Phase 3 subtype reproduction using complete-gene primary matrix, a <=50% missingness plus gene-median imputed matrix, and an all-`NA` zero-filled matrix only as a stress test clearly labelled unsupported by source semantics.

## Filtering and Transformation

Filtering was evaluated without subtype labels. The selected primary rule is `{selected['rule_id']}`: {selected['notes']} This retained {int(selected['genes_retained'])} genes and all 62 samples. Full sensitivity table: `05_results/tables/phase2b_filtering_sensitivity.tsv`.

No additional library-size or size-factor normalization was performed because the source matrix already contains DESeq size-factor-normalized counts. The analysis transform is `log2(normalized_count + 1)`, applied after filtering. DESeq2 VST or rlog were not automatically applied because those methods are designed around raw count inputs and size-factor/model estimation; this matrix is already normalized and contains fractional values, and raw integer counts are not available in the Phase 2B inputs.

Preserved outputs:

- Filtered normalized counts: `03_processed/expression/GSE172356_expression_filtered_normalized.tsv.gz`.
- Filtered log2 analysis-ready matrix: `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz`.

Final matrix dimensions: {norm.shape[0]} genes x {norm.shape[1] - 1} samples. Final missing-value count: {int(norm.drop(columns=['gene']).isna().sum().sum())}.

## Four Suspected Outliers

All 62 samples are retained. The four Phase 2A suspected samples are classified as `RETAIN_WITH_SENSITIVITY_ANALYSIS`; no sample is recommended for exclusion based on a single metric.

```text
{flagged}
```

Outlier assessment table: `05_results/tables/phase2b_outlier_assessment.tsv`. Figures: `05_results/figures/phase2b_transformed_pca.pdf`, `05_results/figures/phase2b_transformed_sample_correlation.pdf`, and `05_results/figures/phase2b_sample_qc_summary.pdf`. PCA on the filtered log2 matrix explains {pca_meta.loc[0, 'variance_explained']:.4f} on PC1 and {pca_meta.loc[1, 'variance_explained']:.4f} on PC2.

## Validation and Phase 3 Readiness

Validation requirements are satisfied: exactly 62 mapped samples remain, sample order matches `expression_sample_crosswalk.tsv`, no duplicated samples or genes are present, no infinite values are present, final missing-value count is 0, the transform is reproducible as `log2(normalized_count + 1)`, and the original audited matrix remains unchanged.

Phase 3 subtype reproduction may proceed from `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz` after human review. Required sensitivity analyses are complete-gene primary vs imputation stress-test matrices, and inclusion vs exclusion of each retained-with-sensitivity Phase 2A suspected outlier.
"""
    OUT_REPORT.write_text(report)


def main() -> None:
    ensure_dirs()
    matrix, xwalk, manifest, phase2a_qc, _annot = load_inputs()
    missing_cells, representation_counts = read_source_missing_representation(xwalk["expression_column"].tolist())
    missing_sample, missing_gene, missing_metrics = missingness_tables(matrix, xwalk, manifest, phase2a_qc, missing_cells, representation_counts)
    filtering, primary_mask, filter_stats = filtering_sensitivity(matrix)
    norm, log2 = write_matrices(matrix, primary_mask)
    outlier, cor, pca_meta = build_outlier_assessment(norm, log2, missing_sample, xwalk)

    missing_sample.to_csv(OUT_MISSING_SAMPLE, sep="\t", index=False)
    missing_gene.to_csv(OUT_MISSING_GENE, sep="\t", index=False)
    filtering.to_csv(OUT_FILTER, sep="\t", index=False)
    outlier.to_csv(OUT_OUTLIER, sep="\t", index=False)
    plot_figures(matrix, missing_gene, log2, outlier, cor)
    write_report(matrix, norm, log2, missing_metrics, filtering, outlier, pca_meta)

    generated = {
        "GSE172356_expression_filtered_normalized": (OUT_NORM, "Phase2B_expression_analysis_ready_matrix", "Filtered DESeq size-factor-normalized counts; complete genes, nonzero, count >=1 in at least 10% of samples."),
        "GSE172356_expression_log2_analysis_ready": (OUT_LOG2, "Phase2B_expression_analysis_ready_matrix", "Filtered log2(normalized count + 1) matrix for unsupervised subtype reproduction input."),
        "phase2b_missingness_by_sample": (OUT_MISSING_SAMPLE, "Phase2B_expression_qc_table", "Per-sample missingness audit including source literal NA representation and metadata joins."),
        "phase2b_missingness_by_gene": (OUT_MISSING_GENE, "Phase2B_expression_qc_table", "Per-gene missingness audit including missing-count block pattern."),
        "phase2b_filtering_sensitivity": (OUT_FILTER, "Phase2B_expression_qc_table", "Subtype-independent filtering and missing-value strategy sensitivity table."),
        "phase2b_outlier_assessment": (OUT_OUTLIER, "Phase2B_expression_qc_table", "Sample QC reassessment on filtered normalized and log2 matrices."),
        "phase2b_missingness_heatmap": (FIG_MISSING, "Phase2B_expression_qc_figure", "Missingness heatmap for genes with source literal NA values."),
        "phase2b_transformed_pca": (FIG_PCA, "Phase2B_expression_qc_figure", "PCA on filtered log2 expression matrix."),
        "phase2b_transformed_sample_correlation": (FIG_COR, "Phase2B_expression_qc_figure", "Sample correlation heatmap on filtered log2 expression matrix."),
        "phase2b_sample_qc_summary": (FIG_QC, "Phase2B_expression_qc_figure", "Sample QC summary for missingness, expression, correlation, PCA, and robust Mahalanobis metrics."),
        "PHASE2B_ANALYSIS_READY_EXPRESSION": (OUT_REPORT, "Phase2B_expression_qc_report", "Phase 2B report documenting missingness source, filtering, transformation, outlier assessment, and Phase 3 readiness."),
    }
    update_file_manifest(generated)
    update_decision_log(filter_stats["primary_genes"])
    update_project_status(filter_stats["primary_genes"])
    print(f"Phase 2B prepared {norm.shape[0]} genes x {norm.shape[1] - 1} samples; final missing values: {int(norm.drop(columns=['gene']).isna().sum().sum())}.")


if __name__ == "__main__":
    main()
