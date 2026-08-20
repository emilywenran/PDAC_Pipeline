#!/usr/bin/env python3
"""Phase 2A audit for the processed GSE172356 expression matrix."""

from __future__ import annotations

import gzip
import hashlib
import math
import os
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "02_data/reference/GSE172356_processed/GSE172356_PDA_gene_expression_matrix.txt.gz"
SAMPLE_MANIFEST = ROOT / "01_metadata/sample_manifest.tsv"
CROSSWALK = ROOT / "01_metadata/rna_microbiome_patient_crosswalk.tsv"
FILE_MANIFEST = ROOT / "01_metadata/file_manifest.tsv"

OUT_MATRIX = ROOT / "03_processed/expression/GSE172356_expression_audited.tsv.gz"
OUT_GENE_ANNOT = ROOT / "03_processed/expression/GSE172356_gene_annotation.tsv"
OUT_SAMPLE_CROSSWALK = ROOT / "01_metadata/expression_sample_crosswalk.tsv"
OUT_SAMPLE_QC = ROOT / "05_results/tables/phase2a_expression_sample_qc.tsv"
OUT_GENE_QC = ROOT / "05_results/tables/phase2a_expression_gene_qc.tsv"
OUT_MAPPING = ROOT / "05_results/tables/phase2a_expression_mapping_summary.tsv"
OUT_REPORT = ROOT / "04_analysis/03_expression_qc/PHASE2A_EXPRESSION_AUDIT.md"
FIG_DIST = ROOT / "05_results/figures/phase2a_expression_distribution.pdf"
FIG_COR = ROOT / "05_results/figures/phase2a_sample_correlation.pdf"
FIG_PCA = ROOT / "05_results/figures/phase2a_expression_pca.pdf"

SOURCE_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE172356&format=file&file="
    "GSE172356%5FPDA%5Fgene%5Fexpression%5Fmatrix%2Etxt%2Egz"
)
GEO_SOFT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172356/soft/GSE172356_family.soft.gz"
TODAY = date.today().isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
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


def ensure_dirs() -> None:
    for path in [
        OUT_MATRIX.parent,
        OUT_GENE_ANNOT.parent,
        OUT_SAMPLE_QC.parent,
        FIG_DIST.parent,
        OUT_REPORT.parent,
    ]:
        path.mkdir(parents=True, exist_ok=True)
    mpl_dir = ROOT / ".matplotlib"
    mpl_dir.mkdir(exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_dir)


def read_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(SOURCE, sep="\t", dtype=str, compression="gzip")
    if raw.columns[0] != "gene":
        raise ValueError(f"Expected first column 'gene', observed {raw.columns[0]!r}")
    genes = raw["gene"].astype(str)
    numeric = raw.drop(columns=["gene"]).apply(pd.to_numeric, errors="coerce")
    audited = pd.concat([genes.rename("gene"), numeric], axis=1)
    return raw, audited


def build_expression_crosswalk(expr_cols: list[str]) -> pd.DataFrame:
    manifest = pd.read_csv(SAMPLE_MANIFEST, sep="\t", dtype=str)
    patient_xwalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)

    alias_map = {}
    for _, row in manifest.iterrows():
        notes = row.get("notes", "")
        alias = "NA"
        marker = "Microbiome alias "
        if marker in notes:
            alias = notes.split(marker, 1)[1].split(";", 1)[0].strip()
        alias_map[alias] = row

    tumor_map = patient_xwalk.set_index("patient_id")["tumor_number"].to_dict()
    records = []
    for col in expr_cols:
        row = alias_map.get(col)
        if row is None:
            records.append(
                {
                    "expression_column": col,
                    "geo_sample_id": "NA",
                    "patient_id": "NA",
                    "tumor_number": "NA",
                    "subtype_original": "NA",
                    "mapping_status": "UNMAPPED",
                    "notes": "No matching microbiome alias in sample_manifest notes.",
                }
            )
        else:
            records.append(
                {
                    "expression_column": col,
                    "geo_sample_id": row["geo_sample_id"],
                    "patient_id": row["patient_id"],
                    "tumor_number": tumor_map.get(row["patient_id"], "NA"),
                    "subtype_original": row["subtype_original"],
                    "mapping_status": "MAPPED",
                    "notes": "Matched expression column to manifest microbiome alias; RNA GEO ID from manifest retained.",
                }
            )
    return pd.DataFrame.from_records(records)


def audit_values(raw: pd.DataFrame, audited: pd.DataFrame) -> dict[str, object]:
    value_strings = raw.drop(columns=["gene"])
    expr = audited.drop(columns=["gene"])
    non_numeric = int(expr.isna().sum().sum() - value_strings.isna().sum().sum())
    finite = np.isfinite(expr.to_numpy(dtype=float))
    return {
        "n_rows": audited.shape[0],
        "n_cols": audited.shape[1],
        "n_genes": audited.shape[0],
        "n_samples": expr.shape[1],
        "gene_id_type": "gene_symbol_or_gene_name",
        "gene_ids_unique": bool(audited["gene"].is_unique),
        "duplicate_gene_rows": int(audited["gene"].duplicated(keep=False).sum()),
        "duplicate_gene_ids": int(audited["gene"].duplicated().sum()),
        "duplicate_samples": int(pd.Index(expr.columns).duplicated().sum()),
        "missing_values": int(expr.isna().sum().sum()),
        "non_numeric_values": non_numeric,
        "infinite_values": int((~finite & ~expr.isna().to_numpy()).sum()),
        "negative_values": int((expr < 0).sum().sum()),
        "zero_values": int((expr == 0).sum().sum()),
        "all_zero_genes": int((expr.fillna(0).sum(axis=1) == 0).sum()),
        "fractional_value_present": bool(((expr.dropna().to_numpy() % 1) != 0).any()),
        "min_value": float(expr.min().min()),
        "max_value": float(expr.max().max()),
    }


def build_qc(audited: pd.DataFrame, xwalk: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr = audited.drop(columns=["gene"])
    genes = audited["gene"]
    sample_qc = pd.DataFrame(
        {
            "expression_column": expr.columns,
            "total_expression": expr.sum(axis=0).values,
            "mean_expression": expr.mean(axis=0).values,
            "median_expression": expr.median(axis=0).values,
            "sd_expression": expr.std(axis=0).values,
            "min_expression": expr.min(axis=0).values,
            "max_expression": expr.max(axis=0).values,
            "detected_genes_gt0": (expr > 0).sum(axis=0).values,
            "zero_gene_count": (expr == 0).sum(axis=0).values,
            "zero_proportion": (expr == 0).mean(axis=0).values,
            "missing_count": expr.isna().sum(axis=0).values,
            "negative_count": (expr < 0).sum(axis=0).values,
            "infinite_count": ((~np.isfinite(expr.to_numpy(dtype=float))) & (~expr.isna().to_numpy())).sum(axis=0),
        }
    )

    cor = expr.corr(method="pearson")
    cor_no_diag = cor.mask(np.eye(cor.shape[0], dtype=bool))
    sample_qc["mean_sample_correlation"] = cor_no_diag.mean(axis=0).reindex(sample_qc["expression_column"]).values
    sample_qc["min_sample_correlation"] = cor_no_diag.min(axis=0).reindex(sample_qc["expression_column"]).values

    sample_qc["total_expression_robust_z"] = robust_z(sample_qc["total_expression"]).values
    sample_qc["detected_genes_robust_z"] = robust_z(sample_qc["detected_genes_gt0"]).values
    q1 = sample_qc["mean_sample_correlation"].quantile(0.25)
    q3 = sample_qc["mean_sample_correlation"].quantile(0.75)
    iqr = q3 - q1
    corr_lower_fence = q1 - 3 * iqr

    pca_input = np.log10(expr.T + 1)
    # Visualization-only imputation: preserve missing values in all tabular outputs.
    pca_input = pca_input.fillna(pca_input.median(axis=0)).fillna(0)
    pca_scaled = StandardScaler().fit_transform(pca_input)
    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(pca_scaled)
    sample_qc["pc1"] = pcs[:, 0]
    sample_qc["pc2"] = pcs[:, 1]
    pc1_z = robust_z(sample_qc["pc1"])
    pc2_z = robust_z(sample_qc["pc2"])
    sample_qc["pc1_robust_z"] = pc1_z.values
    sample_qc["pc2_robust_z"] = pc2_z.values

    flags = []
    for _, row in sample_qc.iterrows():
        reasons = []
        if abs(row["total_expression_robust_z"]) > 3.5:
            reasons.append("total_expression_robust_z_abs_gt_3.5")
        if row["detected_genes_robust_z"] < -3.5:
            reasons.append("detected_genes_robust_z_lt_-3.5")
        if row["mean_sample_correlation"] < corr_lower_fence:
            reasons.append("mean_correlation_below_Q1_minus_3IQR")
        if abs(row["pc1_robust_z"]) > 3.5 or abs(row["pc2_robust_z"]) > 3.5:
            reasons.append("pc1_or_pc2_robust_z_abs_gt_3.5")
        flags.append(";".join(reasons) if reasons else "none")
    sample_qc["outlier_flag"] = flags

    sample_qc = sample_qc.merge(xwalk, on="expression_column", how="left")

    gene_qc = pd.DataFrame(
        {
            "gene": genes,
            "row_index_1based": np.arange(1, len(genes) + 1),
            "mean_expression": expr.mean(axis=1).values,
            "median_expression": expr.median(axis=1).values,
            "sd_expression": expr.std(axis=1).values,
            "min_expression": expr.min(axis=1).values,
            "max_expression": expr.max(axis=1).values,
            "detected_sample_count_gt0": (expr > 0).sum(axis=1).values,
            "zero_sample_count": (expr == 0).sum(axis=1).values,
            "zero_proportion": (expr == 0).mean(axis=1).values,
            "missing_count": expr.isna().sum(axis=1).values,
            "negative_count": (expr < 0).sum(axis=1).values,
            "infinite_count": ((~np.isfinite(expr.to_numpy(dtype=float))) & (~expr.isna().to_numpy())).sum(axis=1),
            "is_all_zero": (expr.fillna(0).sum(axis=1) == 0).values,
            "is_duplicate_gene_id": genes.duplicated(keep=False).values,
        }
    )

    mapping_summary = (
        xwalk.groupby("mapping_status", dropna=False)
        .size()
        .reset_index(name="sample_count")
        .sort_values("mapping_status")
    )
    pca_meta = pd.DataFrame(
        {
            "component": ["PC1", "PC2"],
            "variance_explained": pca.explained_variance_ratio_,
        }
    )
    return sample_qc, gene_qc, mapping_summary, pca_meta


def write_outputs(audited: pd.DataFrame, xwalk: pd.DataFrame, sample_qc: pd.DataFrame, gene_qc: pd.DataFrame, mapping_summary: pd.DataFrame) -> None:
    audited.to_csv(OUT_MATRIX, sep="\t", index=False, compression="gzip")
    xwalk.to_csv(OUT_SAMPLE_CROSSWALK, sep="\t", index=False)
    sample_qc.to_csv(OUT_SAMPLE_QC, sep="\t", index=False)
    gene_qc.to_csv(OUT_GENE_QC, sep="\t", index=False)
    mapping_summary.to_csv(OUT_MAPPING, sep="\t", index=False)

    annot = pd.DataFrame(
        {
            "gene": audited["gene"],
            "row_index_1based": np.arange(1, audited.shape[0] + 1),
            "gene_identifier_type": "gene_symbol_or_gene_name",
            "is_duplicate_gene_id": audited["gene"].duplicated(keep=False),
            "notes": "Identifier taken directly from GEO matrix first column; no external annotation file was required or applied.",
        }
    )
    annot.to_csv(OUT_GENE_ANNOT, sep="\t", index=False)


def plot_figures(audited: pd.DataFrame, sample_qc: pd.DataFrame) -> None:
    expr = audited.drop(columns=["gene"])

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].hist(sample_qc["total_expression"], bins=20, color="#3b6ea8", edgecolor="white")
    axes[0].set_xlabel("Total expression")
    axes[0].set_ylabel("Samples")
    axes[0].set_title("Per-sample total")
    axes[1].hist(sample_qc["zero_proportion"], bins=20, color="#5b8f6b", edgecolor="white")
    axes[1].set_xlabel("Zero proportion")
    axes[1].set_title("Per-sample zeros")
    vals = np.log10(expr.to_numpy(dtype=float).ravel() + 1)
    vals = vals[np.isfinite(vals)]
    axes[2].hist(vals, bins=80, color="#8a6f3d", edgecolor="white")
    axes[2].set_xlabel("log10(expression + 1)")
    axes[2].set_title("Value distribution")
    fig.tight_layout()
    fig.savefig(FIG_DIST)
    plt.close(fig)

    cor = expr.corr(method="pearson")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(cor, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Sample-sample Pearson correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(FIG_COR)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    status_colors = {"none": "#2f6f8f"}
    flagged = sample_qc["outlier_flag"] != "none"
    ax.scatter(sample_qc.loc[~flagged, "pc1"], sample_qc.loc[~flagged, "pc2"], c=status_colors["none"], s=35, label="not flagged")
    if flagged.any():
        ax.scatter(sample_qc.loc[flagged, "pc1"], sample_qc.loc[flagged, "pc2"], c="#b33a3a", s=45, label="flagged")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA QC, log10(expression + 1), gene-scaled")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_PCA)
    plt.close(fig)


def update_file_manifest(paths: list[Path]) -> None:
    manifest = pd.read_csv(FILE_MANIFEST, sep="\t", dtype=str)
    new_file_ids = {
        "GSE172356_PDA_gene_expression_matrix",
        "expression_sample_crosswalk",
        "GSE172356_expression_audited",
        "GSE172356_gene_annotation",
        "phase2a_expression_sample_qc",
        "phase2a_expression_gene_qc",
        "phase2a_expression_mapping_summary",
        "PHASE2A_EXPRESSION_AUDIT",
    }
    manifest = manifest[~manifest["file_id"].isin(new_file_ids)]

    entries = [
        {
            "file_id": "GSE172356_PDA_gene_expression_matrix",
            "dataset": "GSE172356_processed",
            "sample_id": "",
            "data_type": "processed_expression_matrix_source",
            "local_path": str(SOURCE),
            "source_url_or_accession": SOURCE_URL,
            "file_size": str(SOURCE.stat().st_size),
            "md5": f"md5:{md5(SOURCE)};sha256:{sha256(SOURCE)}",
            "download_date": TODAY,
            "processing_status": "downloaded_verified_Phase2A",
            "notes": "Official GEO Series supplementary file; gzip compression; original filename GSE172356_PDA_gene_expression_matrix.txt.gz; GEO SOFT processing evidence reports HTSeq counts normalized by DESeq SizeFactors.",
        }
    ]
    generated = {
        "expression_sample_crosswalk": OUT_SAMPLE_CROSSWALK,
        "GSE172356_expression_audited": OUT_MATRIX,
        "GSE172356_gene_annotation": OUT_GENE_ANNOT,
        "phase2a_expression_sample_qc": OUT_SAMPLE_QC,
        "phase2a_expression_gene_qc": OUT_GENE_QC,
        "phase2a_expression_mapping_summary": OUT_MAPPING,
        "PHASE2A_EXPRESSION_AUDIT": OUT_REPORT,
    }
    for file_id, path in generated.items():
        entries.append(
            {
                "file_id": file_id,
                "dataset": "GSE172356_processed",
                "sample_id": "",
                "data_type": "Phase2A_expression_audit_output",
                "local_path": str(path),
                "source_url_or_accession": "derived_from_GSE172356_PDA_gene_expression_matrix",
                "file_size": str(path.stat().st_size),
                "md5": f"sha256:{sha256(path)}",
                "download_date": TODAY,
                "processing_status": "generated_Phase2A",
                "notes": "Generated by 06_scripts/python/03_phase2a_expression_audit.py; source matrix scale preserved where applicable.",
            }
        )

    manifest = pd.concat([manifest, pd.DataFrame(entries)], ignore_index=True)
    manifest.to_csv(FILE_MANIFEST, sep="\t", index=False)


def write_report(metrics: dict[str, object], sample_qc: pd.DataFrame, mapping_summary: pd.DataFrame, pca_meta: pd.DataFrame) -> None:
    mapped = int(mapping_summary.loc[mapping_summary["mapping_status"] == "MAPPED", "sample_count"].sum())
    flagged = sample_qc[sample_qc["outlier_flag"] != "none"][
        ["expression_column", "geo_sample_id", "patient_id", "outlier_flag"]
    ]
    flagged_text = "None by prespecified criteria."
    if not flagged.empty:
        flagged_text = flagged.to_string(index=False)
    report = f"""# Phase 2A Expression Audit: GSE172356

## Source and Skill Usage

- Required K-Dense skill loaded: `/Users/emily/.agents/skills/exploratory-data-analysis/SKILL.md`.
- Applicable skill workflow: tab-delimited scientific data ingestion with row/column counts, type inference, missingness checks, duplicate checks, outlier detection, and correlation/PCA summaries.
- Official processed expression source: `GSE172356_PDA_gene_expression_matrix.txt.gz`.
- Source URL/accession: `{SOURCE_URL}`; GEO Series accession `GSE172356`.
- GEO SOFT evidence URL: `{GEO_SOFT_URL}`.
- Download date: `{TODAY}`.
- Source file size: `{SOURCE.stat().st_size}` bytes.
- Source MD5: `{md5(SOURCE)}`.
- Source SHA256: `{sha256(SOURCE)}`.
- Compression format: gzip.

## Matrix Structure

- Final audited matrix: `{OUT_MATRIX.relative_to(ROOT)}`.
- Matrix dimensions: {metrics['n_genes']} gene rows x {metrics['n_samples']} expression samples, plus one `gene` identifier column.
- Orientation: genes are rows; samples are columns.
- Sample identifier format: `YX...T` tumor aliases, for example `YX15261T`.
- Gene identifier type: gene symbols/gene names from the GEO matrix first column.
- Gene identifiers unique: {metrics['gene_ids_unique']}.
- Duplicate gene rows: {metrics['duplicate_gene_rows']}; duplicate gene IDs beyond first occurrence: {metrics['duplicate_gene_ids']}.
- Duplicate samples: {metrics['duplicate_samples']}.

## Expression Unit

Expression unit is **DESeq size-factor-normalized counts**. Supporting evidence: official GEO SOFT metadata states reads were aligned to GRCh38 with HISAT2 v2.1.0, raw transcript counts were calculated using HTSeq v0.12.4, and normalization was performed by SizeFactors using DESeq v1.24.0 in R. The matrix contains non-integer fractional values, supporting normalized rather than raw integer counts. It is not FPKM, TPM, CPM, or log-transformed based on available GEO evidence.

## Value Integrity

- Missing values: {metrics['missing_values']}.
- Non-numeric values in expression cells: {metrics['non_numeric_values']}.
- Infinite values: {metrics['infinite_values']}.
- Negative values: {metrics['negative_values']}.
- Zero values: {metrics['zero_values']}.
- All-zero genes: {metrics['all_zero_genes']}.
- Minimum value: {metrics['min_value']}.
- Maximum value: {metrics['max_value']}.

## Sample Mapping

- Expression columns evaluated: {metrics['n_samples']}.
- Successfully mapped samples: {mapped} / 62.
- Mapping output: `{OUT_SAMPLE_CROSSWALK.relative_to(ROOT)}`.
- Mapping rule: expression `YX...T` aliases were matched to the 62-patient manifest aliases recorded in `sample_manifest.tsv`; GEO sample IDs and patient IDs were retained from the finalized manifest.

## Descriptive QC

Outputs:

- Per-sample QC: `{OUT_SAMPLE_QC.relative_to(ROOT)}`.
- Per-gene QC: `{OUT_GENE_QC.relative_to(ROOT)}`.
- Mapping summary: `{OUT_MAPPING.relative_to(ROOT)}`.
- Distribution figure: `{FIG_DIST.relative_to(ROOT)}`.
- Sample correlation figure: `{FIG_COR.relative_to(ROOT)}`.
- PCA figure: `{FIG_PCA.relative_to(ROOT)}`.

PCA was performed for QC visualization only on `log10(expression + 1)` values after visualization-only gene-median imputation for missing cells and gene-wise scaling. This transformed PCA input was not written as the audited matrix and was not used for subtype discovery, feature selection, differential expression, or biological interpretation. PC1 explained {pca_meta.loc[0, 'variance_explained']:.4f} and PC2 explained {pca_meta.loc[1, 'variance_explained']:.4f} of variance in this QC space.

## Prespecified Extreme-Sample Criteria

Samples were flagged if any of the following objective criteria were met:

- absolute robust z-score for total expression > 3.5;
- robust z-score for detected genes < -3.5;
- mean sample correlation below Q1 - 3 x IQR;
- absolute robust z-score for PC1 or PC2 > 3.5 in the QC PCA space.

Suspected outliers:

{flagged_text}

No outliers were removed.

## Transformations and Phase 2B Readiness

The audited matrix preserves the original GEO numerical scale and requires no transformation for archiving or sample mapping. For downstream PCA, clustering, subtype reproduction, and regression diagnostics, variance-stabilizing or log-like transformations should be considered and explicitly scripted in the appropriate later phase because the current normalized-count scale is strongly right-skewed. No normalization, batch correction, sample removal, subtype-driven gene selection, differential expression, or biological interpretation was performed in Phase 2A.

Phase 2B may proceed after human review of the generated files and matrix dimensions. TO_VERIFY: confirm with the project owner that DESeq size-factor-normalized counts are the intended input scale for downstream subtype reproduction, or whether raw-count reprocessing from FASTQ is required in a later phase.
"""
    OUT_REPORT.write_text(report)


def main() -> None:
    ensure_dirs()
    raw, audited = read_matrix()
    expr_cols = audited.columns[1:].tolist()
    xwalk = build_expression_crosswalk(expr_cols)
    metrics = audit_values(raw, audited)
    sample_qc, gene_qc, mapping_summary, pca_meta = build_qc(audited, xwalk)
    write_outputs(audited, xwalk, sample_qc, gene_qc, mapping_summary)
    plot_figures(audited, sample_qc)
    write_report(metrics, sample_qc, mapping_summary, pca_meta)
    update_file_manifest(
        [
            OUT_SAMPLE_CROSSWALK,
            OUT_MATRIX,
            OUT_GENE_ANNOT,
            OUT_SAMPLE_QC,
            OUT_GENE_QC,
            OUT_MAPPING,
            OUT_REPORT,
        ]
    )
    print(f"Audited {metrics['n_genes']} genes x {metrics['n_samples']} samples; mapped {int(mapping_summary.loc[mapping_summary['mapping_status'] == 'MAPPED', 'sample_count'].sum())}/62 samples.")


if __name__ == "__main__":
    main()
