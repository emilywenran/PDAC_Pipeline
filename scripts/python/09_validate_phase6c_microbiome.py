#!/usr/bin/env python3
"""Validate Phase 6C analysis-ready microbiome matrices and reports."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
AUDITED = ROOT / "03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz"
PRIMARY_FILTERED = ROOT / "03_processed/microbiome/PRJNA719915_genus_primary_filtered.tsv.gz"
PRIMARY_CLR = ROOT / "03_processed/microbiome/PRJNA719915_genus_primary_CLR.tsv.gz"
PRIMARY_DIST = ROOT / "03_processed/microbiome/PRJNA719915_primary_aitchison_distance.tsv.gz"
CROSSWALK = ROOT / "01_metadata/microbiome_sample_crosswalk.tsv"
TABLE_DIR = ROOT / "05_results/tables"
FIGURE_DIR = ROOT / "05_results/figures"
REPORT = ROOT / "04_analysis/04_microbiome_qc/PHASE6C_ANALYSIS_READY_MICROBIOME.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    audited_before = sha256(AUDITED)
    audited = pd.read_csv(AUDITED, sep="\t", index_col=0)
    filtered = pd.read_csv(PRIMARY_FILTERED, sep="\t", index_col=0)
    clr = pd.read_csv(PRIMARY_CLR, sep="\t", index_col=0)
    dist = pd.read_csv(PRIMARY_DIST, sep="\t", index_col=0)
    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)

    require(audited.shape == (365, 62), f"Unexpected audited matrix shape: {audited.shape}")
    require(filtered.shape == (122, 62), f"Unexpected primary filtered shape: {filtered.shape}")
    require(clr.shape == (122, 62), f"Unexpected primary CLR shape: {clr.shape}")
    require(dist.shape == (62, 62), f"Unexpected distance shape: {dist.shape}")
    require(list(filtered.columns) == list(audited.columns), "Primary filtered sample ordering changed")
    require(list(clr.columns) == list(audited.columns), "CLR sample ordering changed")
    require(list(dist.index) == list(audited.columns), "Distance row ordering changed")
    require(list(dist.columns) == list(audited.columns), "Distance column ordering changed")
    require(list(crosswalk["microbiome_matrix_sample"]) == list(audited.columns), "Crosswalk order mismatch")
    require(audited.index.is_unique and filtered.index.is_unique and clr.index.is_unique, "Duplicated genera")
    require(audited.columns.is_unique and filtered.columns.is_unique and clr.columns.is_unique, "Duplicated samples")
    require(np.isfinite(clr.to_numpy()).all(), "Missing or infinite CLR values")
    require(np.allclose(clr.sum(axis=0).to_numpy(), 0.0, atol=1e-8), "CLR columns do not sum to zero")
    d = dist.to_numpy(dtype=float)
    require(np.isfinite(d).all(), "Missing or infinite distances")
    require((d >= -1e-12).all(), "Negative distances present")
    require(np.allclose(d, d.T, atol=1e-10), "Distance matrix is not symmetric")
    require(np.allclose(np.diag(d), 0.0, atol=1e-10), "Distance diagonal is not zero")
    require(sha256(AUDITED) == audited_before, "Original audited matrix changed during validation")

    required_tables = [
        "phase6c_retained_taxa_with_contamination_flags.tsv",
        "phase6c_processed_sample_qc.tsv",
        "phase6c_processed_taxon_qc.tsv",
        "phase6c_matrix_inventory.tsv",
        "phase6c_preprocessing_sensitivity_concordance.tsv",
    ]
    required_figures = [
        "phase6c_primary_CLR_heatmap.pdf",
        "phase6c_primary_Aitchison_PCoA.pdf",
        "phase6c_filtering_sensitivity_summary.pdf",
        "phase6c_pseudocount_sensitivity.pdf",
        "phase6c_contamination_flag_summary.pdf",
        "phase6c_sample_depth_proxy.pdf",
    ]
    for name in required_tables:
        path = TABLE_DIR / name
        require(path.exists() and path.stat().st_size > 0, f"Missing/empty table: {name}")
    for name in required_figures:
        path = FIGURE_DIR / name
        require(path.exists() and path.stat().st_size > 0, f"Missing/empty figure: {name}")
    require(REPORT.exists() and REPORT.stat().st_size > 0, "Missing/empty Phase 6C report")

    report_text = REPORT.read_text().lower()
    forbidden = [
        "survival analysis was performed",
        "differential abundance was performed",
        "host-microbiome correlation was performed",
        "pathway analysis was performed",
        "target prioritization was performed",
    ]
    for phrase in forbidden:
        require(phrase not in report_text, f"Forbidden downstream claim in report: {phrase}")

    print("Phase 6C validation passed")
    print(f"Primary matrix: {filtered.shape[0]} genera x {filtered.shape[1]} samples")
    print(f"CLR column max abs sum: {np.abs(clr.sum(axis=0)).max():.3e}")
    print(f"Aitchison distance range: {d.min():.6f} to {d.max():.6f}")


if __name__ == "__main__":
    main()
