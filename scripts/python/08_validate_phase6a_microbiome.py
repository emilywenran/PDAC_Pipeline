#!/usr/bin/env python3
"""Validate Phase 6A PRJNA719915 processed microbiome audit outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
MATRIX = ROOT / "03_processed/microbiome/PRJNA719915_microbiome_abundance_audited.tsv.gz"
CROSSWALK = ROOT / "01_metadata/microbiome_sample_crosswalk.tsv"
PATIENT_CROSSWALK = ROOT / "01_metadata/rna_microbiome_patient_crosswalk.tsv"
TABLES = [
    ROOT / "05_results/tables/phase6a_microbiome_sample_qc.tsv",
    ROOT / "05_results/tables/phase6a_microbiome_taxon_qc.tsv",
    ROOT / "05_results/tables/phase6a_taxon_prevalence.tsv",
    ROOT / "05_results/tables/phase6a_potential_contaminant_flags.tsv",
]
FIGURES = [
    ROOT / "05_results/figures/phase6a_microbiome_abundance_distribution.pdf",
    ROOT / "05_results/figures/phase6a_taxon_prevalence.pdf",
    ROOT / "05_results/figures/phase6a_sample_detection_summary.pdf",
    ROOT / "05_results/figures/phase6a_microbiome_ordination.pdf",
    ROOT / "05_results/figures/phase6a_sample_distance_heatmap.pdf",
]
REPORT = ROOT / "04_analysis/04_microbiome_qc/PHASE6A_MICROBIOME_DATA_AUDIT.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    require(MATRIX.exists(), f"Missing matrix: {MATRIX}")
    matrix = pd.read_csv(MATRIX, sep="\t", index_col=0)
    values = matrix.to_numpy(dtype=float)

    require(matrix.shape == (365, 62), f"Unexpected matrix shape: {matrix.shape}")
    require(matrix.index.is_unique, "Duplicated taxonomic identifiers")
    require(matrix.columns.is_unique, "Duplicated matrix sample columns")
    require(np.isfinite(values).all(), "Missing or infinite values present")
    require(not np.isnan(values).any(), "NaN values present")
    require((values >= 0).all(), "Negative abundance values present")
    require(matrix.index.name == "taxon", "Expected feature rows indexed by taxon")

    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    expected_cols = [
        "microbiome_matrix_sample",
        "patient_id",
        "tumor_number",
        "microbiome_biosample_id",
        "microbiome_run_id",
        "mapping_status",
        "notes",
    ]
    require(list(crosswalk.columns) == expected_cols, "Crosswalk column mismatch")
    require(len(crosswalk) == 62, f"Crosswalk row count is {len(crosswalk)}")
    require(set(crosswalk["microbiome_matrix_sample"]) == set(matrix.columns), "Crosswalk/sample mismatch")
    require(crosswalk["microbiome_matrix_sample"].is_unique, "Duplicated crosswalk matrix samples")
    require(crosswalk["patient_id"].is_unique, "More than one microbiome profile per patient")
    require((crosswalk["mapping_status"] == "VERIFIED").all(), "Unmatched patient mappings present")
    require(crosswalk["patient_id"].nunique() == 62, "Expected 62 unique patients")

    patients = pd.read_csv(PATIENT_CROSSWALK, sep="\t", dtype=str)
    require(set(crosswalk["patient_id"]) == set(patients["patient_id"]), "Not all patients matched")
    require(set(crosswalk["microbiome_run_id"]) == set(patients["microbiome_run_id"]), "Run IDs do not reconcile")

    sample_qc = pd.read_csv(TABLES[0], sep="\t")
    taxon_qc = pd.read_csv(TABLES[1], sep="\t")
    require(len(sample_qc) == 62, "Sample QC row count mismatch")
    require(len(taxon_qc) == 365, "Taxon QC row count mismatch")
    require((sample_qc["detected_taxa"] == (matrix > 0).sum(axis=0).values).all(), "Detected taxa mismatch")
    require(np.allclose(sample_qc["total_abundance"], matrix.sum(axis=0).values), "Sample totals mismatch")
    require(np.allclose(taxon_qc["total_abundance"], matrix.sum(axis=1).values), "Taxon totals mismatch")

    for path in TABLES + FIGURES + [REPORT]:
        require(path.exists(), f"Missing required output: {path}")
        require(path.stat().st_size > 0, f"Empty required output: {path}")

    report_text = REPORT.read_text()
    forbidden = [
        "differential abundance was performed",
        "survival analysis was performed",
        "host-microbiome correlation was performed",
        "decontam prevalence analysis was performed",
    ]
    for phrase in forbidden:
        require(phrase not in report_text, f"Forbidden claim in report: {phrase}")

    print("Phase 6A validation passed")
    print(f"Matrix: {matrix.shape[0]} genera x {matrix.shape[1]} tumor samples")
    print(f"Patients mapped: {crosswalk['patient_id'].nunique()}")
    print(f"Zero fraction: {(values == 0).mean():.4f}")


if __name__ == "__main__":
    main()
