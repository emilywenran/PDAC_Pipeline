#!/usr/bin/env python3
"""Register Phase 9B3R corrected spatial-validation artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/emily/thesis/PDAC")
MANIFEST = ROOT / "01_metadata/file_manifest.tsv"
DATE = "2026-07-04"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def artifact(file_id: str, rel: str, data_type: str, notes: str) -> dict[str, str]:
    return {
        "file_id": file_id,
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": data_type,
        "local_path": str(ROOT / rel),
        "source_url_or_accession": "derived_Phase9B3R_corrected_reanalysis",
        "download_date": DATE,
        "processing_status": "generated_Phase9B3R",
        "notes": notes,
    }


ARTIFACTS = [
    artifact("PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS", "04_analysis/09_external_validation/PHASE9B3R_CORRECTED_SPATIAL_VALIDATION_RESULTS.md", "analysis_report", "Corrected Phase 9B3R spatial validation report."),
    artifact("PHASE9B3R_CORRECTION_LOG", "04_analysis/09_external_validation/PHASE9B3R_CORRECTION_LOG.md", "analysis_report", "Correction log closing Phase 9B3C findings."),
    artifact("PHASE9B3B_SPATIAL_VALIDATION_RESULTS_SUPERSEDED", "04_analysis/09_external_validation/PHASE9B3B_SPATIAL_VALIDATION_RESULTS.md", "analysis_report", "Original Phase 9B3B report marked SUPERSEDED_BY_PHASE9B3R."),
    artifact("phase9b3r_dataset_qc", "05_results/tables/phase9b3r_dataset_qc.tsv", "result_table", "Corrected Phase 9B3R dataset QC table."),
    artifact("phase9b3r_feature_coverage", "05_results/tables/phase9b3r_feature_coverage.tsv", "result_table", "Feature coverage with corrected eligibility gating."),
    artifact("phase9b3r_hwang_naive_models", "05_results/tables/phase9b3r_hwang_naive_models.tsv", "result_table", "Corrected Hwang naive model results."),
    artifact("phase9b3r_hwang_treated_models", "05_results/tables/phase9b3r_hwang_treated_models.tsv", "result_table", "Corrected Hwang treated sensitivity model results."),
    artifact("phase9b3r_moncada_exploratory_results", "05_results/tables/phase9b3r_moncada_exploratory_results.tsv", "result_table", "Moncada exploratory consistency results."),
    artifact("phase9b3r_negative_control_results", "05_results/tables/phase9b3r_negative_control_results.tsv", "result_table", "Corrected negative-control summary table."),
    artifact("phase9b3r_negative_control_null_distributions", "05_results/tables/phase9b3r_negative_control_null_distributions.tsv", "result_table", "Iteration-level null distributions for corrected negative controls."),
    artifact("phase9b3r_spatial_evidence", "05_results/tables/phase9b3r_spatial_evidence.tsv", "result_table", "Programmatically derived corrected spatial evidence categories."),
    artifact("phase9b3r_cross_cohort_synthesis", "05_results/tables/phase9b3r_cross_cohort_synthesis.tsv", "result_table", "Corrected descriptive cross-cohort synthesis."),
    artifact("phase9b3r_runtime_validation", "05_results/tables/phase9b3r_runtime_validation.tsv", "runtime_validation", "Runtime guardrail checks for Phase 9B3R."),
    artifact("phase9b3r_hwang_primary_models", "05_results/figures/phase9b3r_hwang_primary_models.pdf", "figure", "Corrected primary Hwang model coefficient figure excluding invalid model rows."),
    artifact("16_phase9b3r_spatial_validation", "06_scripts/python/16_phase9b3r_spatial_validation.py", "analysis_script", "Corrected Phase 9B3R spatial validation execution script."),
    artifact("16_validate_phase9b3r_spatial", "06_scripts/python/16_validate_phase9b3r_spatial.py", "validation_script", "Phase 9B3R corrected spatial validation validator."),
    artifact("test_phase9b3r_spatial", "06_scripts/python/test_phase9b3r_spatial.py", "unit_test", "Executable pytest tests for Phase 9B3R repair guardrails."),
    artifact("update_phase9b3r_manifest", "06_scripts/python/update_phase9b3r_manifest.py", "manifest_update_script", "Registers Phase 9B3R artifacts in the file manifest."),
    artifact("00_admin__PROJECT_STATUS_md", "00_admin/PROJECT_STATUS.md", "project_documentation", "Project status updated for Phase 9B3R corrected spatial validation."),
    artifact("00_admin__SKILL_USAGE_LOG_tsv", "00_admin/SKILL_USAGE_LOG.tsv", "project_documentation", "Skill usage log updated for Phase 9B3R corrected spatial validation."),
    artifact("09_docs__planning__DECISION_LOG_md", "09_docs/planning/DECISION_LOG.md", "project_documentation", "Decision log updated with Phase 9B3R completion decision."),
]


def main() -> int:
    df = pd.read_csv(MANIFEST, sep="\t")
    for col in df.columns:
        if df[col].dtype != object:
            df[col] = df[col].astype(object)
    added = 0
    updated = 0
    for item in ARTIFACTS:
        path = Path(item["local_path"])
        if not path.exists():
            raise FileNotFoundError(path)
        size = path.stat().st_size
        row = {col: "" for col in df.columns}
        row.update(item)
        row["file_size"] = float(size)
        row["size_bytes"] = size
        row["md5"] = sha256(path)
        row["checksum"] = row["md5"]
        match = df.index[df["local_path"].eq(str(path))].tolist()
        if match:
            for col in df.columns:
                df.at[match[0], col] = row[col]
            updated += 1
        else:
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
            added += 1
    df.fillna("", inplace=True)
    df.to_csv(MANIFEST, sep="\t", index=False)
    print(f"Phase 9B3R manifest updated: {added} added, {updated} updated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
