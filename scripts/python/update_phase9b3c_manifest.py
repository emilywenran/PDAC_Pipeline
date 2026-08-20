#!/usr/bin/env python3
"""Register Phase 9B3C Spatial Validation Independent Review files in the file manifest."""

import os
import hashlib
import pandas as pd
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
MANIFEST_PATH = ROOT / "01_metadata/file_manifest.tsv"

def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()

def main():
    if not MANIFEST_PATH.exists():
        print(f"Error: manifest file not found: {MANIFEST_PATH}")
        return
        
    df = pd.read_csv(MANIFEST_PATH, sep="\t")
    
    new_files = [
        {
            "file_id": "PHASE9B3C_SPATIAL_VALIDATION_INDEPENDENT_REVIEW",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9B3C_SPATIAL_VALIDATION_INDEPENDENT_REVIEW.md"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Independent review report auditing Phase 9B3B spatial validation findings."
        },
        {
            "file_id": "phase9b3c_cohort_count_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_cohort_count_audit.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Cohort count audit checking planned vs actual patients, sections, and segments."
        },
        {
            "file_id": "phase9b3c_roi_pairing_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_roi_pairing_audit.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "ROI pairing audit verifying patient-replicate and pairing preservation rules."
        },
        {
            "file_id": "phase9b3c_model_reproduction_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_model_reproduction_audit.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Model reproduction audit tracking coefficient, convergence, and significance of primary models."
        },
        {
            "file_id": "phase9b3c_feature_eligibility_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_feature_eligibility_audit.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Feature eligibility audit verifying gene coverage and prospective exclusion status."
        },
        {
            "file_id": "phase9b3c_negative_control_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_negative_control_audit.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Negative control audit verifying coordinate permutations, random sets, and leakage checks."
        },
        {
            "file_id": "phase9b3c_evidence_category_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_evidence_category_audit.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Evidence category audit checking classification logic for primary, comparator, and exploratory features."
        },
        {
            "file_id": "phase9b3c_review_findings",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3c_review_findings.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3C",
            "notes": "Review findings table cataloging FIND-01 through FIND-05 with severity, status, and recommended actions."
        }
    ]
    
    rows_added = 0
    rows_updated = 0
    
    for col in df.columns:
        if df[col].dtype == 'float64':
            df[col] = df[col].astype(object)
            
    for f in new_files:
        p = Path(f["local_path"])
        if not p.exists():
            print(f"Warning: File not found: {p}")
            continue
            
        sha = compute_sha256(p)
        size = p.stat().st_size
        
        row_dict = {col: "" for col in df.columns}
        row_dict.update({
            "file_id": f["file_id"],
            "dataset": f["dataset"],
            "sample_id": f["sample_id"],
            "data_type": f["data_type"],
            "local_path": f["local_path"],
            "source_url_or_accession": f["source_url_or_accession"],
            "download_date": f["download_date"],
            "checksum": sha,
            "size_bytes": size,
            "processing_status": f["processing_status"],
            "notes": f["notes"],
            "file_size": float(size),
            "md5": sha
        })
        
        # Check if already exists in df
        match = df[df["local_path"] == f["local_path"]]
        
        if not match.empty:
            idx = match.index[0]
            for col in df.columns:
                df.at[idx, col] = row_dict[col]
            rows_updated += 1
        else:
            df = pd.concat([df, pd.DataFrame([row_dict])], ignore_index=True)
            rows_added += 1
            
    df.fillna("", inplace=True)
    df.to_csv(MANIFEST_PATH, sep="\t", index=False)
    print(f"Manifest updated: {rows_added} files added, {rows_updated} files updated.")

if __name__ == "__main__":
    main()
