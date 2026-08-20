#!/usr/bin/env python3
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
    df = pd.read_csv(MANIFEST_PATH, sep="\t")
    
    # Files to add/update
    new_files = [
        {
            "file_id": "PHASE9B2C2_CORRECTED_SINGLE_CELL_INDEPENDENT_REVIEW",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9B2C2_CORRECTED_SINGLE_CELL_INDEPENDENT_REVIEW.md"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Phase 9B2C2 independent statistical, implementation, annotation, and evidence review of the corrected Phase 9B2R single-cell cellular-source analysis."
        },
        {
            "file_id": "phase9b2c2_correction_verification",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_correction_verification.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Correction verification audit table for Phase 9B2C findings FIND_01, FIND_02, and FIND_03."
        },
        {
            "file_id": "phase9b2c2_provenance_qc_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_provenance_qc_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Provenance and raw data exclusion audit table for PENG_CRA001160."
        },
        {
            "file_id": "phase9b2c2_annotation_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_annotation_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Cell type annotation marker expression audit table."
        },
        {
            "file_id": "phase9b2c2_malignant_cell_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_malignant_cell_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Patient-level malignant cell classification counts."
        },
        {
            "file_id": "phase9b2c2_pseudobulk_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_pseudobulk_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Audit table checking the patient-cell-type pseudobulk inventory eligibility mapping."
        },
        {
            "file_id": "phase9b2c2_feature_eligibility_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_feature_eligibility_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Locked feature coverage and single-cell eligibility audit."
        },
        {
            "file_id": "phase9b2c2_hallmark_tf_results_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_hallmark_tf_results_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Malignant cell basal-classical axis association OLS HC3 results."
        },
        {
            "file_id": "phase9b2c2_composition_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_composition_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Cell fraction composition sensitivity models."
        },
        {
            "file_id": "phase9b2c2_negative_control_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_negative_control_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Audit of negative control execution and statistics."
        },
        {
            "file_id": "phase9b2c2_evidence_category_audit",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_evidence_category_audit.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Final evidence category validation table."
        },
        {
            "file_id": "phase9b2c2_review_findings",
            "dataset": "PENG_CRA001160",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b2c2_review_findings.tsv"),
            "source_url_or_accession": "derived_from_independent_review",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B2C2",
            "notes": "Final table of Phase 9B2C2 review findings, severity, affected features, and correction status."
        }
    ]

    for item in new_files:
        path = Path(item["local_path"])
        if not path.exists():
            print(f"Warning: {path} not found on disk yet (will update hash when generated).")
            # We can create a placeholder size/checksum or skip
            size = 0
            sha = ""
        else:
            size = path.stat().st_size
            sha = compute_sha256(path)
        
        # Check if already in manifest
        idx = df[df["file_id"] == item["file_id"]].index
        if len(idx) > 0:
            df.loc[idx, "file_size"] = size
            df.loc[idx, "md5"] = sha
            df.loc[idx, "local_path"] = item["local_path"]
            print(f"Updated hash for {item['file_id']}: size={size}, sha={sha}")
        else:
            # Append new row
            new_row = {
                "file_id": item["file_id"],
                "dataset": item["dataset"],
                "sample_id": item["sample_id"],
                "data_type": item["data_type"],
                "local_path": item["local_path"],
                "source_url_or_accession": item["source_url_or_accession"],
                "file_size": size,
                "md5": sha,
                "download_date": item["download_date"],
                "processing_status": item["processing_status"],
                "notes": item["notes"]
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
            print(f"Appended new file {item['file_id']} to manifest.")

    # Also update sizes/hashes for administrative logs that might have been modified
    admin_files = [
        {
            "file_id": "09_docs__planning__DECISION_LOG_md",
            "local_path": ROOT / "09_docs/planning/DECISION_LOG.md"
        },
        {
            "file_id": "00_admin__PROJECT_STATUS_md",
            "local_path": ROOT / "00_admin/PROJECT_STATUS.md"
        },
        {
            "file_id": "00_admin__SKILL_USAGE_LOG_tsv",
            "local_path": ROOT / "00_admin/SKILL_USAGE_LOG.tsv"
        }
    ]
    for item in admin_files:
        path = Path(item["local_path"])
        if path.exists():
            size = path.stat().st_size
            sha = compute_sha256(path)
            idx = df[df["file_id"] == item["file_id"]].index
            if len(idx) > 0:
                df.loc[idx, "file_size"] = size
                df.loc[idx, "md5"] = sha
                print(f"Updated hash for admin file {item['file_id']}: size={size}, sha={sha}")

    df.to_csv(MANIFEST_PATH, sep="\t", index=False)
    print("Manifest successfully updated.")

if __name__ == "__main__":
    main()
