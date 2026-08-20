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
    
    # 1. Update existing files in the manifest
    existing_updates = [
        {
            "file_id": "external_validation_dataset_inventory",
            "local_path": ROOT / "01_metadata/external_validation_dataset_inventory.tsv"
        },
        {
            "file_id": "external_validation_parameter_inventory",
            "local_path": ROOT / "01_metadata/external_validation_parameter_inventory.tsv"
        },
        {
            "file_id": "phase9a_external_dataset_shortlist",
            "local_path": ROOT / "05_results/tables/phase9a_external_dataset_shortlist.tsv"
        },
        {
            "file_id": "phase9_external_validation_sources",
            "local_path": ROOT / "09_docs/references/phase9_external_validation_sources.bib"
        },
        {
            "file_id": "phase9_external_validation_source_audit",
            "local_path": ROOT / "09_docs/references/phase9_external_validation_source_audit.tsv"
        },
        {
            "file_id": "PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK",
            "local_path": ROOT / "04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md"
        },
        {
            "file_id": "PDAC_external_validation_protocol",
            "local_path": ROOT / "09_docs/methods/PDAC_external_validation_protocol.md"
        },
        {
            "file_id": "PHASE9A1_SINGLE_CELL_COHORT_RECONCILIATION",
            "local_path": ROOT / "04_analysis/09_external_validation/PHASE9A1_SINGLE_CELL_COHORT_RECONCILIATION.md"
        },
        {
            "file_id": "PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS",
            "local_path": ROOT / "04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md"
        },
        {
            "file_id": "phase9b2_single_cell_dataset_inventory",
            "local_path": ROOT / "01_metadata/phase9b2_single_cell_dataset_inventory.tsv"
        },
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

    for item in existing_updates:
        path = Path(item["local_path"])
        if not path.exists():
            print(f"Warning: {path} not found")
            continue
        size = path.stat().st_size
        sha = compute_sha256(path)
        
        idx = df[df["file_id"] == item["file_id"]].index
        if len(idx) > 0:
            df.loc[idx, "file_size"] = size
            df.loc[idx, "md5"] = sha
            print(f"Updated hash for {item['file_id']}: size={size}, sha={sha}")
        else:
            # If not in manifest but is an existing file, we can also add it
            print(f"File {item['file_id']} not found in manifest, will append later if new.")
            
    # 2. Add new Phase 9A.2 files if not present
    new_files = [
        {
            "file_id": "PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION.md"),
            "source_url_or_accession": "derived_from_planning",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9A2",
            "notes": "Provenance correction and expanded cohort set report for Phase 9B2 single-cell validation."
        },
        {
            "file_id": "phase9a2_single_cell_dataset_provenance_audit",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9a2_single_cell_dataset_provenance_audit.tsv"),
            "source_url_or_accession": "derived_from_dataset_inventory",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9A2",
            "notes": "Audited single-cell and spatial dataset inventory table."
        },
        {
            "file_id": "phase9a2_phase9b2_authoritative_cohort_set",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv"),
            "source_url_or_accession": "derived_from_shortlist_and_method_lock",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9A2",
            "notes": "Authoritative single-cell/spatial execution set for Phase 9B2."
        },
        {
            "file_id": "generate_phase9a2_tables",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_script",
            "local_path": str(ROOT / "06_scripts/python/generate_phase9a2_tables.py"),
            "source_url_or_accession": "generated_manually",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9A2",
            "notes": "Python script to generate corrected Phase 9A.2 metadata tables."
        },
        {
            "file_id": "15_validate_provenance_consistency",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "validation_script",
            "local_path": str(ROOT / "06_scripts/python/15_validate_provenance_consistency.py"),
            "source_url_or_accession": "generated_manually",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9A2",
            "notes": "Python script to validate provenance consistency across all planning files."
        }
    ]
    
    new_rows = []
    for f_info in new_files:
        if f_info["file_id"] in df["file_id"].values:
            # Already present, let's update it instead
            path = Path(f_info["local_path"])
            size = path.stat().st_size
            sha = compute_sha256(path)
            idx = df[df["file_id"] == f_info["file_id"]].index
            df.loc[idx, "file_size"] = size
            df.loc[idx, "md5"] = sha
            print(f"Updated hash for new file {f_info['file_id']}: size={size}, sha={sha}")
        else:
            path = Path(f_info["local_path"])
            size = path.stat().st_size
            sha = compute_sha256(path)
            new_rows.append({
                "file_id": f_info["file_id"],
                "dataset": f_info["dataset"],
                "sample_id": f_info["sample_id"],
                "data_type": f_info["data_type"],
                "local_path": f_info["local_path"],
                "source_url_or_accession": f_info["source_url_or_accession"],
                "file_size": size,
                "md5": sha,
                "download_date": f_info["download_date"],
                "processing_status": f_info["processing_status"],
                "notes": f_info["notes"]
            })
            print(f"Added {f_info['file_id']} to manifest.")
            
    if new_rows:
        df_new = pd.DataFrame(new_rows)
        df = pd.concat([df, df_new], ignore_index=True)
        
    df.to_csv(MANIFEST_PATH, sep="\t", index=False)
    print("Successfully updated and saved file_manifest.tsv")

if __name__ == "__main__":
    main()
