#!/usr/bin/env python3
import pandas as pd
import hashlib
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
MANIFEST_PATH = ROOT / "01_metadata/file_manifest.tsv"

files_to_add = [
    {
        "file_id": "04_analysis__09_external_validation__PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW_md",
        "dataset": "PDAC_Phase9B1C2_independent_review",
        "data_type": "phase9b1c2_review_report",
        "local_path": "/Users/emily/thesis/PDAC/04_analysis/09_external_validation/PHASE9B1C2_CORRECTED_BULK_VALIDATION_INDEPENDENT_REVIEW.md",
        "notes": "Independent review report for corrected Phase 9B1R bulk external validation."
    },
    {
        "file_id": "05_results__tables__phase9b1c2_correction_verification_tsv",
        "dataset": "PDAC_Phase9B1C2_independent_review",
        "data_type": "phase9b1c2_review_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9b1c2_correction_verification.tsv",
        "notes": "Phase 9B1C2 correction verification audit table."
    },
    {
        "file_id": "05_results__tables__phase9b1c2_host_feature_audit_tsv",
        "dataset": "PDAC_Phase9B1C2_independent_review",
        "data_type": "phase9b1c2_review_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9b1c2_host_feature_audit.tsv",
        "notes": "Phase 9B1C2 host feature audit table with reviewer evidence classifications."
    },
    {
        "file_id": "05_results__tables__phase9b1c2_module_coverage_audit_tsv",
        "dataset": "PDAC_Phase9B1C2_independent_review",
        "data_type": "phase9b1c2_review_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9b1c2_module_coverage_audit.tsv",
        "notes": "Phase 9B1C2 module coverage audit table enforcing 80% threshold."
    },
    {
        "file_id": "05_results__tables__phase9b1c2_negative_control_audit_tsv",
        "dataset": "PDAC_Phase9B1C2_independent_review",
        "data_type": "phase9b1c2_review_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9b1c2_negative_control_audit.tsv",
        "notes": "Phase 9B1C2 negative control audit table with iterations and seed."
    },
    {
        "file_id": "05_results__tables__phase9b1c2_review_findings_tsv",
        "dataset": "PDAC_Phase9B1C2_independent_review",
        "data_type": "phase9b1c2_review_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9b1c2_review_findings.tsv",
        "notes": "Phase 9B1C2 findings table with verification status for all implementation issues."
    }
]

def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()

def main():
    # Load manifest
    df = pd.read_csv(MANIFEST_PATH, sep="\t")
    
    # Check if files already present, if so remove them so we can write fresh
    file_ids_to_add = [f["file_id"] for f in files_to_add]
    df = df[~df["file_id"].isin(file_ids_to_add)]
    
    rows = []
    for f in files_to_add:
        path = Path(f["local_path"])
        if not path.exists():
            print(f"Error: {path} does not exist")
            continue
            
        size = path.stat().st_size
        sha = compute_sha256(path)
        
        row = {
            "file_id": f["file_id"],
            "dataset": f["dataset"],
            "sample_id": pd.NA,
            "data_type": f["data_type"],
            "local_path": f["local_path"],
            "source_url_or_accession": "generated_Phase9B1C2_from_corrected_Phase9B1R_results",
            "file_size": size,
            "md5": sha, # Note: column is named md5 but contains sha256 as per other records
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B1C2",
            "notes": f["notes"]
        }
        rows.append(row)
        
    df_new = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
    
    # Save manifest
    df_new.to_csv(MANIFEST_PATH, sep="\t", index=False)
    print(f"Successfully appended {len(rows)} files to manifest")

if __name__ == "__main__":
    main()
