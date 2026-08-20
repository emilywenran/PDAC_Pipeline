#!/usr/bin/env python3
"""Register Phase 9B3R0 files in the file manifest."""

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
            "file_id": "PHASE9B3R0_DUAL_MODEL_REANALYSIS_AUDIT",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9B3R0_DUAL_MODEL_REANALYSIS_AUDIT.md"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3R0",
            "notes": "Dual-model pre-reanalysis audit report outlining implementation failures."
        },
        {
            "file_id": "phase9b3r0_reconciled_repair_specification",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "result_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3r0_reconciled_repair_specification.tsv"),
            "source_url_or_accession": "derived_from_audit",
            "download_date": "2026-07-04",
            "processing_status": "generated_Phase9B3R0",
            "notes": "TSV mapping specific findings to code repairs, statistical repairs, and unit tests."
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
