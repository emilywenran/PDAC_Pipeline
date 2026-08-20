#!/usr/bin/env python3
"""Register Phase 9B3A Spatial Validation Planning files in the file manifest."""

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
    
    # New and updated spatial validation planning files
    new_files = [
        {
            "file_id": "PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9B3A_SPATIAL_VALIDATION_METHOD_LOCK.md"),
            "source_url_or_accession": "derived_from_planning",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Prospective method lock document for Phase 9B3 spatial validation."
        },
        {
            "file_id": "PDAC_spatial_validation_protocol",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "method_protocol",
            "local_path": str(ROOT / "09_docs/methods/PDAC_spatial_validation_protocol.md"),
            "source_url_or_accession": "derived_from_planning",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Execution protocol for spatial validation mapping and statistical models."
        },
        {
            "file_id": "phase9b3_spatial_dataset_inventory",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "metadata_inventory",
            "local_path": str(ROOT / "01_metadata/phase9b3_spatial_dataset_inventory.tsv"),
            "source_url_or_accession": "derived_from_literature_audit",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Detailed inventory of three qualified/split candidate spatial cohorts."
        },
        {
            "file_id": "phase9b3_spatial_parameter_inventory",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "parameter_inventory",
            "local_path": str(ROOT / "01_metadata/phase9b3_spatial_parameter_inventory.tsv"),
            "source_url_or_accession": "derived_from_statistical_lock",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Prospective analysis parameters for spatial deconvolution, scoring, and models."
        },
        {
            "file_id": "phase9b3a_spatial_feature_hierarchy",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3a_spatial_feature_hierarchy.tsv"),
            "source_url_or_accession": "derived_from_single_cell_reanalysis",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Locked feature hierarchy table defining 33 spatial features, roles, and coverage criteria."
        },
        {
            "file_id": "phase9b3a_spatial_dataset_qualification",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3a_spatial_dataset_qualification.tsv"),
            "source_url_or_accession": "derived_from_dataset_inventory",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Qualification, rationale, and exclusion rules applied to 5 candidate spatial datasets."
        },
        {
            "file_id": "phase9b3a_authoritative_spatial_cohort_set",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3a_authoritative_spatial_cohort_set.tsv"),
            "source_url_or_accession": "derived_from_dataset_qualification",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Authoritative cohort set detailing primary, secondary, and sensitivity analysis roles."
        },
        {
            "file_id": "phase9b3a_spatial_resource_estimate",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3a_spatial_resource_estimate.tsv"),
            "source_url_or_accession": "derived_from_cohort_sizes",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Runtime, RAM, storage, and software resource estimates for Apple Silicon and HPC."
        },
        {
            "file_id": "validate_phase9b3a_spatial_plan",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "validation_script",
            "local_path": str(ROOT / "06_scripts/python/16_validate_phase9b3a_spatial_plan.py"),
            "source_url_or_accession": "derived_from_validation_framework",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Python validator enforcing scientific and technical rules for the spatial validation plan."
        },
        {
            "file_id": "generate_phase9b3a_tables",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "helper_script",
            "local_path": str(ROOT / "06_scripts/python/generate_phase9b3a_tables.py"),
            "source_url_or_accession": "derived_from_helper_tools",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A",
            "notes": "Helper script to generate TSV planning tables programmatically."
        },
        # Phase 9B3A.1 specific additions
        {
            "file_id": "PHASE9B3A1_SPATIAL_DESIGN_CONSISTENCY_CORRECTION",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9B3A1_SPATIAL_DESIGN_CONSISTENCY_CORRECTION.md"),
            "source_url_or_accession": "derived_from_planning",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A1",
            "notes": "Design consistency correction document for Phase 9B3 spatial validation."
        },
        {
            "file_id": "phase9b3a1_spatial_analysis_unit_and_models",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3a1_spatial_analysis_unit_and_models.tsv"),
            "source_url_or_accession": "derived_from_statistical_lock",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A1",
            "notes": "Table mapping analysis units, nested hierarchies, LMMs, and reduced-model fallback rules."
        },
        {
            "file_id": "generate_phase9b3a1_tables",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "helper_script",
            "local_path": str(ROOT / "06_scripts/python/generate_phase9b3a1_tables.py"),
            "source_url_or_accession": "derived_from_helper_tools",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A1",
            "notes": "Helper script to generate the models and analysis units TSV table."
        },
        # Phase 9B3A.2 specific additions
        {
            "file_id": "PHASE9B3A2_SPATIAL_HIERARCHY_FINAL_CORRECTION",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "analysis_report",
            "local_path": str(ROOT / "04_analysis/09_external_validation/PHASE9B3A2_SPATIAL_HIERARCHY_FINAL_CORRECTION.md"),
            "source_url_or_accession": "derived_from_planning",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A2",
            "notes": "Final correction report locking paired models and exploratory status."
        },
        {
            "file_id": "phase9b3a2_spatial_model_hierarchy",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "results_table",
            "local_path": str(ROOT / "05_results/tables/phase9b3a2_spatial_model_hierarchy.tsv"),
            "source_url_or_accession": "derived_from_statistical_lock",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A2",
            "notes": "Table locking the specific Model A, B, C hierarchy for GeoMx and exploratory ST protocol."
        },
        {
            "file_id": "generate_phase9b3a2_tables",
            "dataset": "PDAC_Phase9_validation",
            "sample_id": "",
            "data_type": "helper_script",
            "local_path": str(ROOT / "06_scripts/python/generate_phase9b3a2_tables.py"),
            "source_url_or_accession": "derived_from_helper_tools",
            "download_date": "2026-07-03",
            "processing_status": "generated_Phase9B3A2",
            "notes": "Helper script to generate Model A, B, C hierarchy and exploratory ST parameters."
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
            "notes": f["notes"]
        })
        
        # Check if already exists in df
        match = df[df["local_path"] == f["local_path"]]
        
        if not match.empty:
            # Update
            idx = match.index[0]
            for col in df.columns:
                df.at[idx, col] = row_dict[col]
            rows_updated += 1
        else:
            # Append
            df = pd.concat([df, pd.DataFrame([row_dict])], ignore_index=True)
            rows_added += 1
            
    df.fillna("", inplace=True)
    df.to_csv(MANIFEST_PATH, sep="\t", index=False)
    print(f"Manifest updated: {rows_added} files added, {rows_updated} files updated.")

if __name__ == "__main__":
    main()
