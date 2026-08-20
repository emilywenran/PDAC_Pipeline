import os
import hashlib

manifest_file = "/Users/emily/thesis/PDAC/01_metadata/file_manifest.tsv"

new_files = [
    {
        "file_id": "external_validation_dataset_inventory",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "metadata_inventory",
        "local_path": "/Users/emily/thesis/PDAC/01_metadata/external_validation_dataset_inventory.tsv",
        "source_url_or_accession": "derived_from_public_literature_and_database_search",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Inventory of 12 candidate external bulk, single-cell, spatial, and microbiome datasets."
    },
    {
        "file_id": "external_validation_parameter_inventory",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "parameter_inventory",
        "local_path": "/Users/emily/thesis/PDAC/01_metadata/external_validation_parameter_inventory.tsv",
        "source_url_or_accession": "derived_from_locked_validation_framework",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Locked parameter specifications for 10 bulk, single-cell, spatial, and microbiome validation runs."
    },
    {
        "file_id": "phase9a_external_dataset_shortlist",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "results_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9a_external_dataset_shortlist.tsv",
        "source_url_or_accession": "derived_from_dataset_inventory",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Shortlist of 7 PRIORITY_1 datasets for Phase 9B validation."
    },
    {
        "file_id": "phase9a_signature_external_coverage_feasibility",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "results_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9a_signature_external_coverage_feasibility.tsv",
        "source_url_or_accession": "derived_from_expression_manifest_and_external_platforms",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Coverage feasibility checks for Moffitt50, Hallmark pathways, and WGCNA module signatures in external cohorts."
    },
    {
        "file_id": "phase9a_external_analysis_resource_estimate",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "results_table",
        "local_path": "/Users/emily/thesis/PDAC/05_results/tables/phase9a_external_analysis_resource_estimate.tsv",
        "source_url_or_accession": "derived_from_shortlist_sizes_and_requirements",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Estimated download sizes, software, memory, and runtime resource estimates for Phase 9B execution."
    },
    {
        "file_id": "phase9_external_validation_sources",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "reference_bibliography",
        "local_path": "/Users/emily/thesis/PDAC/09_docs/references/phase9_external_validation_sources.bib",
        "source_url_or_accession": "derived_from_citation_management",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "BibTeX bibliography of 12 external validation papers and datasets."
    },
    {
        "file_id": "phase9_external_validation_source_audit",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "reference_audit",
        "local_path": "/Users/emily/thesis/PDAC/09_docs/references/phase9_external_validation_source_audit.tsv",
        "source_url_or_accession": "derived_from_citation_management",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Audit table containing title, authors, year, journal, DOI, PMID, and accession verification for all sources."
    },
    {
        "file_id": "PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "analysis_report",
        "local_path": "/Users/emily/thesis/PDAC/04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md",
        "source_url_or_accession": "derived_from_planning",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Prospective method lock document for Phase 9B external validation."
    },
    {
        "file_id": "PDAC_external_validation_protocol",
        "dataset": "PDAC_Phase9_validation",
        "sample_id": "",
        "data_type": "method_protocol",
        "local_path": "/Users/emily/thesis/PDAC/09_docs/methods/PDAC_external_validation_protocol.md",
        "source_url_or_accession": "derived_from_planning",
        "download_date": "2026-07-02",
        "processing_status": "generated_Phase9A",
        "notes": "Standard operating protocol and step-by-step replication guide for Phase 9B execution."
    }
]

def get_file_stats(path):
    size = os.path.getsize(path)
    sha256_hash = hashlib.sha256()
    with open(path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    checksum = f"sha256:{sha256_hash.hexdigest()}"
    return size, checksum

# Read existing manifest
with open(manifest_file, "r") as f:
    lines = f.read().splitlines()

header = lines[0]
existing_ids = set()
for line in lines[1:]:
    if not line.strip():
        continue
    parts = line.split("\t")
    existing_ids.add(parts[0])

# Append new files
new_lines = []
for f_info in new_files:
    if f_info["file_id"] in existing_ids:
        print(f"Skipping existing file in manifest: {f_info['file_id']}")
        continue
    size, checksum = get_file_stats(f_info["local_path"])
    row = "\t".join([
        f_info["file_id"],
        f_info["dataset"],
        f_info["sample_id"],
        f_info["data_type"],
        f_info["local_path"],
        f_info["source_url_or_accession"],
        str(size),
        checksum,
        f_info["download_date"],
        f_info["processing_status"],
        f_info["notes"]
    ])
    new_lines.append(row)

if new_lines:
    with open(manifest_file, "a") as f:
        f.write("\n".join(new_lines) + "\n")
    print(f"Appended {len(new_lines)} files to file_manifest.tsv")
else:
    print("No new files to append")
