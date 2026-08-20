#!/usr/bin/env python3
"""Provenance and Cohort Set Consistency Validator.

This script validates the metadata, provenance, and cohort definition consistency
across PDAC external validation planning files for Phase 9A.2.
It fails (exits with code 1) if any rules are violated.
"""

from __future__ import annotations
import csv
import sys
from pathlib import Path

ROOT = Path("/Users/emily/thesis/PDAC")
INVENTORY_PATH = ROOT / "01_metadata/external_validation_dataset_inventory.tsv"
SHORTLIST_PATH = ROOT / "05_results/tables/phase9a_external_dataset_shortlist.tsv"
PARAMETER_PATH = ROOT / "01_metadata/external_validation_parameter_inventory.tsv"
AUTHORITATIVE_PATH = ROOT / "05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv"
AUDIT_PATH = ROOT / "05_results/tables/phase9a2_single_cell_dataset_provenance_audit.tsv"

def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        print(f"Error: Required file not found: {path}", file=sys.stderr)
        sys.exit(1)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def main() -> int:
    errors: list[str] = []

    # Read all planning files
    inventory = read_tsv(INVENTORY_PATH)
    shortlist = read_tsv(SHORTLIST_PATH)
    parameters = read_tsv(PARAMETER_PATH)
    authoritative = read_tsv(AUTHORITATIVE_PATH)
    audit = read_tsv(AUDIT_PATH)

    print("Checking database-lookup, citation-management, and experimental-design rules...")

    # 1. Check accession to publication mapping and label errors
    for row in inventory:
        acc = row.get("accession", "")
        pub = row.get("publication", "")
        ds_id = row.get("canonical_dataset_id", "")

        # GSE111672 must not be Peng 2019
        if acc == "GSE111672" and "peng" in pub.lower():
            errors.append(f"Invalid mapping in inventory: GSE111672 is assigned to Peng: {row}")
        # CRA001160 must not be Moncada
        if acc == "CRA001160" and "moncada" in pub.lower():
            errors.append(f"Invalid mapping in inventory: CRA001160 is assigned to Moncada: {row}")

        # CRA001160 must be Peng 2019
        if acc == "CRA001160" and "peng" not in pub.lower():
            errors.append(f"Expected CRA001160 to map to Peng publication, got: {pub}")
        # GSE111672 must be Moncada 2020 (when representing Moncada)
        if ds_id == "MONCADA_GSE111672" and "moncada" not in pub.lower():
            errors.append(f"Expected MONCADA_GSE111672 to map to Moncada publication, got: {pub}")
        # GSE154778 must map to Lin et al. (2020)
        if acc == "GSE154778" and "lin" not in pub.lower():
            errors.append(f"Expected GSE154778 to map to Lin publication, got: {pub}")
        # GSE202051 must map to Hwang or Lin
        if acc == "GSE202051" and "hwang" not in pub.lower() and "lin" not in pub.lower():
            errors.append(f"Expected GSE202051 to map to Hwang/Lin publication, got: {pub}")

    # 2. Check that one accession does not represent two different studies in inventory
    acc_to_ids: dict[str, set[str]] = {}
    for row in inventory:
        acc = row.get("accession", "")
        ds_id = row.get("canonical_dataset_id", "")
        if acc and ds_id:
            if acc not in acc_to_ids:
                acc_to_ids[acc] = set()
            # Wait, TCGA_PAAD and TCGA_PAAD_microbiome might share TCGA-PAAD, which is fine since it's the same cohort
            # GSE111672 and MONCADA_GSE111672 might be the same, but Peng must not share GSE111672.
            # We want to check if the accession belongs to two completely different biological studies
            acc_to_ids[acc].add(ds_id)
    
    for acc, ids in acc_to_ids.items():
        # Check if the IDs represent different authors/studies
        # We know PENG_CRA001160 and MONCADA_GSE111672 must have different accessions
        if "PENG_CRA001160" in ids and "MONCADA_GSE111672" in ids:
            errors.append(f"Accession {acc} represents both Peng and Moncada datasets: {ids}")

    # 3. Patient counts conflict across files
    # We will build a patient count mapping for the four single-cell/spatial datasets:
    # Peng = 24, Moncada = 6, Lin = 10, Hwang = 43
    expected_patients = {
        "PENG_CRA001160": 24,
        "MONCADA_GSE111672": 6,
        "LIN_GSE154778": 10,  # 10 primary tumors
        "HWANG_GSE202051": 43,
    }

    # Verify inventory counts
    for row in inventory:
        ds_id = row.get("canonical_dataset_id", "")
        if ds_id in expected_patients:
            val = int(row.get("PDAC_patients", 0))
            expected = expected_patients[ds_id]
            if val != expected:
                errors.append(f"Patient count conflict in inventory for {ds_id}: expected {expected}, got {val}")

    # Verify shortlist counts
    for row in shortlist:
        ds_id = row.get("dataset_id", "")
        if ds_id in expected_patients:
            val = int(row.get("sample_count", 0))
            expected = expected_patients[ds_id]
            if val != expected:
                errors.append(f"Patient count conflict in shortlist for {ds_id}: expected {expected}, got {val}")

    # Verify authoritative set counts
    for row in authoritative:
        ds_id = row.get("canonical_dataset_id", "")
        if ds_id in expected_patients:
            val = int(row.get("PDAC_patients", 0))
            expected = expected_patients[ds_id]
            if val != expected:
                errors.append(f"Patient count conflict in authoritative cohort set for {ds_id}: expected {expected}, got {val}")

    # Verify audit set counts
    for row in audit:
        ds_id = row.get("canonical_dataset_id", "")
        if ds_id in expected_patients:
            val = int(row.get("PDAC_patients", 0))
            expected = expected_patients[ds_id]
            if val != expected:
                errors.append(f"Patient count conflict in single-cell audit table for {ds_id}: expected {expected}, got {val}")

    # 4. Check that GSE111672 is not labelled as Peng, and CRA001160 is not Moncada
    for row in inventory:
        acc = row.get("accession", "")
        ds_id = row.get("canonical_dataset_id", "")
        if acc == "GSE111672" and ds_id == "PENG_CRA001160":
            errors.append("GSE111672 is labelled as PENG_CRA001160 cohort.")
        if acc == "CRA001160" and ds_id == "MONCADA_GSE111672":
            errors.append("CRA001160 is labelled as MONCADA_GSE111672 cohort.")

    # 5. A dataset is included in Phase 9B2 without matching parameter-inventory rows
    # The included datasets in Phase 9B2 are: PENG_CRA001160, LIN_GSE154778, MONCADA_GSE111672, HWANG_GSE202051
    # Plus bulk cohorts (TCGA_PAAD, GSE71729, GSE62452) and microbiome/spatial cohorts.
    # We will check that every dataset ID in the authoritative cohort table has at least one row in parameter inventory.
    auth_ids = {
        row.get("canonical_dataset_id", "")
        for row in authoritative
        if row.get("included_in_phase9b2_primary") == "TRUE"
        or row.get("included_in_phase9b2_supplementary") == "TRUE"
    }
    param_ids = {row.get("dataset_id", "") for row in parameters}
    
    missing_params = auth_ids - param_ids
    if missing_params:
        errors.append(f"Authoritative dataset(s) included in Phase 9B2 lack parameter-inventory rows: {missing_params}")


    # 6. The authoritative cohort set differs across planning files
    # The set of Phase 9B2 single-cell and spatial cohorts must be exactly:
    # {"PENG_CRA001160", "LIN_GSE154778", "MONCADA_GSE111672", "HWANG_GSE202051"}
    expected_sc_spatial_set = {"PENG_CRA001160", "LIN_GSE154778", "MONCADA_GSE111672", "HWANG_GSE202051"}
    
    auth_sc_spatial = {row.get("canonical_dataset_id", "") for row in authoritative}
    if auth_sc_spatial != expected_sc_spatial_set:
        errors.append(f"Authoritative cohort set differs from expected: got {auth_sc_spatial}, expected {expected_sc_spatial_set}")

    audit_sc_spatial = {row.get("canonical_dataset_id", "") for row in audit}
    if audit_sc_spatial != expected_sc_spatial_set:
        errors.append(f"Audit cohort set differs from expected: got {audit_sc_spatial}, expected {expected_sc_spatial_set}")

    # Report results
    if errors:
        print("Provenance consistency validation FAILED:", file=sys.stderr)
        for err in errors:
            print(f"- {err}", file=sys.stderr)
        return 1

    print("Provenance consistency validation PASSED successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
