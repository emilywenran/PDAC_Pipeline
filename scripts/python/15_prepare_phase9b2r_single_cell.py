#!/usr/bin/env python3
"""Prepare and verify immutable inputs for Phase 9B2R.

This script does not download data. It verifies the PENG_CRA001160-only scope,
processed-file provenance, core Phase 9B2 artifacts that Phase 9B2C found
unaffected, and writes a Phase 9B2R runtime verification table.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results/tables"
DATA_DIR = ROOT / "02_data/external/phase9_single_cell/PENG_CRA001160"
MODEL_DIR = ROOT / "05_results/models/phase9b2"

EXPECTED = {
    "canonical_dataset_id": "PENG_CRA001160",
    "accession": "CRA001160",
    "BioProject": "PRJCA001063",
    "cells": 57530,
    "patients": 35,
    "tumor_patients": 24,
    "control_patients": 11,
    "min_cells_per_patient_celltype": 20,
    "random_seed": 2026,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def add(rows: list[dict[str, object]], item: str, expected: object, observed: object, status: bool, notes: str) -> None:
    rows.append({
        "validation_item": item,
        "expected_value": expected,
        "observed_value": observed,
        "status": "PASS" if status else "FAIL",
        "notes": notes,
    })


def main() -> int:
    rows: list[dict[str, object]] = []

    auth = read_tsv(TABLE_DIR / "phase9a2_phase9b2_authoritative_cohort_set.tsv")
    primary = sorted(r["canonical_dataset_id"] for r in auth if r.get("included_in_phase9b2_primary", "").upper() == "TRUE")
    add(rows, "authoritative_primary_scope", "PENG_CRA001160", ";".join(primary), primary == ["PENG_CRA001160"], "Preserve Phase 9A.3 PENG-only primary scope.")

    inv = read_tsv(ROOT / "01_metadata/external_validation_dataset_inventory.tsv")
    peng = [r for r in inv if r.get("canonical_dataset_id") == "PENG_CRA001160"]
    p = peng[0] if peng else {}
    for field in ("accession", "BioProject"):
        add(rows, f"provenance_{field}", EXPECTED[field], p.get(field, "MISSING"), p.get(field) == EXPECTED[field], "Official Peng cohort identity.")
    add(rows, "provenance_publication", "Peng et al. 2019", p.get("publication", ""), "peng" in p.get("publication", "").lower() and "2019" in p.get("publication", ""), "Publication identity.")
    add(rows, "no_gse111672_alias", "not GSE111672", p.get("accession", ""), p.get("accession") != "GSE111672", "GSE111672 remains Moncada, not Peng.")

    raw_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram", ".sra")
    raw = [str(path.relative_to(ROOT)) for path in DATA_DIR.rglob("*") if path.is_file() and path.name.lower().endswith(raw_suffixes)]
    add(rows, "raw_file_absence", "no FASTQ/BAM/SRA/CRAM", ";".join(raw) or "none", not raw, "Phase 9B2R must use processed files only.")

    for name in ("count-matrix.txt", "all_celltype.txt", "md5sum.txt"):
        path = DATA_DIR / name
        add(rows, f"processed_file_present_{name}", "present", "present" if path.exists() else "missing", path.exists(), "Official processed file required.")
        if path.exists():
            add(rows, f"processed_file_sha256_{name}", "recorded in manifest", sha256(path), True, "Checksum recalculated for audit trail.")

    summary_path = MODEL_DIR / "phase9b2_prepare_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    for key in ("cells", "patients", "tumor_patients", "control_patients", "random_seed", "min_cells_per_patient_celltype"):
        add(rows, f"prepare_summary_{key}", EXPECTED[key], summary.get(key, "MISSING"), summary.get(key) == EXPECTED[key], "Unchanged core prepare summary verification.")

    cohort_qc = read_tsv(TABLE_DIR / "phase9b2_single_cell_cohort_qc.tsv")
    cq = cohort_qc[0] if cohort_qc else {}
    add(rows, "cohort_qc_cell_count", EXPECTED["cells"], cq.get("cells", ""), str(EXPECTED["cells"]) == cq.get("cells"), "Cell count preserved.")
    add(rows, "cohort_qc_patient_count", EXPECTED["patients"], cq.get("patients", ""), str(EXPECTED["patients"]) == cq.get("patients"), "Patient/donor count preserved.")
    add(rows, "patient_independence", "patient-level pseudobulk", cq.get("patient_independence_status", ""), "patient-level pseudobulk" in cq.get("patient_independence_status", ""), "No cell-level pseudoreplication.")

    pb = read_tsv(TABLE_DIR / "phase9b2_pseudobulk_inventory.tsv")
    disagreements = [
        r for r in pb
        if (int(float(r["number_of_cells"])) >= EXPECTED["min_cells_per_patient_celltype"]) != (r["eligibility"] == "ELIGIBLE")
    ]
    add(rows, "pseudobulk_eligibility_rule", ">=20 cells -> ELIGIBLE", len(disagreements), len(disagreements) == 0, "Phase 9B2C pseudobulk audit rule reproduced.")

    module_cov = read_tsv(TABLE_DIR / "phase9b2_module_transfer_coverage.tsv")
    low_modules = [r["module_name"] for r in module_cov if float(r["coverage_fraction"]) < 0.80]
    add(rows, "phase9b2c_module_coverage_reproduced", "all five modules below 0.80", ";".join(low_modules), set(low_modules) == {"MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"}, "Stop if this disagrees with Phase 9B2C.")

    write_tsv(TABLE_DIR / "phase9b2r_runtime_validation.tsv", rows)
    failed = [r for r in rows if r["status"] == "FAIL"]
    if failed:
        print("Phase 9B2R preparation validation FAILED:", file=sys.stderr)
        for row in failed:
            print(f"- {row['validation_item']}: observed {row['observed_value']}", file=sys.stderr)
        return 1
    print("Phase 9B2R preparation validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
