#!/usr/bin/env python3
"""Validate Phase 0 metadata manifests without modifying them."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_manifest": [
        "sample_id",
        "geo_sample_id",
        "biosample_id",
        "sra_run_id",
        "patient_id",
        "expression_available",
        "microbiome_available",
        "subtype_original",
        "subtype_source",
        "survival_available",
        "batch",
        "tumor_purity",
        "notes",
    ],
    "clinical_metadata": [
        "patient_id",
        "age",
        "sex",
        "stage",
        "grade",
        "treatment",
        "overall_survival_time",
        "overall_survival_event",
        "disease_free_survival_time",
        "disease_free_survival_event",
        "source",
        "missingness_note",
    ],
    "file_manifest": [
        "file_id",
        "dataset",
        "sample_id",
        "data_type",
        "local_path",
        "source_url_or_accession",
        "file_size",
        "md5",
        "download_date",
        "processing_status",
        "notes",
    ],
}

PRIMARY_IDS = {
    "sample_manifest": "sample_id",
    "clinical_metadata": "patient_id",
    "file_manifest": "file_id",
}

SECONDARY_IDS = {
    "sample_manifest": ["geo_sample_id", "biosample_id", "sra_run_id"],
    "clinical_metadata": [],
    "file_manifest": [],
}


def read_tsv(path: Path, manifest_name: str) -> tuple[list[str], list[dict[str, str]], list[str]]:
    errors: list[str] = []
    if not path.exists():
        return [], [], [f"{manifest_name}: file not found: {path}"]

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)

    missing = [col for col in REQUIRED_COLUMNS[manifest_name] if col not in header]
    if missing:
        errors.append(f"{manifest_name}: missing required columns: {', '.join(missing)}")

    return header, rows, errors


def nonempty(value: str | None) -> bool:
    return bool((value or "").strip())


def validate_identifiers(manifest_name: str, rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    primary = PRIMARY_IDS[manifest_name]

    empty_rows = [idx for idx, row in enumerate(rows, start=2) if not nonempty(row.get(primary))]
    if empty_rows:
        errors.append(
            f"{manifest_name}: empty required identifier '{primary}' at rows: "
            + ", ".join(map(str, empty_rows))
        )

    primary_counts = Counter((row.get(primary) or "").strip() for row in rows if nonempty(row.get(primary)))
    duplicates = sorted(identifier for identifier, count in primary_counts.items() if count > 1)
    if duplicates:
        errors.append(f"{manifest_name}: duplicate {primary} values: {', '.join(duplicates)}")

    for column in SECONDARY_IDS[manifest_name]:
        counts = Counter((row.get(column) or "").strip() for row in rows if nonempty(row.get(column)))
        duplicate_values = sorted(identifier for identifier, count in counts.items() if count > 1)
        if duplicate_values:
            errors.append(f"{manifest_name}: duplicate {column} values: {', '.join(duplicate_values)}")

    return errors


def validate_cross_manifest_mappings(
    sample_rows: list[dict[str, str]],
    clinical_rows: list[dict[str, str]],
    file_rows: list[dict[str, str]],
) -> list[str]:
    errors: list[str] = []

    sample_ids = {(row.get("sample_id") or "").strip() for row in sample_rows if nonempty(row.get("sample_id"))}
    patient_ids = {
        (row.get("patient_id") or "").strip() for row in clinical_rows if nonempty(row.get("patient_id"))
    }

    sample_patient_map: dict[str, set[str]] = defaultdict(set)
    for row in sample_rows:
        sample_id = (row.get("sample_id") or "").strip()
        patient_id = (row.get("patient_id") or "").strip()
        if sample_id and patient_id:
            sample_patient_map[sample_id].add(patient_id)

    inconsistent_samples = sorted(
        sample_id for sample_id, mapped_patients in sample_patient_map.items() if len(mapped_patients) > 1
    )
    if inconsistent_samples:
        errors.append(
            "sample_manifest: sample_id maps to multiple patient_id values: "
            + ", ".join(inconsistent_samples)
        )

    missing_clinical = sorted(
        patient_id
        for mapped_patients in sample_patient_map.values()
        for patient_id in mapped_patients
        if patient_ids and patient_id not in patient_ids
    )
    if missing_clinical:
        errors.append(
            "sample_manifest: patient_id values absent from clinical_metadata: "
            + ", ".join(sorted(set(missing_clinical)))
        )

    file_sample_ids = sorted(
        {
            (row.get("sample_id") or "").strip()
            for row in file_rows
            if nonempty(row.get("sample_id"))
        }
    )
    missing_samples = [sample_id for sample_id in file_sample_ids if sample_ids and sample_id not in sample_ids]
    if missing_samples:
        errors.append(
            "file_manifest: sample_id values absent from sample_manifest: "
            + ", ".join(sorted(set(missing_samples)))
        )

    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate PDAC Phase 0 metadata manifests.")
    parser.add_argument("--metadata-dir", default="01_metadata", help="Directory containing manifest TSV files.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metadata_dir = Path(args.metadata_dir)
    paths = {
        "sample_manifest": metadata_dir / "sample_manifest.tsv",
        "clinical_metadata": metadata_dir / "clinical_metadata.tsv",
        "file_manifest": metadata_dir / "file_manifest.tsv",
    }

    all_errors: list[str] = []
    rows_by_manifest: dict[str, list[dict[str, str]]] = {}

    for manifest_name, path in paths.items():
        _, rows, errors = read_tsv(path, manifest_name)
        rows_by_manifest[manifest_name] = rows
        all_errors.extend(errors)
        all_errors.extend(validate_identifiers(manifest_name, rows))

    all_errors.extend(
        validate_cross_manifest_mappings(
            rows_by_manifest.get("sample_manifest", []),
            rows_by_manifest.get("clinical_metadata", []),
            rows_by_manifest.get("file_manifest", []),
        )
    )

    if all_errors:
        print("Manifest validation failed:", file=sys.stderr)
        for error in all_errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    total_rows = sum(len(rows) for rows in rows_by_manifest.values())
    print(f"Manifest validation passed for {len(paths)} files and {total_rows} data rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
