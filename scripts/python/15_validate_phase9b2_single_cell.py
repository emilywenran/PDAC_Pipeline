#!/usr/bin/env python3
"""Validate Phase 9B2-primary single-cell execution guardrails and outputs."""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "05_results/tables"
FIG_DIR = ROOT / "05_results/figures"
DATA_DIR = ROOT / "02_data/external/phase9_single_cell/PENG_CRA001160"
AUTHORITATIVE = ROOT / "05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv"
INVENTORY = ROOT / "01_metadata/external_validation_dataset_inventory.tsv"
REPORT = ROOT / "04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md"

EXPECTED_PRIMARY = {"PENG_CRA001160"}
EXPECTED_FILES = {
    "phase9b2_restart_runtime_validation.tsv",
    "phase9b2_single_cell_cohort_qc.tsv",
    "phase9b2_cells_per_patient.tsv",
    "phase9b2_sample_exclusions.tsv",
    "phase9b2_cell_annotation_audit.tsv",
    "phase9b2_malignant_cell_audit.tsv.gz",
    "phase9b2_single_cell_feature_coverage.tsv",
    "phase9b2_pseudobulk_inventory.tsv",
    "phase9b2_patient_celltype_expression_qc.tsv",
    "phase9b2_patient_celltype_state_scores.tsv",
    "phase9b2_cell_state_scores.tsv.gz",
    "phase9b2_patient_celltype_host_program_scores.tsv",
    "phase9b2_cell_host_program_scores.tsv.gz",
    "phase9b2_module_transfer_coverage.tsv",
    "phase9b2_patient_celltype_tf_activity.tsv",
    "phase9b2_cell_tf_activity_scores.tsv.gz",
    "phase9b2_tf_regulon_coverage.tsv",
    "phase9b2_cellular_source_models.tsv",
    "phase9b2_malignant_state_heterogeneity.tsv",
    "phase9b2_malignant_feature_axis_associations.tsv",
    "phase9b2_cell_composition_sensitivity.tsv",
    "phase9b2_tumor_control_descriptive.tsv",
    "phase9b2_negative_control_results.tsv",
    "phase9b2_cellular_source_evidence.tsv",
}
EXPECTED_FIGURES = {
    "phase9b2_cohort_cell_counts.pdf",
    "phase9b2_cell_annotation_markers.pdf",
    "phase9b2_malignant_cell_audit.pdf",
    "phase9b2_moffitt_axis_by_cell_type.pdf",
    "phase9b2_malignant_axis_by_patient.pdf",
    "phase9b2_hallmark_cellular_source.pdf",
    "phase9b2_module_cellular_source.pdf",
    "phase9b2_tf_activity_cellular_source.pdf",
    "phase9b2_malignant_feature_axis_heatmap.pdf",
    "phase9b2_cell_composition_sensitivity.pdf",
    "phase9b2_tumor_control_descriptive.pdf",
    "phase9b2_negative_control_summary.pdf",
    "phase9b2_cellular_source_evidence_summary.pdf",
}
LOCKED_FEATURES = {
    "HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS",
    "MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow",
    "CTCFL", "IRF3", "JUNB", "KLF13", "KLF9", "MNT", "MXI1", "SNAI2",
    "TFAP4", "TP63", "ZBTB7A", "ZNF24",
    "BHLHE40", "E2F6", "ELF1", "GRHL2", "KLF1", "MBD1", "MBD2",
    "OTX2", "SIX5", "SNAPC4", "ZBED1", "ZNF384", "ZNF740",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def main() -> int:
    errors: list[str] = []

    auth = read_tsv(AUTHORITATIVE)
    primary = {r["canonical_dataset_id"] for r in auth if r.get("included_in_phase9b2_primary", "").upper() == "TRUE"}
    if primary != EXPECTED_PRIMARY:
        fail(errors, f"Primary Phase 9B2 set is {primary}, expected {EXPECTED_PRIMARY}.")

    inv = [r for r in read_tsv(INVENTORY) if r.get("canonical_dataset_id") == "PENG_CRA001160"]
    if len(inv) != 1:
        fail(errors, f"Expected one PENG_CRA001160 inventory row, observed {len(inv)}.")
    else:
        row = inv[0]
        for field, expected in {"accession": "CRA001160", "BioProject": "PRJCA001063", "PDAC_patients": "24", "control_patients": "11"}.items():
            if row.get(field) != expected:
                fail(errors, f"PENG_CRA001160 {field}={row.get(field)} expected {expected}.")
        if row.get("accession") == "GSE111672":
            fail(errors, "Peng cohort is incorrectly labelled as GSE111672.")

    raw_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram", ".sra")
    raw = [str(p.relative_to(ROOT)) for p in DATA_DIR.rglob("*") if p.is_file() and p.name.lower().endswith(raw_suffixes)]
    if raw:
        fail(errors, "Raw sequencing files detected: " + "; ".join(raw[:10]))

    for name in EXPECTED_FILES:
        path = TABLE_DIR / name
        if not path.exists() or path.stat().st_size == 0:
            fail(errors, f"Missing or empty required table: {path.relative_to(ROOT)}")
    for name in EXPECTED_FIGURES:
        path = FIG_DIR / name
        if not path.exists() or path.stat().st_size == 0:
            fail(errors, f"Missing or empty required figure: {path.relative_to(ROOT)}")

    feature_files = [
        TABLE_DIR / "phase9b2_cellular_source_models.tsv",
        TABLE_DIR / "phase9b2_malignant_feature_axis_associations.tsv",
        TABLE_DIR / "phase9b2_cellular_source_evidence.tsv",
    ]
    observed = set()
    for path in feature_files:
        if path.exists():
            for r in read_tsv(path):
                if r.get("feature_name"):
                    observed.add(r["feature_name"])
    missing = LOCKED_FEATURES - observed
    forbidden = {"MEred", "MEpurple", "GFI1B", "STAT1", "ZBTB11", "ZNF639", "TWIST1", "FOXK2", "KDM5B", "MAFF", "TEAD4"} & observed
    if missing:
        fail(errors, "Locked features missing from model/evidence outputs: " + ", ".join(sorted(missing)))
    if forbidden:
        fail(errors, "Not-replicated features unexpectedly analyzed: " + ", ".join(sorted(forbidden)))

    tf_path = TABLE_DIR / "phase9b2_patient_celltype_tf_activity.tsv"
    if tf_path.exists():
        tf_text = tf_path.read_text(encoding="utf-8", errors="ignore")[:5000]
        if "expression_proxy" in tf_text.lower() or "tf_symbol_expression" in tf_text.lower():
            fail(errors, "TF output indicates expression-proxy usage.")

    pb = TABLE_DIR / "phase9b2_pseudobulk_inventory.tsv"
    if pb.exists():
        rows = read_tsv(pb)
        if not rows or "patient_id" not in rows[0] or "number_of_cells" not in rows[0]:
            fail(errors, "Pseudobulk inventory lacks patient-level fields.")

    if REPORT.exists():
        text = REPORT.read_text(encoding="utf-8", errors="ignore").lower()
        banned = ["microbiome replication", "ochrobactrum was validated", "causal proof", "causally demonstrates"]
        for term in banned:
            if term in text:
                fail(errors, f"Report contains prohibited interpretation phrase: {term}")
    else:
        fail(errors, "Missing Phase 9B2 report.")

    if errors:
        print("Phase 9B2 validation FAILED:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    print("Phase 9B2 validation passed.")
    print("canonical_dataset_id=PENG_CRA001160")
    print("result=READY_FOR_INDEPENDENT_REVIEW")
    return 0


if __name__ == "__main__":
    sys.exit(main())
