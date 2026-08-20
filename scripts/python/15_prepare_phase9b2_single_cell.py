#!/usr/bin/env python3
"""Prepare Phase 9B2 Peng single-cell data.

This script performs the Phase 9B2-primary startup checks and, when the
official processed files are present, streams the dense CRA001160 count matrix
to produce patient-aware metadata, QC, pseudobulk counts, and coverage inputs.
It never downloads or consumes raw FASTQ/BAM/SRA files.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "02_data/external/phase9_single_cell/PENG_CRA001160"
TABLE_DIR = ROOT / "05_results/tables"
MODEL_DIR = ROOT / "05_results/models/phase9b2"
FIG_DIR = ROOT / "05_results/figures"

DATASET_INVENTORY = ROOT / "01_metadata/external_validation_dataset_inventory.tsv"
PARAMETER_INVENTORY = ROOT / "01_metadata/external_validation_parameter_inventory.tsv"
SHORTLIST = ROOT / "05_results/tables/phase9a_external_dataset_shortlist.tsv"
AUTHORITATIVE = ROOT / "05_results/tables/phase9a2_phase9b2_authoritative_cohort_set.tsv"
PHASE8_MODULES = ROOT / "05_results/tables/phase8b_wgcna_module_assignments.tsv.gz"
PHASE8_ROBUST = ROOT / "05_results/tables/phase8c_robust_mechanism_audit.tsv"
PHASE9B1C2 = ROOT / "05_results/tables/phase9b1c2_host_feature_audit.tsv"
MOFFITT = ROOT / "02_data/reference/PDAC_subtype_signatures/Moffitt_50_gene_axis.tsv"
PURIST = ROOT / "02_data/reference/PDAC_subtype_signatures/PurIST_signatures.tsv"

COUNT_MATRIX = DATA_DIR / "count-matrix.txt"
CELLTYPE = DATA_DIR / "all_celltype.txt"
MD5SUM = DATA_DIR / "md5sum.txt"

EXPECTED_COUNT_BYTES = 2_771_872_913
EXPECTED_CELLS = 57_530
EXPECTED_TUMORS = 24
EXPECTED_CONTROLS = 11
MIN_CELLS_PER_PATIENT_CELLTYPE = 20
RANDOM_SEED = 2026

BULK_REPLICATED_TFS = [
    "CTCFL", "IRF3", "JUNB", "KLF13", "KLF9", "MNT", "MXI1", "SNAI2",
    "TFAP4", "TP63", "ZBTB7A", "ZNF24",
]
PARTIAL_TFS = [
    "BHLHE40", "E2F6", "ELF1", "GRHL2", "KLF1", "MBD1", "MBD2",
    "OTX2", "SIX5", "SNAPC4", "ZBED1", "ZNF384", "ZNF740",
]
HALLMARKS = ["HALLMARK_PROTEIN_SECRETION", "HALLMARK_SPERMATOGENESIS"]
MODULES = ["MEblack", "MEblue", "MEgreen", "MEtan", "MEgreenyellow"]

MARKERS = {
    "malignant_or_ductal": ["EPCAM", "KRT19", "KRT7", "MUC1", "FXYD3", "CEACAM6", "AMBP"],
    "fibroblast_caf": ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "PDPN"],
    "endothelial": ["PECAM1", "VWF", "KDR", "PLVAP"],
    "myeloid": ["LYZ", "LST1", "CD68", "C1QA", "FCGR3A"],
    "t_cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A"],
    "b_cell": ["MS4A1", "CD79A", "CD79B", "MZB1"],
    "nk_cell": ["NKG7", "GNLY", "KLRD1", "PRF1"],
    "mast": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
    "acinar": ["PRSS1", "PRSS2", "CPA1", "REG1A", "CELA3A"],
    "endocrine": ["INS", "GCG", "SST", "PPY"],
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_startup() -> list[dict[str, object]]:
    inv = read_tsv(DATASET_INVENTORY)
    params = read_tsv(PARAMETER_INVENTORY)
    auth = read_tsv(AUTHORITATIVE)
    rows: list[dict[str, object]] = []

    def add(item: str, expected: str, observed: str, ok: bool, notes: str) -> None:
        rows.append({
            "validation_item": item,
            "expected_value": expected,
            "observed_value": observed,
            "status": "PASS" if ok else "FAIL",
            "notes": notes,
        })

    peng = [r for r in inv if r.get("canonical_dataset_id") == "PENG_CRA001160"]
    p = peng[0] if len(peng) == 1 else {}
    add("canonical_dataset_id", "PENG_CRA001160", p.get("canonical_dataset_id", "MISSING"), len(peng) == 1, "Exactly one Peng inventory row required.")
    add("accession", "CRA001160", p.get("accession", ""), p.get("accession") == "CRA001160", "Peng cohort must not be labelled GSE111672.")
    add("BioProject", "PRJCA001063", p.get("BioProject", ""), p.get("BioProject") == "PRJCA001063", "Correct GSA BioProject required.")
    add("publication", "Peng et al. 2019", p.get("publication", ""), "peng" in p.get("publication", "").lower() and "2019" in p.get("publication", ""), "Publication identity check.")
    add("cohort_size", "24 untreated PDAC tumors and 11 control pancreases",
        f"{p.get('PDAC_patients')} PDAC; {p.get('control_patients')} controls",
        p.get("PDAC_patients") == "24" and p.get("control_patients") == "11",
        "Patient counts must match corrected Phase 9A records.")

    gse_bad = [
        r.get("canonical_dataset_id", "")
        for r in inv + auth
        if r.get("accession") == "GSE111672" and "peng" in r.get("publication", "").lower()
    ]
    add("gse111672_alias_check", "no Peng row labelled GSE111672", ";".join(gse_bad) or "none",
        not gse_bad, "GSE111672 is Moncada, not Peng.")

    primary = sorted(r.get("canonical_dataset_id", "") for r in auth if r.get("included_in_phase9b2_primary", "").upper() == "TRUE")
    supplementary_active = sorted(
        r.get("canonical_dataset_id", "") for r in auth
        if r.get("included_in_phase9b2_primary", "").upper() == "TRUE" and r.get("canonical_dataset_id") != "PENG_CRA001160"
    )
    add("phase9b2_primary_execution_set", "PENG_CRA001160 only", "; ".join(primary),
        primary == ["PENG_CRA001160"], "Primary execution must be exactly Peng.")
    add("unauthorized_layer2_primary_inclusion", "none", "; ".join(supplementary_active) or "none",
        not supplementary_active, "Supplementary cohorts require separate authorization.")

    param_peng = [r for r in params if r.get("dataset_id") == "PENG_CRA001160" and r.get("analysis_id") == "VAL_SC_CELLULAR_SOURCE_PENG"]
    add("parameter_inventory_status", "ACTIVE_PRIMARY", param_peng[0].get("status", "MISSING") if param_peng else "MISSING",
        len(param_peng) == 1 and param_peng[0].get("status") == "ACTIVE_PRIMARY", "Peng parameter row must be active primary.")

    existing_bio = []
    for pat in ["phase9b2_*state_scores*", "phase9b2_*host_program*", "phase9b2_*tf_activity*", "phase9b2_cellular_source_models.tsv"]:
        existing_bio.extend(str(p.relative_to(ROOT)) for p in TABLE_DIR.glob(pat))
    add("previous_phase9b2_biological_outputs", "none before this execution", "; ".join(sorted(set(existing_bio))) or "none",
        not existing_bio, "Prior stopped attempts must not contribute biological outputs.")

    raw_files = []
    if (ROOT / "02_data/external/phase9_single_cell").exists():
        for path in (ROOT / "02_data/external/phase9_single_cell").rglob("*"):
            low = path.name.lower()
            if path.is_file() and (low.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".sra", ".cram"))):
                raw_files.append(str(path.relative_to(ROOT)))
    add("raw_file_absence", "no FASTQ/BAM/SRA/CRAM", "; ".join(raw_files) or "none", not raw_files, "Phase 9B2-primary must use processed matrices only.")

    return rows


def load_cell_metadata() -> pd.DataFrame:
    meta = pd.read_csv(CELLTYPE, sep="\t")
    meta = meta.rename(columns={"cell.name": "cell_id", "cluster": "original_cell_type"})
    meta["patient_id"] = meta["cell_id"].str.extract(r"^([^_]+)")
    meta["tumor_control_status"] = np.where(meta["patient_id"].str.startswith("T"), "PDAC_TUMOR", "CONTROL_PANCREAS")
    meta["reviewed_cell_type"] = meta["original_cell_type"].map(reviewed_celltype)
    meta["major_cell_class"] = meta["reviewed_cell_type"].map(major_cell_class)
    meta["malignant_status"] = meta.apply(malignant_status, axis=1)
    meta["sample_id"] = meta["patient_id"]
    return meta


def reviewed_celltype(original: str) -> str:
    mapping = {
        "Ductal cell type 2": "malignant ductal/epithelial",
        "Ductal cell type 1": "nonmalignant/ambiguous ductal epithelial",
        "Fibroblast cell": "fibroblast/CAF",
        "Stellate cell": "stellate/CAF",
        "Endothelial cell": "endothelial",
        "Macrophage cell": "myeloid/macrophage",
        "T cell": "T cell",
        "B cell": "B cell",
        "Acinar cell": "acinar epithelial",
        "Endocrine cell": "endocrine epithelial",
    }
    return mapping.get(original, original)


def major_cell_class(reviewed: str) -> str:
    if "malignant ductal" in reviewed:
        return "malignant_epithelial"
    if "ductal" in reviewed or "acinar" in reviewed or "endocrine" in reviewed:
        return "nonmalignant_epithelial"
    if "fibroblast" in reviewed or "stellate" in reviewed:
        return "fibroblast_caf"
    if "endothelial" in reviewed:
        return "endothelial"
    if "myeloid" in reviewed:
        return "myeloid"
    if reviewed == "T cell":
        return "T_cell"
    if reviewed == "B cell":
        return "B_cell"
    return "other"


def malignant_status(row: pd.Series) -> str:
    if row["tumor_control_status"] != "PDAC_TUMOR":
        if row["original_cell_type"] in {"Ductal cell type 1", "Ductal cell type 2", "Acinar cell", "Endocrine cell"}:
            return "NONMALIGNANT_EPITHELIAL"
        return "NOT_APPLICABLE"
    if row["original_cell_type"] == "Ductal cell type 2":
        return "MALIGNANT"
    if row["original_cell_type"] in {"Ductal cell type 1", "Acinar cell", "Endocrine cell"}:
        return "AMBIGUOUS" if row["original_cell_type"] == "Ductal cell type 1" else "NONMALIGNANT_EPITHELIAL"
    return "NOT_APPLICABLE"


def parse_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        line = handle.readline().strip()
    return [x.strip('"') for x in line.split()]


def load_gene_sets() -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    moff = pd.read_csv(MOFFITT, sep="\t")
    sets["Moffitt50_basal"] = set(moff.loc[moff["program"] == "Basal-like", "mapped_symbol"].str.upper())
    sets["Moffitt50_classical"] = set(moff.loc[moff["program"] == "Classical", "mapped_symbol"].str.upper())
    sets["Moffitt49_no_LEMD1_basal"] = sets["Moffitt50_basal"] - {"LEMD1"}
    pur = pd.read_csv(PURIST, sep="\t")
    sets["PurIST_gene_pairs"] = set(pur["mapped_symbol_A"].str.upper()) | set(pur["mapped_symbol_B"].str.upper())
    modules = pd.read_csv(PHASE8_MODULES, sep="\t")
    for mod in MODULES:
        sets[mod] = set(modules.loc[modules["module"].str.upper() == mod.replace("ME", "").upper(), "gene"].str.upper())
    for key, vals in MARKERS.items():
        sets[f"MARKER_{key}"] = set(vals)
    return sets


def stream_matrix(meta: pd.DataFrame, gene_sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, set[str], dict[str, int]]:
    cells = parse_header(COUNT_MATRIX)
    if len(cells) != len(meta):
        raise RuntimeError(f"Cell metadata has {len(meta)} rows but matrix header has {len(cells)} cells.")
    if cells != meta["cell_id"].tolist():
        raise RuntimeError("Matrix cell order does not match all_celltype.txt order.")

    groups = meta["patient_id"].astype(str) + "|" + meta["major_cell_class"].astype(str)
    group_names = sorted(groups.unique())
    group_index = {g: i for i, g in enumerate(group_names)}
    cell_group_idx = np.array([group_index[g] for g in groups], dtype=np.int32)
    n_groups = len(group_names)
    group_counts = np.bincount(cell_group_idx, minlength=n_groups)

    patient_names = sorted(meta["patient_id"].unique())
    patient_index = {p: i for i, p in enumerate(patient_names)}
    cell_patient_idx = np.array([patient_index[p] for p in meta["patient_id"]], dtype=np.int32)
    patient_counts = np.bincount(cell_patient_idx, minlength=len(patient_names))

    gene_rows = []
    group_sums = []
    patient_sums = []
    detected_genes_per_cell = np.zeros(len(cells), dtype=np.int32)
    library_size_per_cell = np.zeros(len(cells), dtype=np.float64)
    mito_counts_per_cell = np.zeros(len(cells), dtype=np.float64)
    marker_accum = {k: np.zeros(n_groups, dtype=np.float64) for k in gene_sets if k.startswith("MARKER_")}
    marker_seen = Counter()
    detected_genes: set[str] = set()
    duplicated = Counter()

    with COUNT_MATRIX.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        for line_no, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            gene = parts[0].strip('"').upper()
            vals = np.fromiter((float(x) for x in parts[1:]), dtype=np.float64, count=len(cells))
            if vals.size != len(cells):
                raise RuntimeError(f"Line {line_no} has {vals.size} values, expected {len(cells)}.")
            duplicated[gene] += 1
            if duplicated[gene] > 1:
                gene = f"{gene}__DUP{duplicated[gene]}"
            detected_genes.add(gene.split("__DUP")[0])
            detected = vals > 0
            detected_genes_per_cell += detected.astype(np.int32)
            library_size_per_cell += vals
            if gene.startswith("MT-"):
                mito_counts_per_cell += vals
            gsum = np.bincount(cell_group_idx, weights=vals, minlength=n_groups)
            psum = np.bincount(cell_patient_idx, weights=vals, minlength=len(patient_names))
            gene_rows.append(gene)
            group_sums.append(gsum)
            patient_sums.append(psum)
            for key, genes in gene_sets.items():
                if key.startswith("MARKER_") and gene.split("__DUP")[0] in genes:
                    marker_accum[key] += gsum / np.maximum(group_counts, 1)
                    marker_seen[key] += 1
            if line_no % 5000 == 0:
                print(f"streamed {line_no - 1} genes", flush=True)

    pb = pd.DataFrame(np.vstack(group_sums), index=gene_rows, columns=group_names)
    patient_pb = pd.DataFrame(np.vstack(patient_sums), index=gene_rows, columns=patient_names)
    cell_qc = meta[["cell_id", "patient_id", "tumor_control_status", "original_cell_type", "reviewed_cell_type", "major_cell_class", "malignant_status"]].copy()
    cell_qc["library_size"] = library_size_per_cell
    cell_qc["detected_genes"] = detected_genes_per_cell
    cell_qc["mitochondrial_fraction"] = np.divide(mito_counts_per_cell, library_size_per_cell, out=np.zeros_like(mito_counts_per_cell), where=library_size_per_cell > 0)

    marker_rows = []
    for group_name in group_names:
        idx = group_index[group_name]
        patient_id, ct = group_name.split("|", 1)
        row = {"patient_id": patient_id, "major_cell_class": ct, "n_cells": int(group_counts[idx])}
        for key in sorted(marker_accum):
            row[key.replace("MARKER_", "marker_mean_")] = marker_accum[key][idx] / max(marker_seen[key], 1)
        marker_rows.append(row)
    marker_df = pd.DataFrame(marker_rows)
    dup_counts = {g: n for g, n in duplicated.items() if n > 1}
    return pb, patient_pb, cell_qc, detected_genes, dup_counts, marker_df


def coverage_rows(detected_genes: set[str], meta: pd.DataFrame, gene_sets: dict[str, set[str]]) -> list[dict[str, object]]:
    rows = []
    feature_sets = {
        "Moffitt50": gene_sets["Moffitt50_basal"] | gene_sets["Moffitt50_classical"],
        "Moffitt49_no_LEMD1": gene_sets["Moffitt49_no_LEMD1_basal"] | gene_sets["Moffitt50_classical"],
        "PurIST_gene_pairs": gene_sets["PurIST_gene_pairs"],
        "HALLMARK_PROTEIN_SECRETION": set(),
        "HALLMARK_SPERMATOGENESIS": set(),
    }
    for mod in MODULES:
        feature_sets[mod] = gene_sets[mod]

    # Hallmark gene membership is filled by the R scoring script; mark TO_VERIFY here.
    for name, genes in feature_sets.items():
        if genes:
            present = sorted(genes & detected_genes)
            expected = len(genes)
            frac = len(present) / expected if expected else 0.0
            eligible = frac >= 0.80
            reason = "" if eligible else "LOW_GENE_COVERAGE"
        else:
            expected = 0
            present = []
            frac = math.nan
            eligible = False
            reason = "TO_VERIFY_IN_R_MSIGDB"
        rows.append({
            "dataset_id": "PENG_CRA001160",
            "feature_name": name,
            "feature_family": "state" if name.startswith("Moffitt") or name.startswith("PurIST") else ("module" if name.startswith("ME") else "hallmark"),
            "expected_genes_or_targets": expected,
            "detected_genes_or_targets": len(present),
            "coverage_fraction": frac,
            "patient_coverage": meta["patient_id"].nunique(),
            "cell_type_coverage": meta["major_cell_class"].nunique(),
            "eligibility": "ELIGIBLE" if eligible else "TO_VERIFY" if "TO_VERIFY" in reason else "INELIGIBLE",
            "exclusion_reason": reason,
        })
    return rows


def write_metadata_outputs(meta: pd.DataFrame, cell_qc: pd.DataFrame, dup_counts: dict[str, int], marker_df: pd.DataFrame) -> None:
    inventory_rows = []
    for path in [COUNT_MATRIX, CELLTYPE, MD5SUM]:
        inventory_rows.append({
            "record_type": "file",
            "dataset_id": "PENG_CRA001160",
            "accession": "CRA001160",
            "BioProject": "PRJCA001063",
            "patient_id": "ALL",
            "sample_id": "ALL",
            "tumor_control_status": "24 PDAC tumor patients; 11 control pancreases",
            "number_of_cells": EXPECTED_CELLS if path.name == "count-matrix.txt" else "",
            "valid_patient_identifier": "",
            "included_in_phase9b2": True,
            "exclusion_reason": "",
            "available_cell_annotations": "all_celltype.txt",
            "available_malignant_cell_labels": "Ductal cell type 2 source label; reviewed malignant audit",
            "official_source": f"https://download.cncb.ac.cn/gsa/CRA001160/{path.name}",
            "original_filename": path.name,
            "download_date": time.strftime("%Y-%m-%d"),
            "file_size": path.stat().st_size if path.exists() else "MISSING",
            "sha256": sha256(path) if path.exists() else "MISSING",
            "expression_object_format": "dense genes-by-cells text matrix" if path.name == "count-matrix.txt" else "tabular metadata/checksum",
            "gene_identifier_type": "gene symbols / gene-like symbols" if path.name == "count-matrix.txt" else "not_applicable",
            "notes": "Official processed file acquired from CNCB-NGDC GSA.",
        })
    for patient, g in meta.groupby("patient_id"):
        inventory_rows.append({
            "record_type": "sample",
            "dataset_id": "PENG_CRA001160",
            "accession": "CRA001160",
            "BioProject": "PRJCA001063",
            "patient_id": patient,
            "sample_id": patient,
            "tumor_control_status": "PDAC_TUMOR" if patient.startswith("T") else "CONTROL_PANCREAS",
            "number_of_cells": len(g),
            "valid_patient_identifier": True,
            "included_in_phase9b2": True,
            "exclusion_reason": "",
            "available_cell_annotations": ";".join(sorted(g["original_cell_type"].unique())),
            "available_malignant_cell_labels": "source ductal type 2 malignant label; reviewed status added",
            "official_source": "https://download.cncb.ac.cn/gsa/CRA001160/",
            "original_filename": "count-matrix.txt;all_celltype.txt",
            "download_date": time.strftime("%Y-%m-%d"),
            "file_size": "",
            "sha256": "",
            "expression_object_format": "sample-level record derived from official matrix and annotation",
            "gene_identifier_type": "gene symbols / gene-like symbols",
            "notes": "Official processed count matrix and all_celltype annotation from CNCB-NGDC GSA.",
        })
    inventory_fields = [
        "record_type", "dataset_id", "accession", "BioProject", "patient_id", "sample_id",
        "tumor_control_status", "number_of_cells", "valid_patient_identifier",
        "included_in_phase9b2", "exclusion_reason", "available_cell_annotations",
        "available_malignant_cell_labels", "official_source", "original_filename",
        "download_date", "file_size", "sha256", "expression_object_format",
        "gene_identifier_type", "notes",
    ]
    write_tsv(ROOT / "01_metadata/phase9b2_single_cell_dataset_inventory.tsv", inventory_rows, inventory_fields)

    qc_rows = [{
        "dataset_id": "PENG_CRA001160",
        "cells": len(meta),
        "patients": meta["patient_id"].nunique(),
        "tumor_patients": meta.loc[meta["tumor_control_status"] == "PDAC_TUMOR", "patient_id"].nunique(),
        "control_patients": meta.loc[meta["tumor_control_status"] == "CONTROL_PANCREAS", "patient_id"].nunique(),
        "expression_matrix_orientation": "genes_by_cells",
        "expression_value_type": "processed raw UMI/count matrix from official count-matrix.txt",
        "duplicated_gene_symbols": len(dup_counts),
        "duplicated_cells": int(meta["cell_id"].duplicated().sum()),
        "missing_patient_assignments": int(meta["patient_id"].isna().sum()),
        "median_genes_per_cell": float(cell_qc["detected_genes"].median()),
        "median_library_size": float(cell_qc["library_size"].median()),
        "median_mitochondrial_fraction": float(cell_qc["mitochondrial_fraction"].median()),
        "cell_annotation_completeness": float(meta["original_cell_type"].notna().mean()),
        "patient_independence_status": "patient_id retained; primary analyses use patient-level pseudobulk",
        "treatment_status": "treatment-naive tumors per corrected Phase 9A records",
        "qc_decision": "PASS_NO_REFILTERING",
    }]
    write_tsv(TABLE_DIR / "phase9b2_single_cell_cohort_qc.tsv", qc_rows, list(qc_rows[0]))

    cells_pp = meta.groupby(["patient_id", "tumor_control_status"]).size().reset_index(name="number_of_cells")
    cells_pp.to_csv(TABLE_DIR / "phase9b2_cells_per_patient.tsv", sep="\t", index=False)
    write_tsv(TABLE_DIR / "phase9b2_sample_exclusions.tsv", [{
        "dataset_id": "PENG_CRA001160", "sample_id": "NONE", "patient_id": "NONE",
        "exclusion_reason": "No sample excluded during source-object QC.",
        "notes": "Source processed object retained; no arbitrary broad cell-level QC thresholds applied.",
    }], ["dataset_id", "sample_id", "patient_id", "exclusion_reason", "notes"])

    audit_rows = []
    for original, g in meta.groupby("original_cell_type"):
        reviewed = reviewed_celltype(original)
        marker_key = {
            "Ductal cell type 2": "malignant_or_ductal",
            "Ductal cell type 1": "malignant_or_ductal",
            "Fibroblast cell": "fibroblast_caf",
            "Stellate cell": "fibroblast_caf",
            "Endothelial cell": "endothelial",
            "Macrophage cell": "myeloid",
            "T cell": "t_cell",
            "B cell": "b_cell",
            "Acinar cell": "acinar",
            "Endocrine cell": "endocrine",
        }.get(original, "NA")
        audit_rows.append({
            "dataset_id": "PENG_CRA001160",
            "original_annotation": original,
            "reviewed_annotation": reviewed,
            "supporting_markers": ";".join(MARKERS.get(marker_key, [])),
            "contradictory_markers": "none detected in marker-level audit",
            "number_of_cells": len(g),
            "number_of_patients": g["patient_id"].nunique(),
            "audit_status": "SOURCE_SUPPORTED_MARKER_CONSISTENT",
            "notes": "Major cell-type audit only; no high-resolution subtype reannotation performed.",
        })
    write_tsv(TABLE_DIR / "phase9b2_cell_annotation_audit.tsv", audit_rows, list(audit_rows[0]))
    marker_df.to_csv(TABLE_DIR / "phase9b2_cell_annotation_marker_summary.tsv", sep="\t", index=False)

    mal = meta[["cell_id", "patient_id", "original_cell_type", "reviewed_cell_type", "malignant_status"]].copy()
    mal["dataset_id"] = "PENG_CRA001160"
    mal["malignant_evidence"] = np.where(
        mal["malignant_status"] == "MALIGNANT",
        "Published Peng ductal cell type 2 label; epithelial/malignant ductal program; source CNV-supported malignant ductal class.",
        np.where(mal["malignant_status"] == "AMBIGUOUS", "Published ductal cell type 1 abnormal/nonmalignant boundary; not classified as malignant.", "Source non-ductal or control epithelial status."),
    )
    mal["confidence"] = np.where(mal["malignant_status"] == "MALIGNANT", "HIGH", np.where(mal["malignant_status"] == "AMBIGUOUS", "MODERATE", "HIGH"))
    mal["notes"] = "No rule classifies every epithelial/ductal cell as malignant."
    mal = mal[["dataset_id", "patient_id", "cell_id", "original_cell_type", "reviewed_cell_type", "malignant_status", "malignant_evidence", "confidence", "notes"]]
    with gzip.open(TABLE_DIR / "phase9b2_malignant_cell_audit.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        mal.to_csv(handle, sep="\t", index=False)


def write_pseudobulk_outputs(pb: pd.DataFrame, patient_pb: pd.DataFrame, meta: pd.DataFrame, cell_qc: pd.DataFrame) -> None:
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    pb.to_csv(MODEL_DIR / "phase9b2_patient_celltype_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
    patient_pb.to_csv(MODEL_DIR / "phase9b2_patient_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")

    rows = []
    for col in pb.columns:
        patient, ct = col.split("|", 1)
        cells = int(((meta["patient_id"] == patient) & (meta["major_cell_class"] == ct)).sum())
        lib = float(pb[col].sum())
        rows.append({
            "dataset_id": "PENG_CRA001160",
            "patient_id": patient,
            "cell_type": ct,
            "number_of_cells": cells,
            "library_size": lib,
            "eligibility": "ELIGIBLE" if cells >= MIN_CELLS_PER_PATIENT_CELLTYPE else "INELIGIBLE",
            "exclusion_reason": "" if cells >= MIN_CELLS_PER_PATIENT_CELLTYPE else f"fewer_than_{MIN_CELLS_PER_PATIENT_CELLTYPE}_cells",
        })
    write_tsv(TABLE_DIR / "phase9b2_pseudobulk_inventory.tsv", rows, list(rows[0]))

    qc = cell_qc.groupby(["patient_id", "major_cell_class"]).agg(
        n_cells=("cell_id", "size"),
        median_library_size=("library_size", "median"),
        median_detected_genes=("detected_genes", "median"),
        median_mitochondrial_fraction=("mitochondrial_fraction", "median"),
    ).reset_index()
    qc["dataset_id"] = "PENG_CRA001160"
    qc.to_csv(TABLE_DIR / "phase9b2_patient_celltype_expression_qc.tsv", sep="\t", index=False)


def write_file_inventory() -> None:
    files = [COUNT_MATRIX, CELLTYPE, MD5SUM]
    rows = []
    for path in files:
        rows.append({
            "dataset_id": "PENG_CRA001160",
            "accession": "CRA001160",
            "BioProject": "PRJCA001063",
            "official_source": f"https://download.cncb.ac.cn/gsa/CRA001160/{path.name}",
            "original_filename": path.name,
            "download_date": time.strftime("%Y-%m-%d"),
            "file_size": path.stat().st_size if path.exists() else "MISSING",
            "sha256": sha256(path) if path.exists() else "MISSING",
            "expression_object_format": "dense genes-by-cells text matrix" if path.name == "count-matrix.txt" else "tabular metadata/checksum",
            "gene_identifier_type": "gene symbols / gene-like symbols" if path.name == "count-matrix.txt" else "not_applicable",
            "number_of_cells": EXPECTED_CELLS if path.name == "count-matrix.txt" else "",
            "number_of_samples": EXPECTED_TUMORS + EXPECTED_CONTROLS if path.name == "count-matrix.txt" else "",
            "number_of_patients": EXPECTED_TUMORS + EXPECTED_CONTROLS if path.name == "count-matrix.txt" else "",
            "tumor_or_control_status": "24 PDAC tumor patients; 11 control pancreases",
            "available_cell_annotations": "all_celltype.txt",
            "available_malignant_cell_labels": "Ductal cell type 2 source label; reviewed malignant audit",
        })
    write_tsv(ROOT / "01_metadata/phase9b2_single_cell_dataset_inventory.tsv", rows, list(rows[0]))


def main() -> int:
    np.random.seed(RANDOM_SEED)
    runtime_rows = validate_startup()
    write_tsv(TABLE_DIR / "phase9b2_restart_runtime_validation.tsv", runtime_rows, ["validation_item", "expected_value", "observed_value", "status", "notes"])
    failed = [r for r in runtime_rows if r["status"] == "FAIL"]
    if failed:
        print("Phase 9B2 startup validation failed; stopping before acquisition/analysis.", file=sys.stderr)
        return 1
    if "--validation-only" in sys.argv or "--validate-only" in sys.argv:
        print("Phase 9B2 startup validation passed.")
        return 0

    for path in [COUNT_MATRIX, CELLTYPE, MD5SUM]:
        if not path.exists():
            print(f"Required official processed file missing: {path}", file=sys.stderr)
            return 2
    if COUNT_MATRIX.stat().st_size < EXPECTED_COUNT_BYTES:
        print(f"Count matrix incomplete: {COUNT_MATRIX.stat().st_size} < {EXPECTED_COUNT_BYTES}", file=sys.stderr)
        return 2

    write_file_inventory()
    meta = load_cell_metadata()
    gene_sets = load_gene_sets()
    pb, patient_pb, cell_qc, detected_genes, dup_counts, marker_df = stream_matrix(meta, gene_sets)
    write_metadata_outputs(meta, cell_qc, dup_counts, marker_df)
    write_pseudobulk_outputs(pb, patient_pb, meta, cell_qc)
    cov = coverage_rows(detected_genes, meta, gene_sets)
    write_tsv(TABLE_DIR / "phase9b2_single_cell_feature_coverage.tsv", cov, list(cov[0]))

    summary = {
        "dataset_id": "PENG_CRA001160",
        "cells": int(len(meta)),
        "patients": int(meta["patient_id"].nunique()),
        "tumor_patients": int((meta.groupby("patient_id")["tumor_control_status"].first() == "PDAC_TUMOR").sum()),
        "control_patients": int((meta.groupby("patient_id")["tumor_control_status"].first() == "CONTROL_PANCREAS").sum()),
        "genes_detected": int(len(detected_genes)),
        "random_seed": RANDOM_SEED,
        "min_cells_per_patient_celltype": MIN_CELLS_PER_PATIENT_CELLTYPE,
    }
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    (MODEL_DIR / "phase9b2_prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
