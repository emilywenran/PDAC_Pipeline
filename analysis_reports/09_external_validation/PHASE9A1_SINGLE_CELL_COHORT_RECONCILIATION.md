# Phase 9A.1: Layer 2 Single-Cell Validation Cohort Reconciliation Report

> [!WARNING]
> **SUPERSEDED:** The conclusions of this Phase 9A.1 report have been fully superseded by the Phase 9A.2 provenance correction and cohort set expansion detailed in [PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION.md](file://~/thesis/PDAC/04_analysis/09_external_validation/PHASE9A2_SINGLE_CELL_DATASET_PROVENANCE_CORRECTION.md).

## 1. Executive Summary

This report completes Phase 9A.1: reconciliation of the locked `PRIORITY_1` single-cell cohort set prior to Phase 9B2 execution. The Phase 9B2 startup audit identified a record disagreement in the locked Layer 2 cohorts across three primary Phase 9A administrative and results files:
- `01_metadata/external_validation_dataset_inventory.tsv` (GSE111672, GSE154778, GSE202051)
- `05_results/tables/phase9a_external_dataset_shortlist.tsv` (GSE111672 only)
- `01_metadata/external_validation_parameter_inventory.tsv` (GSE111672 only)

This reconciliation report independently re-audits the three candidate datasets, assigns authoritative suitability statuses, resolves the database record inconsistencies, and documents the Phase 9B2 startup-stop event.

**Final Readiness Status:** **`SUPERSEDED_BY_PHASE9A2_PROVENANCE_CORRECTION`**

---

## 2. Inconsistency Analysis & Root Cause

The root cause of the cohort record disagreement is twofold:
1. **Broader Candidate Inventory Rows:** The full dataset inventory (`external_validation_dataset_inventory.tsv`) contained candidate datasets that were initially reviewed during the discovery phase but were not selected or parameter-locked for the final execution set in Phase 9A. These were incorrectly retained as `PRIORITY_1` for Layer 2 in the full inventory.
2. **Metadata and Accession Mapping Errors:**
   - **GSE111672 (Peng et al. 2019 / GSA CRA001160):** The dataset from Peng et al. (2019) is a foundational single-cell cohort of 24 patients, but its primary data files are hosted under the National Genomics Data Center (NGDC) Genome Sequence Archive (GSA) accession `CRA001160` (BioProject `PRJCA001063`). In the project metadata, it was assigned the accession `GSE111672`, which in GEO actually belongs to Moncada et al. (2020) (a spatial/single-cell dataset of only 2 patient tumors).
   - **GSE154778 (Lin et al. 2020 / Steele et al. 2020):** In GEO, `GSE154778` belongs to Lin et al. (2020) (16 samples). However, the inventory listed it as Steele et al. (2020) (Nature Cancer, PMID: 35122046). Steele et al. (2020) deposited their data in dbGaP under accession `phs002071`, which is controlled-access and flow-sorted for CD45+ immune cells.

---

## 3. Task 1 & 2: Authoritative Audit and Suitability Assignments

The three candidate datasets have been independently re-audited according to the 14 prospective criteria:

| Audit Parameter | GSE111672 (Peng / GSA CRA001160) | GSE154778 (Lin / Steele 2020) | GSE202051 (Lin et al. 2023) |
| :--- | :--- | :--- | :--- |
| **Primary Human PDAC Tissue** | Yes | Yes (Lin 2020) / Yes (Steele 2020) | Yes |
| **Data Modality** | Single-cell RNA-seq (scRNA-seq) | Single-cell RNA-seq (scRNA-seq) | Single-nucleus RNA-seq (snRNA-seq) |
| **Independent Patients** | Yes | Yes (16 samples) | Yes (43 patients) |
| **PDAC Patients / Tumors** | 24 patients / 24 tumors | 16 samples / 16 tumors | 43 patients / 43 tumors |
| **Processed Matrix** | Yes (CNCB-NGDC GSA CRA001160) | Yes (for Lin 2020) / dbGaP restricted (Steele 2020) | Yes (GEO GSE202051) |
| **Metadata Available** | Yes | Yes | Yes |
| **Patient Identifiers** | Yes | Yes | Yes |
| **Cell Annotations** | Yes | Yes | Yes |
| **Malignant-Cell Labels** | Yes | Yes (Lin 2020) / No (Steele 2020 is CD45+ sorted) | Yes |
| **Treatment Status** | Treatment-naïve | Mixed | Mixed (21 naïve, 22 neoadjuvant-treated) |
| **Discovery Overlap** | No | No | No |
| **Patient-Aware Pseudobulk** | Yes (suitable) | Yes (Lin 2020) / No (Steele 2020) | Yes (suitable) |
| **Locked 9B2 Host-Features** | Yes (feasible on corrected CRA001160) | No (Steele 2020 lacks epithelial cells) | Yes (feasible but redundant with Layer 3) |
| **Authoritative Status** | **`PRIORITY_1`** | **`NOT_SUITABLE`** | **`NOT_SUITABLE`** (for Layer 2) |
| **Exclusion/Downgrade Reason** | None. Retained as primary execution cohort (data files mapped via CRA001160). | Controlled-access dbGaP restriction prevents public download. Sort strategy (CD45+) eliminates malignant epithelial cells. | Already locked as primary Layer 3 Spatial dataset. Excluded from Layer 2 to avoid duplicate testing and neoadjuvant confounding. |

### Exclusion Rationale Summary
- **GSE154778:** Downgraded to `NOT_SUITABLE`. The publication Steele et al. (2020) mapped in the inventory requires dbGaP authorized access (`phs002071`) which prohibits automated public data acquisition. Operationally, the study sorted for CD45+ immune cells, which means it contains almost no malignant epithelial cells. This makes cell-type enrichment contrasts and basal-classical scoring of epithelial cells impossible.
- **GSE202051:** Classified as `PRIORITY_1` for **Layer 3 (Spatial)** but downgraded to `NOT_SUITABLE` for **Layer 2 (Single-Cell)**. Because it contains both single-nucleus and spatial transcriptomics, its single-nucleus portion is technically a candidate for Layer 2. However, it was already shortlisted and parameter-locked under Layer 3. Evaluating it in Layer 2 would result in duplicate testing of the same patient series. Additionally, 22 out of 43 patients received neoadjuvant chemotherapy/radiotherapy, which confounds baseline host expression profiles. Single-nucleus sequencing also introduces differences in transcript capture compared to single-cell RNA-seq.

---

## 4. Task 3: Authoritative Phase 9B2 Cohort Set Definition

The final, reconciled Phase 9B2 Layer 2 cohort set is defined as containing exactly one dataset: **`GSE111672`** (representing Peng et al. 2019, accessed via CNCB-NGDC GSA `CRA001160` / BioProject `PRJCA001063` files to maintain cohort integrity).

This design is now programmatically and textually aligned across all 5 validation files:
1. `01_metadata/external_validation_dataset_inventory.tsv` (Row updated: GSE154778 downgraded; GSE202051 layer reassigned to Layer 3).
2. `05_results/tables/phase9a_external_dataset_shortlist.tsv` (Remains GSE111672 only for Layer 2; GSE202051 remains Layer 3).
3. `01_metadata/external_validation_parameter_inventory.tsv` (Only `VAL_SC_CELLULAR_SOURCE_PENG` is locked for Layer 2).
4. `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md` (Updated Section 8 to lock GSE111672 only).
5. `09_docs/methods/PDAC_external_validation_protocol.md` (Updated Section 3 to specify GSE111672 as sole cohort).

---

## 5. Task 4: Phase 9B2 Startup-Stop Audit Trail

The Phase 9B2 startup stop event is a valid methodological guardrail. It occurred before any single-cell data acquisition, extraction, or scoring took place, and has been preserved in:
- `09_docs/planning/DECISION_LOG.md` (Decision D-33)
- `04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md`

### Stopped State Verification
- **Status:** `STOPPED_BEFORE_DATA_ACQUISITION`
- **Reason:** Inconsistent locked single-cell cohort definitions across Phase 9A records.
- **Acquisition Check:** No processed single-cell matrices, cell metadata, or Seurat/Scanpy objects were downloaded. No raw FASTQ, BAM, or SRA files were downloaded.
- **Analysis Check:** No biological analysis, cell-type scoring (ssGSEA/decoupleR), TF activity inference (VIPER), pseudobulking, patient-aware mixed-effects modeling, negative-control permutation, spatial mapping, or microbiome association was performed.

---

## 6. Implementation Verification

A Python verification run was performed:
- `06_scripts/python/00_validate_manifests.py` passed successfully.
- `06_scripts/python/15_validate_phase9b2_single_cell.py` confirmed that the stopped-state guardrails were fully respected.
- The new consistency check confirming that all Phase 9A records contain the same Layer 2 cohort set now passes.

With the metadata records fully reconciled, Phase 9B2 may safely restart using the corrected single-cohort design.

*Date: 2026-07-03*
*Reviewer: Antigravity*
