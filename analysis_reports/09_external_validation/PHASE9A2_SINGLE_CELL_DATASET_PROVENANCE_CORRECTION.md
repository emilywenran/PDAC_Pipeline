# Phase 9A.2: Single-Cell Dataset Provenance Correction and Phase 9B2 Cohort Reconciliation

## 1. Executive Summary

This report documents the resolution of metadata and planning inconsistencies in the Layer 2 single-cell and Layer 3 spatial validation cohorts for Phase 9B2. The previous Phase 9A.1 reconciliation process contained errors regarding dataset accessions, publications, controlled-access flags, and compartment availability. 

This correction establishes a rigorous, unified, and scientifically defensible database of candidate cohorts, ensuring complete traceability. 

**Final Readiness Status:** **`READY_WITH_DATASET_LIMITATIONS`** (Phase 9B2 may proceed using the corrected single-cell and spatial validation configurations, acknowledging the specific limitations of each cohort).

---

## 2. Inconsistency Analysis & Core Corrections

### Correction 1: Separation of Peng et al. (2019) and Moncada et al. (2020)

In legacy planning files, the GEO accession `GSE111672` was incorrectly used as a shared identifier or alias for the Peng et al. (2019) cohort. This is biochemically and operationally incorrect. 

We have formally separated these two cohorts into distinct records:
*   **Peng et al. 2019:**
    *   **canonical_dataset_id:** `PENG_CRA001160`
    *   **accession:** `CRA001160` (GSA-CNCB)
    *   **BioProject:** `PRJCA001063`
    *   **Cohort Description:** 24 untreated primary PDAC tumors and 11 control pancreases (~57,000 cells).
    *   **Analysis Role:** Layer 2 primary cellular source validation (treated as `PRIORITY_1`).
*   **Moncada et al. 2020:**
    *   **canonical_dataset_id:** `MONCADA_GSE111672`
    *   **accession:** `GSE111672` (NCBI GEO)
    *   **BioProject:** `PRJNA437847`
    *   **Cohort Description:** Single-cell RNA-seq and spatial transcriptomics from six pancreatic cancer patients.
    *   **Analysis Role:** Layer 3 spatial validation (`PRIORITY_1`) and Layer 2 exploratory single-cell validation (`PRIORITY_2`).

### Correction 2: Re-Audit of GSE154778 (Lin et al. 2020)

The previous reconciliation report incorrectly claimed that `GSE154778` was controlled-access, CD45-sorted-only, and lacked malignant cells. This arose from confusing GSE154778 (Lin et al. 2020 Genome Medicine, PMID: 32988401) with Steele et al. (2020) Nature Cancer (which is controlled-access under dbGaP `phs002071` and sorted for immune cells).

Our re-audit of the official GEO record and original publication for **GSE154778** confirms:
*   **Public processed-data availability:** Digital gene expression matrices are publicly available via `GSE154778_dgeMtx.csv.gz` and MTX files under `GSE154778_RAW.tar` in GEO.
*   **Cohort composition:** 10 primary pancreatic tumors and 6 metastatic biopsies.
*   **Technology:** 10x Genomics Chromium (scRNA-seq).
*   **Compartments:** Profiles both tumor and stromal compartments. Malignant epithelial cells are present and fully available.
*   **Identifiers:** Patient-level identifiers are present, allowing patient-aware pseudobulk.
*   **Sub-cohort analysis:** The 10 primary tumors are from treatment-naïve patients and can be analyzed separately from the 6 metastatic lesions.
*   **Suitability Status:** Assigned as **`PRIORITY_2`** for Layer 2 cellular source validation (included in Phase 9B2).

### Correction 3: Role Classification for GSE202051 (Hwang et al. 2022)

The dataset GSE202051 (PMID: 37277650) includes both single-nucleus RNA-seq and Visium spatial profiling across 43 patients. Rather than excluding it from Layer 2 to avoid duplicate analysis, we have assigned it separate, complementary roles:
1.  **Layer 3 Spatial Validation (`PRIORITY_1`):** Primary cohort to map spatial coordinates and localize continuous host programs.
2.  **Layer 2 Secondary/Treatment-Sensitivity Cellular Source Evaluation (`PRIORITY_2`):** Evaluates cellular programs in a larger cohort that includes neoadjuvant-treated samples (22 treated, 21 naïve), enabling comparative analysis of treatment sensitivity.

---

## 3. Authoritative Phase 9B2 Cohort Set

Based on these audits, the final Phase 9B2 execution cohort set comprises the following single-cell and spatial datasets:

| Canonical Dataset ID | Accession | Publication | PDAC Patients | Primary Tumors | Metastatic Samples | Technology | Validation Layer | Suitability Status | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **PENG_CRA001160** | CRA001160 | Peng et al. (2019) | 24 | 24 | 0 | scRNA-seq (10x) | Layer 2 | `PRIORITY_1` | Primary single-cell cohort; treatment-naïve primary tumors. |
| **LIN_GSE154778** | GSE154778 | Lin et al. (2020) | 10 | 10 | 6 | scRNA-seq (10x) | Layer 2 | `PRIORITY_2` | Secondary single-cell cohort; primary tumors analyzed separately. |
| **MONCADA_GSE111672** | GSE111672 | Moncada et al. (2020) | 6 | 6 | 0 | scRNA-seq & ST | Layer 3 & Layer 2 | `PRIORITY_1` (Spatial) / `PRIORITY_2` (SC) | Smaller cohort with paired spatial and single-cell data. |
| **HWANG_GSE202051** | GSE202051 | Hwang et al. (2022) | 43 | 43 | 0 | snRNA-seq & Visium | Layer 3 & Layer 2 | `PRIORITY_1` (Spatial) / `PRIORITY_2` (SC) | Large dual-modality cohort; includes neoadjuvant-treated samples. |

---

## 4. Planning & Validation File Reconciliations

We have updated and aligned all authoritative records:
1.  `00_admin/PROJECT_STATUS.md`: Updated to record Phase 9A.2 provenance corrections and cohort set expansion.
2.  `01_metadata/external_validation_dataset_inventory.tsv`: Rebuilt with the canonical dataset ID schema and corrected single-cell fields.
3.  `05_results/tables/phase9a_external_dataset_shortlist.tsv`: Updated to list the canonical single-cell IDs and correct patient counts.
4.  `01_metadata/external_validation_parameter_inventory.tsv`: Added distinct analysis parameters for the new Phase 9B2 cohorts (`VAL_SC_CELLULAR_SOURCE_LIN`, `VAL_SC_CELLULAR_SOURCE_MONCADA`, `VAL_SC_CELLULAR_SOURCE_HWANG`, `VAL_SPATIAL_CO_LOCAL_MONCADA`, `VAL_SPATIAL_CO_LOCAL_ZHANG`, etc.).
5.  `04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md`: Updated Section 8 to lock the expanded cohort set.
6.  `09_docs/methods/PDAC_external_validation_protocol.md`: Updated Section 3 to specify the multi-cohort single-cell validation strategy.
7.  `09_docs/references/phase9_external_validation_source_audit.tsv`: Corrected accessions and mappings (assigned CRA001160 to Peng and GSE111672 to Moncada).
8.  `09_docs/references/phase9_external_validation_sources.bib`: Updated citation keys and metadata.
9.  `09_docs/planning/DECISION_LOG.md`: Appended decision entry **D-32**.

---

## 5. Verification Results

We executed the following programmatic validators:
1.  **Provenance Consistency Validator (`06_scripts/python/15_validate_provenance_consistency.py`):** **`PASSED`**. Confirmed that:
    *   No planning file equates GSE111672 with Peng et al. (2019).
    *   No accession is mapped to the wrong publication or represents multiple distinct studies.
    *   Patient counts match exactly across all files.
    *   All cohorts included in Phase 9B2 have corresponding parameter-inventory rows.
    *   Authoritative cohort lists match perfectly across planning files.
2.  **General Manifest Validator (`06_scripts/python/00_validate_manifests.py`):** **`PASSED`**. Confirmed that all file manifest records are complete and checksums are verified.
3.  **Phase 9B2 Startup Script (`06_scripts/python/15_prepare_phase9b2_single_cell.py`):** **`PASSED`**. Successfully verified record agreement and wrote the block dataset inventory.
4.  **R Validation Script (`06_scripts/R/15_phase9b2_single_cell_validation.R`):** **`PASSED`**. Verified and exited successfully without performing biological analyses.

---

## 6. Final Readiness Decision

### **`READY_WITH_DATASET_LIMITATIONS`**

The single-cell and spatial validation records are now fully reconciled, consistent, and scientifically sound. Phase 9B2 may proceed with the expanded, multi-cohort validation set.

*Locked on: 2026-07-03*
*Reviewer: Antigravity*
