# Phase 9A.3: Phase 9B2 Execution-Scope Correction Report

## 1. Executive Summary

This report documents the Phase 9A.3 execution-scope correction for the Layer 2 single-cell validation phase (Phase 9B2). In the Phase 9B2 restart, the execution was stopped correctly before data acquisition because of a conflict: the authoritative cohort table marked four datasets as active, while the approved execution contract specified **`PENG_CRA001160`** only.

This correction resolves the conflict by separating single-cell dataset suitability from execution authorization. We have updated the metadata, parameters, and startup validators to restrict the active primary execution set to `PENG_CRA001160` only. The other three cohorts (`LIN_GSE154778`, `MONCADA_GSE111672`, and `HWANG_GSE202051`) remain validated and suitable but are designated as supplementary, requiring subsequent separate authorization for execution.

No biological data acquisition, scoring, or downstream analyses have been performed. All validators have been successfully run in validation-only mode and pass.

**Final Readiness Status:** **`READY_FOR_PHASE9B2_PRIMARY_EXECUTION`**

---

## 2. Inconsistency Analysis & Scope Separation

The previous Phase 9A.2 correction expanded the single-cell validation cohort set to four datasets. However, the primary execution contract requires staging these validations to manage complexity, focusing first on the primary dataset.

We have resolved this by replacing the single, ambiguous `included_in_phase9b2` field in the authoritative cohort set with explicit, fine-grained fields:
1. `layer2_suitability_status`
2. `layer2_analysis_role`
3. `included_in_phase9b2_primary`
4. `included_in_phase9b2_supplementary`
5. `supplementary_execution_status`
6. `layer3_spatial_role`

This allows us to clearly distinguish:
- Datasets suitable for validation.
- Datasets authorized for the *current* primary execution (`included_in_phase9b2_primary == TRUE`).
- Datasets reserved for subsequent supplementary executions (`included_in_phase9b2_supplementary == TRUE` and `supplementary_execution_status == NOT_YET_AUTHORIZED`).

> [!IMPORTANT]
> This sequencing decision is purely operational to manage complexity and prioritize initial validation efforts. It does not imply that the supplementary datasets are scientifically or biologically unsuitable.

---

## 3. Updated Authoritative Execution Scope

The final resolved execution scope for the single-cell cohorts is documented in [phase9a3_phase9b2_execution_scope.tsv](file://~/thesis/PDAC/05_results/tables/phase9a3_phase9b2_execution_scope.tsv):

| Canonical Dataset ID | Accession | Layer 2 Suitability | Layer 2 Analysis Role | Included in Primary | Included in Supplementary | Supplementary Status | Layer 3 Spatial Role | Current Execution Authorized |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **PENG_CRA001160** | CRA001160 | PRIORITY_1 | PRIMARY | TRUE | FALSE | NOT_APPLICABLE | NOT_APPLICABLE | **TRUE** |
| **LIN_GSE154778** | GSE154778 | PRIORITY_2 | SUPPLEMENTARY | FALSE | TRUE | NOT_YET_AUTHORIZED | NOT_APPLICABLE | **FALSE** |
| **MONCADA_GSE111672** | GSE111672 | PRIORITY_2 | EXPLORATORY | FALSE | TRUE | NOT_YET_AUTHORIZED | PRIORITY_1 | **FALSE** |
| **HWANG_GSE202051** | GSE202051 | PRIORITY_2 | TREATMENT_SENSITIVITY | FALSE | TRUE | NOT_YET_AUTHORIZED | PRIORITY_1 | **FALSE** |

---

## 4. Parameter Inventory Alignment

The parameter inventory [external_validation_parameter_inventory.tsv](file://~/thesis/PDAC/01_metadata/external_validation_parameter_inventory.tsv) has been updated to align the execution status columns:
- `VAL_SC_CELLULAR_SOURCE_PENG` is marked as **`ACTIVE_PRIMARY`**.
- `VAL_SC_CELLULAR_SOURCE_LIN` is marked as **`PLANNED_SUPPLEMENTARY`**.
- `VAL_SC_CELLULAR_SOURCE_MONCADA` is marked as **`PLANNED_EXPLORATORY`**.
- `VAL_SC_CELLULAR_SOURCE_HWANG` is marked as **`PLANNED_TREATMENT_SENSITIVITY`**.
- All non-Phase 9B2 rows (bulk host, spatial, and microbiome parameters) are marked as **`NOT_APPLICABLE`** for the single-cell validation execution.

---

## 5. Audit Trail of Stopped Phase 9B2 Attempts

We preserve the audit trail of both previous stopped Phase 9B2 attempts to maintain complete traceability:

1. **First Attempt (Phase 9B2 Initial Startup - 2026-07-03):**
   - **Trigger:** Disagreement of single-cell cohort definitions across full inventory, shortlist, and parameter tables.
   - **Action:** Startup validator returned exit code 2. The run stopped correctly before data acquisition. No files were downloaded.
   - **Outputs:** [phase9b2_single_cell_dataset_inventory.tsv](file://~/thesis/PDAC/01_metadata/phase9b2_single_cell_dataset_inventory.tsv) (blocked placeholder state) and [PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md](file://~/thesis/PDAC/04_analysis/09_external_validation/PHASE9B2_SINGLE_CELL_CELLULAR_SOURCE_RESULTS.md) (documented `STOPPED_BEFORE_DATA_ACQUISITION`).
2. **Second Attempt (Phase 9B2 Restart - 2026-07-03):**
   - **Trigger:** Conflict between the PENG-only primary execution contract and the authoritative Phase 9A.2 cohort set, which marked all four single-cell datasets as `included_in_phase9b2=yes`.
   - **Action:** Restart-specific startup validator reported `FAIL` for the execution set check and stopped before data acquisition. No files were downloaded.
   - **Outputs:** [phase9b2_restart_runtime_validation.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2_restart_runtime_validation.tsv) (recorded `FAIL` on `phase9b2_included_set` and `primary_execution_set`) and updated stopped results report.

---

## 6. Programmatic Verification Results

We executed the complete validator suite to verify the consistency of the Phase 9A.3 correction:

1. **Provenance Consistency Validator (`06_scripts/python/15_validate_provenance_consistency.py`):**
   - **Command:** `python3 06_scripts/python/15_validate_provenance_consistency.py`
   - **Result:** `PASSED`
   - **Output:**
     ```
     Checking database-lookup, citation-management, and experimental-design rules...
     Provenance consistency validation PASSED successfully.
     ```
2. **Preparation Script in Validation-Only Mode (`06_scripts/python/15_prepare_phase9b2_single_cell.py`):**
   - **Command:** `python3 06_scripts/python/15_prepare_phase9b2_single_cell.py --validation-only`
   - **Result:** `READY_FOR_DATA_ACQUISITION`
   - **Output:**
     ```
     expected primary set: PENG_CRA001160
     observed primary set: PENG_CRA001160
     supplementary planned datasets: HWANG_GSE202051, LIN_GSE154778, MONCADA_GSE111672
     result: READY_FOR_DATA_ACQUISITION
     ```
3. **General Manifest Validator (`06_scripts/python/00_validate_manifests.py`):**
   - **Command:** `python3 06_scripts/python/00_validate_manifests.py`
   - **Result:** `PASSED`
   - **Output:**
     ```
     Manifest validation passed for 3 files and 524 data rows.
     ```
4. **Phase 9B2 Startup Validator (`06_scripts/python/15_validate_phase9b2_single_cell.py`):**
   - **Command:** `python3 06_scripts/python/15_validate_phase9b2_single_cell.py`
   - **Result:** `PASSED`
   - **Output:**
     ```
     Phase 9B2 corrected startup validation passed successfully.
     observed primary set: PENG_CRA001160
     result: READY_FOR_PHASE9B2_PRIMARY_EXECUTION
     ```

---

## 7. Final Decision

Based on the complete reconciliation of metadata, parameter tables, and validator scripts:

### **`READY_FOR_PHASE9B2_PRIMARY_EXECUTION`**

Phase 9B2 is ready to proceed with data acquisition and cellular source validation restricted to **`PENG_CRA001160`** only.

*Report compiled on: 2026-07-03*
*Agent: Antigravity*
