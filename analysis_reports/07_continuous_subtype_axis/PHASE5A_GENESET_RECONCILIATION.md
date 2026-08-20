# Phase 5A Gene-Set Reconciliation Report

## 1. Overview

This report reconciles the Moffitt gene-set definitions for the continuous basal–classical transcriptomic-axis analysis of the pancreatic ductal adenocarcinoma (PDAC) transcriptomics in the GSE172356 cohort. The reconciliation addresses and resolves the inconsistency between the 49 active Moffitt genes described in the Phase 5A method lock and the 50 genes in the verified reference signature file.

## 2. Reconciled Signature Files

Two separate, explicit signature files have been created in the reference directory `02_data/reference/PDAC_subtype_signatures/`:

### A. Primary Analysis Signature
* **Filename**: [Moffitt_50_gene_axis.tsv](file://~/thesis/PDAC/02_data/reference/PDAC_subtype_signatures/Moffitt_50_gene_axis.tsv)
* **Description**: Canonical verified Moffitt 50-gene signature with `LEMD1` retained.
* **Expected Row Count**: 50 active rows (plus header).
* **Expected Program Counts**:
  * **Basal-like**: 25 genes
  * **Classical**: 25 genes
* **SHA256 Checksum**: `3fa1790ff692898d01e2f4f8058d438c1263245c0d5316afab9840c968a2b72f`

### B. Sensitivity Analysis Signature (LEMD1-Excluded)
* **Filename**: [Moffitt_49_gene_axis_no_LEMD1.tsv](file://~/thesis/PDAC/02_data/reference/PDAC_subtype_signatures/Moffitt_49_gene_axis_no_LEMD1.tsv)
* **Description**: Prespecified sensitivity signature excluding `LEMD1` to assess impact of LEMD1 removal.
* **Expected Row Count**: 49 active rows (plus header).
* **Expected Program Counts**:
  * **Basal-like**: 24 genes
  * **Classical**: 25 genes
* **SHA256 Checksum**: `65cadb4c059a4b5be81efe03b8be1b5a6fc88937fd3eadf46f399ee007f1fc61`

---

## 3. Verification Findings

A programmatic validation script was executed to audit the properties of the signature files:

1. **One-Gene Difference**:
   * The files differ **only** by the presence of `LEMD1`.
   * `LEMD1` is present only in the Basal-like program of `Moffitt_50_gene_axis.tsv`.
   * No other gene symbol, mapping, program assignment, or column content differs between the two files.

2. **Expression-Matrix Coverage**:
   * All 50 mapped symbols in the primary signature and all 49 mapped symbols in the sensitivity signature were cross-referenced against the Phase 2B expression matrix `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz` (containing 42,654 genes).
   * **Result**: **100% Coverage**. Zero genes from either signature file are missing in the processed expression matrix.

---

## 4. Parameter Inventory Integration

The continuous axis parameters inventory at `01_metadata/continuous_axis_parameter_inventory.tsv` was updated to define separate and distinct analysis IDs:

* **`AXIS_MOFFITT50_PRIMARY`**: Row-scaled signature-mean contrast score using the 50-gene signature.
* **`AXIS_MOFFITT49_NO_LEMD1_SENSITIVITY`**: Row-scaled signature-mean contrast score using the 49-gene signature (LEMD1-excluded).

---

## 5. Phase 5B Readiness

* **Calculation Status**: Programmatic verification confirmed that no continuous scores, downstream models, or score-based files have been generated.
* **Recommendation**: With the signature files verified and parameter inventory locked, the cohort and analytical workflows are **fully ready** to proceed to the Phase 5B continuous scoring execution.
