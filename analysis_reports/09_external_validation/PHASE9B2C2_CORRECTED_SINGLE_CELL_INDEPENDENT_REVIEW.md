# Phase 9B2C2 Independent Review of Corrected Single-Cell Cellular-Source Analysis

This document contains the independent statistical, implementation, annotation, and evidence review of the corrected Phase 9B2R primary single-cell cellular-source analysis. The audit was conducted using the locked criteria in [PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md](file://~/thesis/PDAC/04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md) and the standard protocol in [PDAC_external_validation_protocol.md](file://~/thesis/PDAC/09_docs/methods/PDAC_external_validation_protocol.md).

---

## 1. Executive Summary & Review Decision

*   **FINAL REVIEW DECISION:** **`PASS`**
*   **DECISION RATIONALE:** All critical, major, and minor implementation errors identified during Phase 9B2C have been resolved:
    1.  **FIND_01 (CRITICAL - Corrected):** The mandatory negative controls have been executed. Permutations (1,000 iterations for patient and cell-type labels), unrelated Hallmark pathways, and nonselected TF regulons are computed, yielding valid statistical distributions and empirical P-values. Size/expression-matched module controls are correctly documented as technically inapplicable because all target modules were excluded.
    2.  **FIND_02 (MAJOR - Corrected):** The locked 80% coverage threshold is enforced. All five transferred WGCNA modules (MEblack, MEblue, MEgreen, MEtan, MEgreenyellow) have coverages $< 49\%$ and are classified as `INSUFFICIENT_SINGLE_CELL_DATA`, excluding them from formal scoring, models, figures, and reporting.
    3.  **FIND_03 (MINOR - Corrected):** The hardcoded blanket `TO_VERIFY` categories have been removed. All TF evidence categories are programmatically derived using locked hierarchical rules based on target coverage, pseudobulk source models, malignant-axis associations, and composition sensitivity.
*   **SPATIAL VALIDATION PLANNING STATUS:** **`READY_TO_PROCEED`**. Phase 9B3 spatial-validation planning may begin immediately as no critical, major, or moderate findings remain.

---

## 2. Answers to the 13 Final Review Questions

1.  **Whether all Phase 9B2C findings are corrected:**
    *   **Yes**. All three findings (FIND_01, FIND_02, and FIND_03) have been corrected in Phase 9B2R.
2.  **Whether cohort provenance and counts are verified:**
    *   **Yes**. Canonical dataset `PENG_CRA001160` (accession CRA001160, BioProject PRJCA001063) is verified. It consists of 24 untreated primary PDAC tumor patients, 11 control pancreas donors, and exactly 57,530 analyzed cells. No GSE111672 aliasing was observed, and no raw files were downloaded.
3.  **Whether annotation and malignant-cell classification are verified:**
    *   **Yes**. The broad-class cell annotations and the conservative malignant classification (Ductal type 2 as `MALIGNANT`, Ductal type 1 as `AMBIGUOUS` in tumor patients, and control/acinar/endocrine epithelial as `NONMALIGNANT_EPITHELIAL`) are verified.
4.  **Whether patient-aware pseudobulk inference is verified:**
    *   **Yes**. The patient is the biological unit of replication. Pseudobulk profiles were aggregated requiring $\ge 20$ cells per patient-cell-type combination, and individual cells were not treated as independent observations.
5.  **Whether the five WGCNA modules are correctly excluded:**
    *   **Yes**. All five modules are correctly classified as `INSUFFICIENT_SINGLE_CELL_DATA` due to low coverage ($< 80\%$) and do not appear in formal models or figures.
6.  **Whether Hallmark scoring is verified:**
    *   **Yes**. ssGSEA via decoupleR was correctly run on the full Hallmark sets.
7.  **Whether TF activity scoring and classification are verified:**
    *   **Yes**. TF activities were computed via decoupleR/VIPER on DoRothEA A/B/C regulons, and classifications are programmatically derived. No expression proxies were used.
8.  **Whether all mandatory negative controls are complete:**
    *   **Yes**. Negative-control tables are complete: 64 rows are marked `EXECUTED` with nonmissing statistics, and 10 rows are marked `TECHNICALLY_INAPPLICABLE` due to module coverage exclusion.
9.  **Whether HALLMARK_PROTEIN_SECRETION is the sole malignant-cell-intrinsic supported feature:**
    *   **Yes**. It is the only feature classified as `MALIGNANT_CELL_INTRINSIC_SUPPORT`.
10. **Whether q = 0.0336 is independently reproduced:**
    *   **Yes**. The malignant-axis OLS HC3 regression for `HALLMARK_PROTEIN_SECRETION` yields a positive coefficient of 0.00400 (SE = 0.00155, p = 0.0168) and BH adjusted $q = 0.03361$, which matches the reanalysis report.
11. **The exact cellular-source category counts:**
    *   `CELL_COMPOSITION_EXPLAINED`: 20
    - `INSUFFICIENT_SINGLE_CELL_DATA`: 5
    - `STROMAL_OR_IMMUNE_SOURCE_SUPPORTED`: 4
    - `MALIGNANT_CELL_INTRINSIC_SUPPORT`: 1
    - `PARTIAL_CELLULAR_SUPPORT`: 1
    - `NOT_SUPPORTED_AT_CELLULAR_LEVEL`: 1
12. **Whether any CRITICAL, MAJOR, or MODERATE issue remains:**
    *   **No**. Zero critical, major, or moderate issues remain.
13. **Whether Phase 9B3 spatial-validation planning may begin:**
    *   **Yes**. Spatial validation planning is fully authorized.

---

## 3. Detailed Audit Findings and Verification Results

### Task 1: Complete Provenance and Input Audit
*   Official CNCB GSA processed count matrix (`count-matrix.txt`), published cell annotations (`all_celltype.txt`), and official checksums are verified.
*   The patient-level mapping (24 tumors, 11 controls) and total cells (57,530) are reproduced exactly. No duplicated cells or patient mappings exist.
*   Inclusion is strictly restricted to `PENG_CRA001160`; no FASTQ or BAM files were downloaded.

### Task 2: Correction Audit for FIND_01, FIND_02, and FIND_03
*   Verification details are recorded in [phase9b2c2_correction_verification.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_correction_verification.tsv). All findings are closed.

### Task 3: Cell Annotation Audit
*   Canonical cell-type markers (e.g., *PRSS1* for acinar, *MS4A1* for B cells, *EPCAM*/*KRT19* for ductal, *PECAM1* for endothelial, *COL1A1* for fibroblasts, *CD68* for myeloid, *CD3D* for T cells) were verified. Broad cell annotations are supported across multiple donors. Audit table: [phase9b2c2_annotation_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_annotation_audit.tsv).

### Task 4: Malignant-Cell Classification Audit
*   Tumor Ductal type 2 (11,315 cells) is confirmed as `MALIGNANT`. Ductal type 1 (3,117 cells in tumor patients) is confirmed as `AMBIGUOUS` and excluded from primary malignant models. Acinar, endocrine, and control ductal cells are confirmed as `NONMALIGNANT_EPITHELIAL`. Patient-level counts are verified in [phase9b2c2_malignant_cell_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_malignant_cell_audit.tsv).

### Task 5: Biological-Replicate and Pseudobulk Audit
*   Patient-cell-type pseudobulk counts (using $\ge 20$ cells as the cutoff) are verified. 121 combinations were excluded as ineligible. Patient-level models with patient blocking were used. Audit table: [phase9b2c2_pseudobulk_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_pseudobulk_audit.tsv).

### Task 6: Feature Inventory and Eligibility Audit
*   The eligibility of the 32 locked features is verified in [phase9b2c2_feature_eligibility_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_feature_eligibility_audit.tsv). The 5 WGCNA modules are ineligible, while the 2 Hallmark pathways and 25 TF regulons are eligible.

### Task 7: WGCNA Module Exclusion Audit
*   All five modules have single-cell coverages ranging from 25.2% to 48.5%, below the 80% threshold. They are correctly classified as `INSUFFICIENT_SINGLE_CELL_DATA` and do not appear in malignant-axis models or figures.

### Task 8: Moffitt and PurIST Scoring Audit
*   Moffitt50 basal-classical contrast scoring and Moffitt49 sensitivity scoring are verified. PurIST scoring was verified as ineligible due to single-cell data properties but not processed.

### Task 9: Hallmark Scoring Audit
*   ssGSEA GSVA scoring was performed correctly on the full hallmark sets. `HALLMARK_PROTEIN_SECRETION` (96 genes, coverage 97.9%) and `HALLMARK_SPERMATOGENESIS` (45 genes, coverage 82.2%) are eligible.

### Task 10: TF Activity Audit
*   DoRothEA VIPER scoring was run on eligible regulons with no expression proxies. No TF passed the axis q < 0.10 threshold. The null associations represent valid null results rather than coverage or implementation errors.

### Task 11: Malignant-Cell Axis Audit
*   OLS HC3 models for axis associations were audited. `HALLMARK_PROTEIN_SECRETION` is the sole axis-associated program ($q = 0.03361 < 0.10$, positive coefficient = 0.00400). All modules were successfully excluded. Verification table: [phase9b2c2_hallmark_tf_results_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_hallmark_tf_results_audit.tsv).

### Task 12: Cellular-Source Models Audit
*   Cellular source models were verified. TF localization results are correct: 19 composition-explained, 4 stromal/immune (ELF1, MBD2, ZBTB7A, ZNF384), 1 partial (ZNF740), 1 not-supported (KLF13).

### Task 13: Complete Evidence Classification Audit
*   Programmatic classifications match the locked rules: 20 composition-explained, 5 insufficient data, 4 stromal/immune, 1 malignant-intrinsic, 1 partial, 1 not-supported. Audit table: [phase9b2c2_evidence_category_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_evidence_category_audit.tsv).

### Task 14: Composition-Sensitivity Audit
*   Regressions on individual fractions (omitting control epithelial fraction) were verified. 20 features are composition-sensitive at $q < 0.10$. `HALLMARK_PROTEIN_SECRETION` is composition-sensitive (associated with malignant and lymphoid fractions) but remains supported as malignant-intrinsic because it is also axis-associated and localized to malignant cells. Audit table: [phase9b2c2_composition_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_composition_audit.tsv).

### Task 15: Negative-Control Audit
*   64 rows are marked `EXECUTED` with valid empirical null statistics. The 10 inapplicable rows are correctly marked `TECHNICALLY_INAPPLICABLE` due to module coverage exclusion. Audit table: [phase9b2c2_negative_control_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_negative_control_audit.tsv).

### Task 16: Multiple-Testing Audit
*   BH correction families were correctly separated. No nominal P-value was promoted.

### Task 17: Tumor-Control Comparison Audit
*   Control donors were used only for contextual comparison and did not redefine models or thresholds.

### Task 18: Heterogeneity Analysis Audit
*   Malignant cell heterogeneity is descriptive-only and cells were not treated as independent observations.

### Task 19: Figures and Report Audit
*   All `phase9b2r_` figures and reports were audited. They correctly represent the corrected evidence, exclude modules from formal support, and distinguish coverage from biological support.

### Task 20: Implementation and Validator Audit
*   The python validator script correctly checks the substance of the results (permutation iterations, coverage checks, negative control executions) and passes.

---

## 4. Audit Tables Summary

The 11 generated audit tables are:
1.  **Correction Verification Table:** [phase9b2c2_correction_verification.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_correction_verification.tsv)
2.  **Provenance QC Audit Table:** [phase9b2c2_provenance_qc_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_provenance_qc_audit.tsv)
3.  **Annotation Audit Table:** [phase9b2c2_annotation_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_annotation_audit.tsv)
4.  **Malignant Cell Audit Table:** [phase9b2c2_malignant_cell_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_malignant_cell_audit.tsv)
5.  **Pseudobulk Audit Table:** [phase9b2c2_pseudobulk_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_pseudobulk_audit.tsv)
6.  **Feature Eligibility Audit Table:** [phase9b2c2_feature_eligibility_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_feature_eligibility_audit.tsv)
7.  **Hallmark TF Results Audit Table:** [phase9b2c2_hallmark_tf_results_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_hallmark_tf_results_audit.tsv)
8.  **Composition Audit Table:** [phase9b2c2_composition_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_composition_audit.tsv)
9.  **Negative Control Audit Table:** [phase9b2c2_negative_control_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_negative_control_audit.tsv)
10. **Evidence Category Audit Table:** [phase9b2c2_evidence_category_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_evidence_category_audit.tsv)
11. **Review Findings Table:** [phase9b2c2_review_findings.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c2_review_findings.tsv)

*Date: 2026-07-03*
*Reviewer Agent: Antigravity*
