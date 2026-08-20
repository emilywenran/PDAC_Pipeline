# Phase 9B2C Independent Review of single-cell cellular-source analysis

This document contains the independent statistical, implementation, annotation, and evidence review of the completed Phase 9B2 primary single-cell cellular-source analysis. The audit was conducted using the locked criteria in [PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md](file://~/thesis/PDAC/04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md) and the standard protocol in [PDAC_external_validation_protocol.md](file://~/thesis/PDAC/09_docs/methods/PDAC_external_validation_protocol.md).

---

## 1. Executive Summary & Review Decision

*   **FINAL REVIEW DECISION:** **`FAIL_REQUIRES_REANALYSIS`**
*   **DECISION RATIONALE:** The independent review identified two significant execution failures:
    1.  **CRITICAL: Negative Controls Not Actually Executed (FIND_01):** The primary execution script did not execute the mandatory negative control and falsification analyses (patient-label permutations, cell-type-label permutations, expression-matched randomized modules, and unrelated Hallmark scoring). Instead, it generated hardcoded placeholder rows in the results table `phase9b2_negative_control_results.tsv` and reported them as `PASS_DESCRIPTIVE_CONTROL_GENERATED` or `TO_VERIFY`. This directly violates the prospective method lock which requires negative controls.
    2.  **MAJOR: Coverage Threshold Violated for WGCNA Modules (FIND_02):** The R script violated the prospective 80% coverage threshold. All five co-expression modules had coverages $< 49\%$ (ranging from 25.1% to 48.5%) but were still scored, evaluated in association models, and reported as having `MALIGNANT_CELL_INTRINSIC_SUPPORT` or `STROMAL_OR_IMMUNE_SOURCE_SUPPORTED` instead of being excluded and classified as `INSUFFICIENT_SINGLE_CELL_DATA`.
*   **SPATIAL VALIDATION PLANNING STATUS:** **`BLOCKED`**. Spatial validation planning may NOT begin because critical and major implementation errors remain unresolved, and required negative controls are incomplete.

---

## 2. Answers to the 12 Final Review Questions

1.  **Whether PENG_CRA001160 provenance is verified:**
    *   **Yes**. The cohort identity is verified: canonical_dataset_id `PENG_CRA001160`, accession `CRA001160`, BioProject `PRJCA001063`, Peng et al. (2019). Processed files `count-matrix.txt`, `all_celltype.txt`, and `md5sum.txt` match the CNCB GSA records. No GSE111672 aliasing was observed, and no raw files were downloaded.
2.  **Whether 57,530 cells and 35 patients/donors are verified:**
    *   **Yes**. The analysis processed exactly 57,530 cells across 35 individuals (24 untreated primary PDAC tumor patients, 11 control pancreas donors). All patient and cell-level counts are reproduced and verified.
3.  **Whether malignant-cell classification is verified:**
    *   **Yes**. Tumor `Ductal type 2` cells were correctly classified as `MALIGNANT` (11,315 cells). Ductal type 1 cells were classified as `AMBIGUOUS` (3,117 cells in tumor patients) and were not silently counted as malignant. Acinar/endocrine and control epithelial cells were classified as `NONMALIGNANT_EPITHELIAL`. Non-epithelial cells were classified as `NOT_APPLICABLE`.
4.  **Whether patient-aware pseudobulk inference is verified:**
    *   **Yes**. The patient was treated as the biological replicate. The cell-by-gene matrix was aggregated into 225 eligible patient-cell-type pseudobulk columns (requiring $\ge 20$ cells per combination). Cells were not treated as independent observations in the models.
5.  **Whether Hallmark and WGCNA module scoring is verified:**
    *   **Partial**. Hallmark pathway scoring using decoupleR ssGSEA is verified. However, WGCNA module scoring is **FAILED** because all five modules had coverages far below the 80% threshold (ranging from 25.1% to 48.5%) and should have been excluded.
6.  **Whether TF regulon activity scoring is verified:**
    *   **Yes**. TF activities were correctly computed at the pseudobulk level using VIPER and DoRothEA A/B/C regulons. No TF-symbol gene expression proxy was used, and regulons with $< 15$ targets or $< 80\%$ coverage were excluded.
7.  **Whether the three malignant-cell-intrinsic features are verified:**
    *   **No**. Only `HALLMARK_PROTEIN_SECRETION` has verified malignant-intrinsic support. The co-expression modules `MEblue` and `MEgreen` failed the coverage threshold and must be classified as `INSUFFICIENT_SINGLE_CELL_DATA`.
8.  **Whether the five malignant-axis associations are verified:**
    *   **No**. Only `HALLMARK_PROTEIN_SECRETION` has a verified positive association with the malignant Moffitt50 axis (q = 0.0336 < 0.10). The other four modules (`MEblack`, `MEblue`, `MEgreen`, `MEtan`) failed coverage and should have been excluded from association models.
9.  **Which features are composition-explained:**
    *   The features verified as composition-explained (significantly associated with cell fractions at q < 0.10) are: `MEgreenyellow` (ineligible due to coverage), `CTCFL`, `JUNB`, `KLF9`, `TFAP4`, `TP63`, `GRHL2`, `KLF1`, `SIX5`, and `SNAPC4`.
10. **Whether required negative controls are complete:**
    *   **No**. The negative controls are incomplete. None of the five required controls were actually executed; the script wrote placeholder results instead.
11. **Whether Phase 9B2 evidence classifications are verified:**
    *   **No**. Evidence classifications are invalidated due to the WGCNA coverage violations and hardcoded negative controls.
12. **Whether spatial-validation planning may begin:**
    *   **No**. Spatial validation planning is blocked until a re-analysis completes all negative controls and enforces the coverage thresholds.

---

## 3. Detailed Task-by-Task Audit Findings

### Task 1: Data Provenance and Cohort Audit
*   Official CNCB GSA processed files were used: `count-matrix.txt`, `all_celltype.txt`, and `md5sum.txt`.
*   MD5 checksums and sizes are correct:
    *   `count-matrix.txt`: 2,771,872,913 bytes
    *   `all_celltype.txt`: 2,101,436 bytes
*   Patient mapping and counts are correct: 24 tumors, 11 controls.
*   No duplicated cells or patients are present.
*   No raw sequence files (FASTQ, BAM, SRA) were downloaded.
*   Inference scope was preserved (restricted to PENG_CRA001160).
*   No GSE111672 aliasing was observed.

### Task 2: Expression Object and QC Audit
*   Expression value type is verified as raw processed count counts, correctly converted to log2-CPM for pseudobulk profile scoring.
*   Gene symbol duplicate handling (Probe aggregation and highest mean expression) is verified.
*   Processed-source QC was respected; no arbitrary cell-filtering thresholds were added.
*   Donor and batch structure are complete and patient assignments are complete.

### Task 3: Major Cell-Type Annotation Audit
*   Broad cell type annotations were verified. The reviewed broad classes are: acinar epithelial (1,935 cells), B cell (2,447 cells), nonmalignant ductal epithelial (10,317 cells), malignant ductal/epithelial (11,315 cells), endocrine epithelial (729 cells), endothelial (9,117 cells), fibroblast/CAF (12,649 cells, combining fibroblast and stellate cells), and myeloid/macrophage (5,361 cells).
*   NK cells and mast cells are absent in the source dataset.
*   Supporting markers are canonical and verified. Annotation status is PASS. All reviewed annotations are supported across multiple patients and not driven by a single patient or single marker.

### Task 4: Malignant-Cell Classification Audit
*   Tumor `Ductal type 2` cells were classified as `MALIGNANT`.
*   Tumor `Ductal type 1` cells were retained as `AMBIGUOUS`.
*   Control epithelial and non-malignant epithelial (acinar, endocrine) cells were classified as `NONMALIGNANT_EPITHELIAL`.
*   No ambiguous epithelial cells were silently counted as malignant. Patient-level counts are recorded in [phase9b2c_malignant_cell_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_malignant_cell_audit.tsv).

### Task 5: Biological-Replicate Audit
*   The patient was correctly treated as the biological replicate.
*   Cell-level observations were not treated as independent replicates.
*   Pre-specified threshold of $\ge 20$ cells per patient-cell-type combination was enforced.
*   121 ineligible patient-cell-type combinations were excluded and documented in the pseudobulk inventory.
*   Patient-level fixed effects were used in cellular source models to approximate repeated-measures.
*   Pseudobulk counts were successfully verified and reproduced in the audit.

### Task 6: Feature-List and Coverage Audit
*   The 32 locked features were present in the input.
*   No post hoc feature additions occurred.
*   WGCNA modules blue and green failed the 80% coverage threshold. Enrichment models and classifications on WGCNA modules represent a major implementation violation because these features should have been excluded.
*   All TF target coverages were checked and low-coverage regulons were correctly excluded.

### Task 7: Host-State Scoring Audit
*   Moffitt50 basal and classical centroids were scored on log2-CPM profiles using mean z-scores.
*   The contrast direction (basal - classical) matches discovery.
*   The Moffitt49 no-LEMD1 sensitivity score was calculated correctly.
*   PurIST implementation is verified: includes beta0 intercept, all 8 gene pairs, and logistic transform.
*   Moffitt50 patient scores were successfully verified.

### Task 8: Hallmark and Module Scoring Audit
*   ssGSEA via decoupleR was run on the full Hallmark sets.
*   No reduced proxy gene set was used.
*   WGCNA modules MEblack, MEblue, MEgreen, MEtan, and MEgreenyellow were scored on rank matrices using discovery gene lists. No new single-cell WGCNA network was built.
*   Scoring was correct, but their reporting violates the coverage rule.

### Task 9: TF Activity Audit
*   DoRothEA A/B/C regulons were scored using decoupleR VIPER.
*   No TF expression proxy was used.
*   Regulon target coverage was checked.
*   Reported result that no TF regulon passed the malignant axis threshold (q < 0.10) is verified.
*   This represents a true null result rather than a power or implementation failure, as the statistical tests were correctly powered at the pseudobulk level but failed to show strong correlation with Moffitt50.

### Task 10: Cellular-Source Classification Audit
*   Single-cell evidence categories reported in `phase9b2_cellular_source_evidence.tsv` were audited.
*   The three features reported as malignant-intrinsic were `HALLMARK_PROTEIN_SECRETION`, `MEblue`, and `MEgreen`.
*   Since `MEblue` and `MEgreen` failed coverage, their classification must be updated to `INSUFFICIENT_SINGLE_CELL_DATA`.
*   Only `HALLMARK_PROTEIN_SECRETION` remains validated as malignant-cell intrinsic.

### Task 11: Malignant-Cell Axis Association Audit
*   Five features were reported as associated with the Moffitt50 contrast in malignant cells at q < 0.10: `HALLMARK_PROTEIN_SECRETION` (positive), `MEblack` (positive), `MEblue` (positive), `MEgreen` (negative), and `MEtan` (negative).
*   Since the WGCNA modules failed coverage, their associations are invalid. Only `HALLMARK_PROTEIN_SECRETION` is a verified association.
*   The OLS models were fit correctly using robust standard errors, but the threshold of q < 0.10 was not prospectively locked in Phase 9A, which represents a minor design/reporting issue.

### Task 12: Composition-Sensitivity Audit
*   Cell-fraction sensitivity models were evaluated. Predictors included malignant, endothelial, myeloid, fibroblast/CAF, and lymphoid fractions.
*   Multicollinearity exists among fractions because they sum to 1. Regressions were run on single fractions separately, which avoids multicollinearity but requires care in interpretation.
*   The reference/omitted cell type is `nonmalignant_epithelial`.
*   Ten features are classified as `CELL_COMPOSITION_EXPLAINED`.

### Task 13: Malignant-State Heterogeneity Audit
*   Within-patient malignant cell heterogeneity was analyzed descriptively. Cells were not treated as independent patients.
*   No new Hybrid cutoff was optimized.

### Task 14: Tumor-Control Comparison Audit
*   The 11 control pancreases were used only for contextual comparison.
*   Control samples did not redefine feature lists or evidence thresholds.

### Task 15: Negative-Control Audit
*   All negative controls are incomplete/unresolved. The R script wrote placeholder rows and did not run the simulations.

### Task 16: Reporting-Language Audit
*   The report did not claim Ochrobactrum was tested, nor did it infer microbial localization or causality.
*   All reporting language remains associative.

---

## 4. Audit Tables Summary

The six generated audit tables are:
1.  **Annotation Audit Table:** [phase9b2c_annotation_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_annotation_audit.tsv)
2.  **Malignant Cell Audit Table:** [phase9b2c_malignant_cell_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_malignant_cell_audit.tsv)
3.  **Pseudobulk Audit Table:** [phase9b2c_pseudobulk_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_pseudobulk_audit.tsv)
4.  **Feature Evidence Audit Table:** [phase9b2c_feature_evidence_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_feature_evidence_audit.tsv)
5.  **Negative Control Audit Table:** [phase9b2c_negative_control_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_negative_control_audit.tsv)
6.  **Review Findings Table:** [phase9b2c_review_findings.tsv](file://~/thesis/PDAC/05_results/tables/phase9b2c_review_findings.tsv)

*Date: 2026-07-03*
*Reviewer Agent: Antigravity*
