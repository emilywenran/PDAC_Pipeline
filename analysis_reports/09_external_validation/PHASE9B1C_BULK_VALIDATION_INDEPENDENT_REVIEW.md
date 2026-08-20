# Phase 9B1C: Independent Statistical, Implementation, and Evidence Review

This document contains the independent review of the Phase 9B1 bulk-transcriptome external validation. The audit was conducted using the locked criteria in [PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md](file://~/thesis/PDAC/04_analysis/09_external_validation/PHASE9A_EXTERNAL_VALIDATION_METHOD_LOCK.md) and the standard protocol in [PDAC_external_validation_protocol.md](file://~/thesis/PDAC/09_docs/methods/PDAC_external_validation_protocol.md).

---

## 1. Executive Summary & Review Decision

*   **FINAL REVIEW DECISION:** **`FAIL_REQUIRES_REANALYSIS`**
*   **DECISION RATIONALE:** The Phase 9B1 validation implementation contains several critical and major implementation errors that violate the prospectively locked Phase 9A method plan and parameter inventory. Specifically:
    1.  **PurIST Intercept Omitted:** The PurIST classifier omitted the model intercept ($-6.815$) in the logistic link function, rendering all calculated probabilities invalid ($\ge 0.5$ for all samples) and skewing predictions.
    2.  **Missing-Gene Policy Violations:** WGCNA module eigengene-like scores were computed on microarray datasets (GSE71729 and GSE62452) where gene coverage was between 16% and 40% (violating the $80\%$ minimum threshold). These cohorts should have been excluded under the missing-gene policy, which downgrades the module evidence classifications from `EXTERNALLY_REPLICATED` to `PARTIALLY_REPLICATED` (replicated in TCGA_PAAD only). Additionally, the PurIST score for GSE62452 was incorrectly set to NaN even though coverage was $87.5\% \ge 80\%$.
    3.  **Pathway Scoring Implementation Failure:** Instead of running the locked ssGSEA algorithm via decoupleR on the full MSigDB Hallmark gene sets (96 and 135 genes, respectively), the script used a simple rank-percentile mean over a custom 15-gene proxy list.
    4.  **TF Activity Proxy Error:** Single-gene expression was used as a proxy for TF activity without running decoupleR/VIPER, and the results were inappropriately classified as `NOT_REPLICATED` or `PARTIALLY_REPLICATED` instead of remaining `TO_VERIFY`.
    5.  **Negative Control Omission:** Unrelated pathway control statistics were hardcoded to `NaN` instead of being calculated.
*   **COHORTS AUDITED AND VERIFIED:** Yes, the independence of the three cohorts (TCGA_PAAD: 178 tumor samples; GSE71729: 145 tumor samples; GSE62452: 69 tumor samples) from GSE172356 (discovery) and from each other has been verified.
*   **HALLMARK NULL RESULT VERIFIED:** **`YES`** (Under both the proxy scoring and correct rules, the Hallmark pathway association fails replication due to opposite effect direction for `HALLMARK_SPERMATOGENESIS` and lack of significance in microarrays for `HALLMARK_PROTEIN_SECRETION`).
*   **MODULE REPLICATION STATUS:**
    *   **`MEblack`**: `PARTIALLY_REPLICATED_HOST_FEATURE` (replicated in 1 cohort: TCGA_PAAD; GSE71729 and GSE62452 failed coverage).
    *   **`MEblue`**: `PARTIALLY_REPLICATED_HOST_FEATURE` (replicated in 1 cohort: TCGA_PAAD; GSE71729 and GSE62452 failed coverage).
    *   **`MEgreen`**: `PARTIALLY_REPLICATED_HOST_FEATURE` (replicated in 1 cohort: TCGA_PAAD; GSE71729 and GSE62452 failed coverage).
    *   **`MEred`**: `NOT_REPLICATED` (fails replication in TCGA_PAAD; GSE71729 and GSE62452 failed coverage).
*   **TF VALIDATION STATUS:** Remains **`TO_VERIFY`** and must be deferred until decoupleR/VIPER can be executed in the active project-local `renv` environment.
*   **PHASE 9B2 SINGLE-CELL VALIDATION STATUS:** **`DEFERRED`** until Phase 9B1 re-analysis is executed and passed.

---

## 2. Detailed Task-by-Task Audit Findings

### Task 1: Cohort and Input Audit
*   **Accessions & Sources:** TCGA-PAAD (STAR counts / UCSC Xena), GSE71729 (NCBI GEO Custom Agilent), GSE62452 (NCBI GEO Affymetrix GPL6244).
*   **Sample Counts:** Retained primary tumor counts are exactly 178 (TCGA), 145 (GSE71729), and 69 (GSE62452).
*   **Sample Exclusions:** TCGA normal samples (barcode 11) and GSE normal tissue control samples were correctly excluded.
*   **Gene-expression Scale:** TCGA is log2(FPKM-UQ+1); GSE71729 is log-ratio; GSE62452 is log2 RMA intensity.
*   **Duplicates & Mapping:** Row duplicates were resolved correctly using the highest mean probe. Probe-to-gene symbols mapped via official SOFT platform files.
*   **Cohort Independence:** Verified that there is no patient overlap between GSE172356, TCGA-PAAD, GSE71729, and GSE62452.
*   **Missingness & Provenance:** Written metadata and acquisition manifests (`phase9b1_bulk_data_acquisition_manifest.tsv`) match the downloaded records and SHA256 hashes.

### Task 2: Basal–Classical Scoring Audit
*   **Moffitt50 & Moffitt49:** Centroid-anchored scoring was correctly implemented using the z-scored expression matrices. Note: the script used row-mean z-scoring instead of the row-median centering specified in the parameters.
*   **PurIST Implementation Errors:**
    *   **Critical Error (FIND_01):** The model intercept of $-6.815$ was omitted, causing all log-odds to be non-negative and all probabilities to be $\ge 0.5$ (mean probability 0.74 across samples), invalidating the PurIST predictions.
    *   **Major Error (FIND_02):** In GSE62452, PurIST coverage was 7/8 pairs (87.5%), which is above the 80% threshold. The script returned all NaNs instead of calculating and rescaling, violating the missing-gene policy.

### Task 3: Hallmark Validation Audit
*   **Major Error (FIND_04):** The script substituted the locked decoupleR ssGSEA method (minsize=15) with a simple average rank-percentile score over a custom 15-gene proxy list for `HALLMARK_PROTEIN_SECRETION` and `HALLMARK_SPERMATOGENESIS`.
*   **Null Result Verification:**
    *   `HALLMARK_PROTEIN_SECRETION` (expected negative association) only reached significance in TCGA_PAAD ($\beta = -0.0013, q = 0.039$) but was non-significant in GSE62452 ($q = 0.8169$) and showed opposite direction in GSE71729.
    *   `HALLMARK_SPERMATOGENESIS` showed a statistically significant positive association across all cohorts (opposite to the discovery negative association). It is classified as `NOT_REPLICATED`.
*   The null result for pathways reflects the locked replication rules (opposite sign and lack of consistent significance across multiple cohorts) rather than just proxy scoring, but the scoring itself remains invalid.

### Task 4: Transferred WGCNA Module Audit
*   **Major Error (FIND_03):** WGCNA module genes contain a large fraction of non-coding RNA genes that are absent from microarray platforms. Gene coverage for GSE71729 and GSE62452 was only 16%-40%, violating the locked 80% minimum coverage threshold. Scores should not have been calculated or pooled.
*   **Evidence Category Downgrading:** Enforcing the 80% coverage threshold leaves only TCGA_PAAD as a valid validation cohort. Consequently, no module can be classified as `EXTERNALLY_REPLICATED_HOST_FEATURE` (requires consistent replication in $\ge 2$ cohorts).
    *   `MEblack`, `MEblue`, and `MEgreen` are downgraded to `PARTIALLY_REPLICATED_HOST_FEATURE` (replicated in TCGA_PAAD only).
    *   `MEred`, `MEgreenyellow`, `MEpurple` are `NOT_REPLICATED` (failed to replicate in TCGA_PAAD).

### Task 5: TF Validation Status Audit
*   **Major Error (FIND_05):** Single-gene expression was used as a proxy for TF activity instead of running the VIPER algorithm via decoupleR. In the evidence table, TFs were incorrectly classified as `NOT_REPLICATED` or `PARTIALLY_REPLICATED` instead of remaining `TO_VERIFY`.
*   **Recommendation:** TF activity validation must remain `TO_VERIFY` and be deferred. Substituting proxy gene expression is a direct protocol violation.

### Task 6: Cross-Cohort Synthesis Audit
*   Platform-incompatible WGCNA scores (calculated on 20% vs. 95% signature genes) were inappropriately pooled in the random-effects meta-analysis (`phase9b1_cross_cohort_synthesis.tsv`). Enforcing the coverage threshold renders meta-analysis inapplicable due to having only one valid cohort.

### Task 7: Negative-Control Audit
*   **Major Error (FIND_06):** Unrelated pathway negative controls (5 Hallmark pathways) were set to `NaN` in the script rather than being calculated, leaving the negative-control audit incomplete. Patient-label permutation and size-matched random modules were completed successfully.

### Task 8: Evidence-Category Verification
The verified reviewer classifications (correcting for coverage filters and TF proxies) are stored in [phase9b1c_host_feature_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b1c_host_feature_audit.tsv). Key changes include:
*   Downgrading `MEblack`, `MEblue`, and `MEgreen` from `EXTERNALLY_REPLICATED` to `PARTIALLY_REPLICATED_HOST_FEATURE` (replicated in TCGA_PAAD only).
*   Reclassifying all 34 TFs and `HALLMARK_PROTEIN_SECRETION` to `TO_VERIFY` due to implementation/scoring failures.

---

## 3. Audit Tables Summary

*   **Review Findings:** [phase9b1c_review_findings.tsv](file://~/thesis/PDAC/05_results/tables/phase9b1c_review_findings.tsv) (Lists 6 implementation findings: 1 Critical, 4 Major, 1 Moderate).
*   **Module Replication Audit:** [phase9b1c_module_replication_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b1c_module_replication_audit.tsv) (Details module coverage and statistics with and without coverage filtering).
*   **Host Feature Audit:** [phase9b1c_host_feature_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase9b1c_host_feature_audit.tsv) (Contains verified replication statistics and reviewer classifications for all 43 host features).

---

## 4. Final Review Decision and Answers

1.  **Are TCGA_PAAD, GSE71729, and GSE62452 verified?** Yes, their patient sample independence is verified.
2.  **Is the Hallmark null result verified?** Yes, both Hallmark pathways fail to meet external replication criteria.
3.  **Are MEblack, MEblue, MEgreen, and MEred externally replicated?** No. Under correct coverage filters ($\ge 80\%$), they fail multi-cohort replication. `MEblack`, `MEblue`, and `MEgreen` are replicated in TCGA_PAAD only (`PARTIALLY_REPLICATED_HOST_FEATURE`), and `MEred` is not replicated (`NOT_REPLICATED`).
4.  **How many cohorts support each module?**
    *   `MEblack`: 1 cohort (TCGA_PAAD)
    *   `MEblue`: 1 cohort (TCGA_PAAD)
    *   `MEgreen`: 1 cohort (TCGA_PAAD)
    *   `MEred`: 0 cohorts (non-significant in TCGA_PAAD)
    *   `MEgreenyellow`: 0 cohorts
    *   `MEpurple`: 0 cohorts
    *   `MEtan`: 1 cohort (TCGA_PAAD)
5.  **Does TF validation remain TO_VERIFY?** Yes, TF validation remains `TO_VERIFY` and is deferred until decoupleR/VIPER can be executed.
6.  **May Phase 9B2 single-cell validation proceed?** **`NO`**. Phase 9B2 is deferred until the Phase 9B1 bulk-host validation is re-analyzed and corrected.

*Date: 2026-07-03*
*Reviewer Agent: Antigravity*
