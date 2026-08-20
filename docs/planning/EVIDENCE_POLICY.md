# Evidence Policy

This document defines the evidence tiers utilized to grade the strength of all scientific findings, statistical associations, and biomarker candidates identified in this project.

---

## 1. Overview of Evidence Tiers

All biological findings and statistical associations must be classified into one of the following four tiers. No claim may be presented in a manuscript draft or presentation without its corresponding evidence tier assignment.

```
┌────────────────────────────────────────────────────────┐
│  Tier A: Replicated in Independent Cohorts             │
├────────────────────────────────────────────────────────┤
│  Tier B: Internally Robust under Resampling            │
├────────────────────────────────────────────────────────┤
│  Tier C: Exploratory Association                       │
├────────────────────────────────────────────────────────┤
│  Tier D: Literature-Supported Hypothesis Only          │
└────────────────────────────────────────────────────────┘
```

---

## 2. Detailed Tier Definitions and Criteria

### Tier A: Replicated in Independent Cohorts
*   **Definition:** Findings that are statistically significant in the primary cohort and successfully replicated in an independent external dataset with consistent effect directions.
*   **Application:** Primarily applied to host transcriptomic signatures, continuous subtype scoring behavior, and cell-type composition changes (e.g., validated in TCGA-PAAD). Any TCGA-derived microbiome signal must be treated as exploratory due to extensive contamination and batch concerns, and we do not assume that an identical external paired transcriptome-microbiome cohort exists.
*   **Mandatory Criteria:**
    1.  **Primary Dataset (GSE172356 / PRJNA719915):** Statistical significance achieved at Benjamini-Hochberg False Discovery Rate (FDR) $q < 0.05$.
    2.  **External Validation Dataset (e.g., TCGA-PAAD, ICGC, scRNA-seq):** Independent statistical significance achieved at $p < 0.05$ (or FDR $q < 0.05$ where appropriate).
    3.  **Directional Consistency:** The effect size must have the same sign (e.g., positive correlation in both, or overexpression in the same subtype).
    4.  **Confounder Control:** The association must remain statistically significant after adjusting for tumor purity (ESTIMATE score), immune/stromal infiltration, and technical covariates.

### Tier B: Internally Robust under Resampling and Sensitivity Analysis
*   **Definition:** Findings that are statistically significant in the primary dataset and demonstrate high stability under resampling, cross-validation, and sensitivity analysis, but lack independent external cohort replication.
*   **Application:** Subtype classification stability, internal cross-validation of classifiers, and primary cohort host-microbiome associations.
*   **Mandatory Criteria:**
    1.  **Primary Significance:** FDR $q < 0.05$ in the main unresampled model.
    2.  **Sensitivity Testing:** The association or classification must remain statistically significant ($p < 0.05$) in at least $95\%$ of bootstrap resampling iterations ($1000$ iterations) or cross-validation folds.
    3.  **Non-circularity:** The finding is not an artifact of circular prediction (e.g., predicting transcriptome-derived labels using overlapping transcriptomic features).
    4.  **Robustness to Confounders:** The association remains significant after adding tumor purity, library size, and batch effects as covariates in linear regression (e.g. MaAsLin2).

### Tier C: Exploratory Association
*   **Definition:** Findings that are statistically significant in the primary cohort under standard FDR controls but have not been validated through resampling/sensitivity analysis or replicated in external cohorts.
*   **Application:** Novel host-microbiome associations, pathway enrichments, or clinical correlations in the GSE172356/PRJNA719915 dataset.
*   **Mandatory Criteria:**
    1.  **Significance:** Benjamini-Hochberg FDR $q < 0.05$ in the main unresampled model.
    2.  **Adjustment:** Must include basic adjustment for sequencing depth and demographics (age, sex).
    3.  **Reporting:** In any manuscript, abstract, or presentation, these findings must be explicitly qualified as "exploratory," "candidate," or "requiring validation."

### Tier D: Literature-Supported Hypothesis Only
*   **Definition:** Mechanistic models, biological claims, or pathway mappings that are proposed based on external literature or predictive models, but lack statistically significant support ($q \ge 0.05$) in our data.
*   **Application:** Inferred or low-coverage functional profiling pathways (from shotgun metagenomic data without transcript/metabolite confirmation), downstream functional mechanisms of specific bacteria, or clinical survival correlations that fail significance tests.
*   **Mandatory Criteria:**
    1.  **Primary Dataset:** Non-significant test results ($q \ge 0.05$) or data not directly measurable in our assays.
    2.  **External Peer-Reviewed Evidence:** Supported by at least two independent peer-reviewed publications.
    3.  **Reporting:** Must be explicitly labeled as "speculative," "hypothesis-generating," or "literature-supported only" in all project outputs.

---

## 3. Discrepancy Reporting Policy
If an association is significant in the primary cohort (Tier C) but fails replication in an external cohort, it must not be reported as a validated finding. The discrepancy must be recorded in the final report, with a statistical comparison of the cohort characteristics (e.g., tumor purity differences, stage distribution, or sequencing platforms) to explain the lack of replication.
