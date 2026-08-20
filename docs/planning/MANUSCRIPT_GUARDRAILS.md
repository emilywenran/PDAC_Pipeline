# Manuscript Guardrails

This document establishes strict writing guidelines, language constraints, and reporting standards to ensure scientific integrity and prevent overinterpretation of results in any manuscript drafts, reports, or presentations arising from this project.

---

## 1. Language Constraints: Differentiating Association from Causation

Because the data in GSE172356 and PRJNA719915 are cross-sectional (sampled at a single clinical timepoint post-resection), it is impossible to determine the temporal sequence of events or isolate causal mechanisms. Therefore, causal language is strictly prohibited.

*   **Forbidden Causal Verbs:** Do *not* use words like:
    *   *drives, causes, triggers, induces, promotes, accelerates, dictates, determines, mediates, shapes.*
    *   *Example of a prohibited statement:* "The high abundance of *Acinetobacter* drives the basal-like subtype phenotype."
*   **Approved Associative Verbs:** Use objective, descriptive terms such as:
    *   *is associated with, correlates with, co-occurs with, is enriched in, exhibits co-abundance with, is linked to.*
    *   *Example of an approved statement:* "The abundance of *Acinetobacter* was positively correlated with the host continuous basal-like score."

---

## 2. Biomarker Classification Rules

Feature selection algorithms (e.g., Lasso, ElasticNet, random forest feature importances, or differential expression analysis) run on our primary cohort identify candidate variables, not verified clinical biomarkers.

*   **Rule:** Do *not* describe internally selected features as "validated biomarkers" or "diagnostic biomarkers."
*   **Terminology:** Refer to them as "candidate features," "exploratory predictors," "candidate transcriptional markers," or "subtype-associated features."
*   **Upgrade Criteria:** A feature may only be described as a "validated biomarker" if it achieves **Tier A** evidence (statistically significant validation in an independent external cohort with controlled covariates).

---

## 3. Validation Reporting Standards

Cross-validation (including k-fold, leave-one-out, or nested cross-validation) is a method for estimating model generalization and tuning hyperparameters. It is part of model training and does *not* constitute independent validation.

*   **Rule:** Do *not* report cross-validation results as "independent validation" or "external validation."
*   **Reporting:** Clearly label cross-validation metrics as "internal cross-validation performance" or "internal generalization estimate."
*   **Validation Standard:** Reserve the term "independent validation" or "external validation" exclusively for models trained on the primary dataset and tested on a completely separate, independent cohort that was not accessed during feature selection, normalization tuning, or model fitting.
*   **External Validation Guidelines:**
    1.  TCGA-PAAD should primarily support host transcriptomic validation.
    2.  Any TCGA-derived microbiome signal must be treated as exploratory because of extensive contamination and batch concerns.
    3.  Do not assume that an identical external paired transcriptome–microbiome cohort exists for direct validation of host-microbiome interactions.

---

## 4. Citation and Citation-Isolation Rules

To ensure scientific independence and prevent the dissemination of unverified claims:

*   **Rule:** Do *not* cite, discuss, or reference any user-provided unpublished manuscripts, drafts, or private communications.
*   **Source Integrity:** Base all manuscript statements and biological claims exclusively on:
    1.  The primary statistical and bioinformatic analyses performed in this project.
    2.  Independently verified, peer-reviewed public literature (indexed in NCBI PubMed, Europe PMC) or officially deposited preprints (arXiv, bioRxiv).
*   All external claims must have verifiable DOIs.

---

## 5. Reporting Negative and Non-Significant Results

Failing to reject the null hypothesis is a core scientific finding that prevents publication bias.

*   **Rule:** All negative, non-significant, or neutral results must be reported with the same prominence and detail as positive findings.
*   **Implementation:**
    *   If the primary candidate taxa (*Acinetobacter*, *Pseudomonas*, *Sphingopyxis*) show no significant association with subtypes after adjusting for tumor purity, this must be explicitly stated in the abstract, results, and discussion.
    *   Include full tables of non-significant tests (with exact p-values, FDR values, and effect sizes) in the supplementary materials.

---

## 6. Separating Prespecified and Exploratory Analyses

To prevent "p-hacking" and data dredging, all analyses must be clearly categorized:

*   **Prespecified Analyses:** Those defined in the primary hypotheses and analysis plans (e.g., reproducing Guo et al. subtypes, testing continuous scores, running MaAsLin2 on the targeted taxa, validating host signatures in TCGA-PAAD). These must be reported first in the manuscript.
*   **Exploratory Analyses:** Post-hoc analyses prompted by patterns observed during the study (e.g., testing unexpected pathway interactions or searching for novel microbial biomarkers). These must be clearly segregated into a section titled "Exploratory Post-Hoc Analyses" and qualified as hypothesis-generating only.
