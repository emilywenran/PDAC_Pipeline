# Phase 8A Host Mechanism Method Lock

This document locks the prospective statistical and bioinformatics framework for identifying host-mechanism associations linked to tumor-microbiome features. No association tests, enrichment analyses, pathway scoring, or network constructions are executed in this phase.

---

## 1. Executive Summary & Readiness Decision

*   **OVERALL READINESS DECISION:** **`READY_WITH_TRANSFORMATION_LIMITATIONS`**
*   **RATIONALE:** 8 of the 9 robust tumor-associated genera (`Azoarcus`, `Candida`, `Ensifer`, `Cutibacterium`, `Chryseobacterium`, `Burkholderia`, `Rhizobium`, `Herbaspirillum`) exhibit direction reversals under the locked robust Centered Log-Ratio (rCLR) sensitivity representation. Downstream mechanisms are feasible to test, but biological interpretations must report both primary CLR and rCLR results, and any link where the interpretation reverses under rCLR is locked as a transformation-sensitive mechanism.
*   **PRIMARY COHORT:** $n=62$ PDAC patients with matched host transcriptomics and microbiome metagenomics.

---

## 2. Taxon Candidate Hierarchy

Downstream analyses will be conducted strictly using the taxonomic classifications locked in Phase 7C:

*   **PRIMARY TAXA (Robust Associations):**
    *   *Azoarcus*
    *   *Candida*
    *   *Ensifer*
    *   *Cutibacterium*
    *   *Chryseobacterium*
    *   *Ochrobactrum*
    *   *Burkholderia*
    *   *Rhizobium*
    *   *Herbaspirillum*
*   **SECONDARY TAXA (Suggestive Associations):**
    *   *Staphylococcus* (downgraded in Phase 7C due to sample leverage concern)
    *   *Citrobacter* (nominal q-value between 0.05 and 0.10)

*   **GUARDRAIL:** No method-sensitive or contamination-sensitive genera from Phase 7C will be promoted to the primary mechanism analysis.

---

## 3. Five-Layer Host Feature Hierarchy

### Layer 1: Predefined Pathway Activity
*   **Primary Collections:** MSigDB Hallmark (50 gene sets), PROGENy pathway activities (14 pathways).
*   **Secondary Collections:** Reactome, KEGG (only where gene-set interpretation is biologically appropriate).
*   **Constraint:** Screening an unlimited number of overlapping pathway databases and reporting only favorable results is strictly prohibited.

### Layer 2: Transcription-Factor Activity
*   **Primary Framework:** DoRothEA/VIPER regulon activities.
*   **Locked Parameters:**
    *   Regulon confidence levels restricted to high-confidence categories (A, B, and C).
    *   Minimum target-gene coverage per regulon set to $\ge 15$ genes.
    *   Input expression scale: log2 analysis-ready expression.
    *   Activity-score direction: VIPER normalized enrichment score (NES) representing positive/negative transcription-factor activity.
    *   Multiple-testing correction: Benjamini-Hochberg FDR applied separately within each transcription-factor family.

### Layer 3: Tumor Microenvironment (TME) Programs
*   **Primary Scores:** Inferred tumor purity, immune score, stromal score, and ESTIMATE score.
*   **Source:** Official MD Anderson ESTIMATE package scores validated in Phase 7A.5.
*   **Secondary Frameworks:** Alternate deconvolution methods (e.g., CIBERSORTx, EPIC) are classified as secondary exploratory analyses only and must be verified before execution.

### Layer 4: Data-Driven Co-Expression Modules
*   **Primary Framework:** Weighted Gene Co-expression Network Analysis (WGCNA).
*   **Locked Parameters (Prior to Execution):**
    *   Gene filtering: Keep top 25% most variable genes based on median absolute deviation (MAD).
    *   Transformation: Log2 analysis-ready expression values.
    *   Soft-threshold selection rule: Scale-free topology fit index $R^2 \ge 0.85$, selecting the lowest power meeting this criterion.
    *   Minimum module size: 30 genes.
    *   Module merging threshold: Merge modules with eigengene correlation $\ge 0.80$ (dendrogram cut height 0.20).
    *   Preservation and robustness checks: Module preservation Z-summary score ($Z_{summary} > 10$ for high preservation, $2 < Z_{summary} < 10$ for moderate) calculated using random split-cohort bootstrap resampling.
*   **Constraint:** Public subtype labels and microbiome association results must not be used to tune network parameters or select soft-thresholds.

### Layer 5: Taxon-Associated Host Genes
*   **Model Formula:**
    $$host\_gene\_expression \sim standardized\_CLR\_genus$$
*   **Execution Policy:** Genome-wide host-gene association is exploratory and high-dimensional ($42,654$ genes).
*   **Linear Modeling:** Ordinary least squares (OLS) with HC3 robust standard errors, or empirical Bayes moderation (limma-style) if variance shrinkage is required.
*   **Covariates:** Separate sensitivity models run adjusting for inferred tumor purity, immune score, and stromal score.
*   **Correction:** Benjamini-Hochberg FDR correction applied separately for each taxon. All effect sizes and 95% confidence intervals must be reported.
*   **Constraint:** Combining all taxa and genes into a single undifferentiated testing family is prohibited.

---

## 4. Circularity Safeguards

1.  **Moffitt50 Independence Guardrail:** Because the Moffitt50 signature genes were used to define the continous transcriptional axis in Phase 5B and evaluate microbiome associations in Phase 7B, correlations with Moffitt signature genes cannot be reported as independent biological validation.
2.  **Analysis Reporting Policy:** All downstream host-mechanism analyses (pathway, TF, and gene-level) must be reported in two formats:
    *   Including all host genes/pathways.
    *   Excluding the Moffitt50 signature genes (where technically applicable) to verify that findings are not driven solely by the direct reuse of the subtyping genes.
3.  **Distinguish Evidence Sources:** The final report must clearly distinguish:
    *   Direct reuse of subtype-signature genes.
    *   Independent pathway or TF evidence.
    *   Confounding or composition effects driven by immune/stromal infiltration.
4.  **No Post-hoc Gene-Set Construction:** Constructing custom pathway gene sets from Phase 7B-associated genes and testing them on the same PDAC cohort is exploratory and must be explicitly documented as such.

---

## 5. Primary Statistical Model & Sensitivities

### Primary Model Formula
For every predefined pathway, TF activity, TME score, or WGCNA module eigengene:
$$host\_activity \sim standardized\_CLR\_genus$$

*   **Test Type:** Two-sided hypothesis tests.
*   **Standard Errors:** HC3 heteroscedasticity-robust standard errors.
*   **Effect Size:** Standardized regression coefficient ($\beta$) and its 95% confidence interval.
*   **Multiple-Testing Correction:** Benjamini-Hochberg (BH) FDR applied within each predefined host-feature collection (e.g., MSigDB Hallmark, DoRothEA) separately for each taxon.

### Sensitivity Models
Each primary model must be evaluated against the following locked sensitivity adjustments:
1.  **Tumor Purity Adjustment:** `host_activity ~ standardized_CLR_genus + inferred_tumor_purity`
2.  **Immune Score Adjustment:** `host_activity ~ standardized_CLR_genus + immune_score`
3.  **Stromal Score Adjustment:** `host_activity ~ standardized_CLR_genus + stromal_score`
4.  **rCLR Representation:** Substitute primary CLR with the robust CLR (rCLR) genus representation.
5.  **Technical Outlier Exclusion:** Re-estimate model parameters after removing the three extreme samples (`Basal-like1`, `Hybrid18`, `Hybrid23`).
6.  **Contaminant Exclusion:** Re-estimate models using taxonomic abundances computed from contaminant-exclusion representations (locked in Phase 6C).
7.  **Leave-One-Sample-Out (LOO):** Iteratively drop one sample at a time to verify that associations are not driven by single patients.

*   **Constraint:** Inferred tumor purity, immune score, stromal score, and ESTIMATE score must never be included in the same regression model due to extreme collinearity.

---

## 6. Multi-Taxon Shared-Mechanism Analysis

To distinguish shared biological mechanisms from compositional correlation:
1.  **Cross-Taxon Sign-Consistency:** Summarize shared pathways/TFs by checking whether coefficients across multiple genera have consistent directional signs.
2.  **Multivariable Taxon Models:** Run joint regressions:
    $$host\_activity \sim CLR\_genus\_A + CLR\_genus\_B$$
    only if collinearity (VIF $< 5$, condition number $< 30$) and sample size criteria are satisfied.
3.  **Exploratory Sparse Multivariate Methods:** Apply Lasso/ElasticNet regression only as exploratory analyses to select independent microbial predictors of host programs.
4.  **Biological Interpretation Constraint:** Compositionally correlated taxa must not be interpreted as independent biological exposures.

---

## 7. Host-Mechanism Evidence Categories

Every candidate association must be classified into one of the following locked evidence categories:

*   **`ROBUST_HOST_MECHANISM`:**
    *   Pathway/TF activity passes the locked FDR threshold ($q < 0.05$).
    *   95% confidence interval excludes zero.
    *   Direction is stable across all LOO sensitivity runs.
    *   Effect is not driven by one or two extreme samples.
    *   Association remains significant ($P < 0.05$) after adjusting for tumor purity, immune score, or stromal score.
    *   rCLR taxon transformation does not reverse the biological sign or eliminate the association.
*   **`TRANSFORMATION_SENSITIVE_MECHANISM`:**
    *   Primary CLR association meets robust criteria, but rCLR transformation reverses the coefficient sign or removes statistical support.
*   **`COMPOSITION_SENSITIVE_MECHANISM`:**
    *   Primary CLR association is supported, but coefficient becomes non-significant after adjustment for tumor purity, immune, or stromal scores.
*   **`SAMPLE_SENSITIVE_MECHANISM`:**
    *   Association significance or coefficient direction depends on the inclusion of one or two specific patients.
*   **`EXPLORATORY_HOST_MECHANISM`:**
    *   Nominal significance ($P < 0.05$, $q \ge 0.05$) or secondary pathway support without satisfying all robustness criteria.
*   **`NO_SUPPORTED_MECHANISM`:**
    *   Insufficient evidence after multiple-testing correction and sensitivity analyses.
*   **`TO_VERIFY`:**
    *   Incomplete package execution, missing inputs, or validation script failures.

---

## 8. Causal Interpretation and Literature Guardrails

### Mechanistic Interpretation Limits
*   **Permitted Verbs:** "associated with", "consistent with", "suggests a potential link".
*   **Prohibited Verbs:** "drives", "causes", "mediates", "induces".
*   **Constraint:** Causal mediation analysis is not permitted in this project phase due to the cross-sectional, observational nature of the dataset.

### Literature Use Policy
*   **Sequencing:** Literature search via `citation-management` must be performed **only after** computational candidate mechanisms are defined and categorized.
*   **Literature Categorization:** For each citation, the discussion must explicitly distinguish:
    *   Evidence generated in this project.
    *   Published PDAC evidence.
    *   Published evidence from other cancer types.
    *   Speculative microbial functional mechanisms.
*   **Rescue Prohibition:** Literature citations must never be used to rescue or support computationally non-significant or unstable findings.
*   **Draft Manuscripts:** Citing or discussing any user-provided draft manuscript is strictly prohibited.

---

## 9. Computational Defaults

Unless an existing locked protocol requires otherwise, the following defaults are locked:
*   **Random Seed:** 2026.
*   **Hypothesis Tests:** Two-sided.
*   **FDR Method:** Benjamini-Hochberg (BH).
*   **Bootstrap Iterations:** 2,000 for selected candidate mechanisms to compute robust confidence intervals.
*   **Threshold Optimization:** No data-driven threshold optimization (e.g., grid-searching p-value thresholds or pathway cutoffs) is permitted.
