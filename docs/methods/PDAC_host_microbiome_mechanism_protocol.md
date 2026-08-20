# Standard Operating Procedure: PDAC Host-Microbiome Mechanism Protocol

This document defines the step-by-step execution protocol for analyzing host pathways, transcription-factor activities, tumor microenvironment (TME) programs, data-driven co-expression modules, and genome-wide host-gene associations linked to verified tumor-associated microbial taxa.

---

## 1. Prerequisites and Environmental Setup

### 1.1 Python and R Environments
*   **Python Version:** $\ge 3.10$ (Python 3.13.9 base conda environment recommended).
*   **R Version:** $\ge 4.0$ (R 4.5.3 recommended).
*   **Required Python Packages (Pre-execution audit):**
    *   `statsmodels` ($\ge 0.14.5$, available) — OLS regression and HC3 standard error estimation.
    *   `pandas`, `numpy`, `scipy` — data handling and basic statistics.
    *   `decoupler` (unavailable) — python interface for PROGENy, DoRothEA, and MSigDB.
    *   `gseapy` (unavailable) — python interface for pathway enrichment.
*   **Required R Packages (Pre-execution audit):**
    *   `limma` (available) — empirical Bayes moderated linear models.
    *   `estimate` (available) — ESTIMATE score calculation.
    *   `decoupleR` (unavailable) — R pathway and regulon activity calculation.
    *   `WGCNA` (unavailable) — weighted gene co-expression network analysis.

### 1.2 Input File Requirements
All input files must reside under their designated folders and have their MD5/SHA256 checksums verified before run:
1.  **Host Expression Matrix:** `03_processed/expression/GSE172356_expression_log2_analysis_ready.tsv.gz` (dimensions: $42,654$ genes $\times 62$ samples).
2.  **Microbiome Abundance Transform:** Centered Log-Ratio (CLR) and robust Centered Log-Ratio (rCLR) matrices under `03_processed/microbiome/`.
3.  **TME Covariates:** `01_metadata/host_tme_covariates.tsv` containing ESTIMATE stromal, immune, and inferred tumor-purity scores for all 62 patients.
4.  **Sample Manifest:** `01_metadata/sample_manifest.tsv` (62 patient rows).

---

## 2. Step-by-Step Execution Protocol

### Step 2.1: Predefined Pathway and TF Activity Calculation
*Once the required `decoupleR` or `gseapy` packages are installed/made available:*

1.  **MSigDB Hallmark (Layer 1):**
    *   Retrieve the 50 Hallmark gene sets from MSigDB.
    *   Use the Single-Sample GSEA (ssGSEA) algorithm or decoupleR's `run_mlm` (multivariate linear model) to calculate sample-wise pathway scores.
    *   Save scores in `03_processed/expression/host_pathway_hallmark_scores.tsv`.
2.  **PROGENy Pathway Activities (Layer 1):**
    *   Compute progeny activities using decoupleR/progeny for the 14 standard pathways.
    *   Save scores in `03_processed/expression/host_pathway_progeny_scores.tsv`.
3.  **DoRothEA Transcription-Factor Activities (Layer 2):**
    *   Load DoRothEA regulons. Filter to confidence levels A, B, and C only.
    *   Calculate activities using the VIPER algorithm (`run_viper` in decoupleR).
    *   Save scores in `03_processed/expression/host_tf_dorothea_activities.tsv`.

### Step 2.2: Tumor Microenvironment Score Integration (Layer 3)
1.  Verify that stromal, immune, ESTIMATE, and inferred tumor-purity scores are correctly loaded from `01_metadata/host_tme_covariates.tsv`.
2.  No further deconvolution is performed in the primary pipeline. Secondary methods (such as CIBERSORTx) must be validated separately and logged.

### Step 2.3: Data-Driven Co-Expression Network Analysis (Layer 4)
*Once `WGCNA` (R) or `PyWGCNA` (Python) is installed/made available:*

1.  **Filter Genes:** Keep the top 25% most variable genes based on Median Absolute Deviation (MAD).
2.  **Verify No Missingness:** Confirm no missing values in the filtered matrix.
3.  **Select Power:** Calculate scale-free topology fit index ($R^2$) for powers $1$ to $30$. Select the lowest power $\beta$ where $R^2 \ge 0.85$. Do not optimize this power using public subtype labels or microbiome associations.
4.  **Construct Network:** Generate a signed hybrid network using the selected soft-threshold power.
5.  **Identify Modules:** Use hierarchical clustering and the Dynamic Tree Cut algorithm with a minimum module size of 30 genes.
6.  **Merge Modules:** Calculate module eigengenes (MEs) and merge modules with eigengene correlation $\ge 0.80$ (dendrogram height threshold 0.20).
7.  **Preservation:** Run module preservation checks using 100 permutation rounds, and output Z-summary statistics.
8.  **Save MEs:** Write module eigengenes for all samples to `03_processed/expression/host_wgcna_module_eigengenes.tsv`.

### Step 2.4: Genome-wide Host-Gene Association (Layer 5)
1.  Loop through each primary candidate taxon.
2.  For each of the $42,654$ host genes, fit the OLS regression:
    $$host\_gene\_expression \sim standardized\_CLR\_genus$$
3.  Calculate p-values and apply Benjamini-Hochberg FDR correction separately within the gene family for that specific taxon.
4.  Save coefficients, robust SEs, raw p-values, q-values, and 95% CIs in `05_results/tables/phase8b_gene_associations_[taxon].tsv`.

---

## 3. Statistical Testing and Model Specifications

For each host feature (Hallmark score, PROGENy score, TF activity, TME score, or WGCNA module eigengene):

### 3.1 Primary Model
Fit the ordinary least squares (OLS) regression:
$$host\_activity \sim standardized\_CLR\_genus$$

*   **Standard Errors:** Heteroscedasticity-robust standard errors (HC3).
*   **Hypothesis Testing:** Two-sided tests.
*   **Correction:** Benjamini-Hochberg (BH) FDR correction, applied within each feature collection (family) separately for each taxon.
*   **FDR Alpha:** $q < 0.05$ defines statistical significance.

### 3.2 Sensitivity Adjustments
For each significant association ($q < 0.05$), run the following sensitivity models:
1.  **Tumor Purity Control:** Add ESTIMATE inferred purity as a covariate:
    $$host\_activity \sim standardized\_CLR\_genus + inferred\_tumor\_purity$$
2.  **Immune Control:** Add ESTIMATE immune score as a covariate:
    $$host\_activity \sim standardized\_CLR\_genus + immune\_score$$
3.  **Stromal Control:** Add ESTIMATE stromal score as a covariate:
    $$host\_activity \sim standardized\_CLR\_genus + stromal\_score$$
4.  **Composition Transformation Sensitivity (rCLR):** Re-run the model replacing primary CLR with rCLR genus representation. Flag any results as `RCLR_DIRECTION_SENSITIVE` if the sign of the coefficient reverses under rCLR.
5.  **Technical Outliers:** Remove the three technical extreme samples (`Basal-like1`, `Hybrid18`, `Hybrid23`) and re-estimate.
6.  **Contamination Control:** Re-run model using CLR calculated from contaminant-excluded abundance matrices.
7.  **Sample Influence (LOO):** Run 62 separate regressions dropping one sample at a time. Record sign changes or loss of significance ($P \ge 0.05$).

*   **Collinearity Guardrail:** Never include inferred tumor purity, immune score, stromal score, and ESTIMATE score in the same model.

---

## 4. Circularity Safeguards and Control Analyses

1.  **Moffitt50 Gene Exclusion:**
    *   To verify that pathway/TF associations are not driven by the direct reuse of the subtype-defining genes, compile the list of 50 Moffitt signature genes.
    *   For any significant pathway or TF activity, re-calculate the activity scores after removing all Moffitt50 genes from the underlying gene sets.
    *   Re-run the OLS regression and compare coefficients. If an association loses significance, classify it as subtype-dependent and document it in the report.
2.  **No Custom Pathway Enrichment:** Do not construct custom pathway definitions from Phase 7B-associated genes and then run enrichment tests on the same cohort without clearly documenting the analysis as exploratory.

---

## 5. Multi-Taxon Integration

1.  **Sign Consistency Audit:** Summarize the overlap of host pathways and TF activities across all 9 robust genera. Build a directional matrix mapping each pathway to the sign of its coefficient (+/-) for each genus.
2.  **Collinearity Screen:** If multiple genera are associated with the same pathway, test their collinearity:
    *   Compute VIF and condition numbers for a model containing both genera.
    *   If VIF $\ge 5$ or condition number $\ge 30$, do not run the multivariable model. Describe the association as compositionally correlated.
3.  **Multivariable OLS (Only if permitted by VIF):**
    $$host\_activity \sim CLR\_genus\_A + CLR\_genus\_B$$

---

## 6. Classification of Evidence

Classify every candidate host-mechanism association into one of the following locked evidence categories:

1.  **`ROBUST_HOST_MECHANISM`:** Pass primary BH FDR ($q < 0.05$), CI excludes zero, stable direction under LOO, not driven by technical outliers, remains significant after purity/immune/stromal adjustments, and rCLR representation does not reverse the biological sign.
2.  **`TRANSFORMATION_SENSITIVE_MECHANISM`:** Meets all primary criteria, but rCLR transformation reverses the coefficient sign or eliminates significance.
3.  **`COMPOSITION_SENSITIVE_MECHANISM`:** Primary association is significant, but becomes non-significant ($P \ge 0.05$) when adjusting for tumor purity, immune, or stromal scores.
4.  **`SAMPLE_SENSITIVE_MECHANISM`:** Association significance or coefficient direction depends on the inclusion of one or two specific samples.
5.  **`EXPLORATORY_HOST_MECHANISM`:** Nominal significance ($P < 0.05$, $q \ge 0.05$) or secondary pathway support without satisfying all robustness criteria.
6.  **`NO_SUPPORTED_MECHANISM`:** Insufficient evidence after multiple-testing correction and sensitivity analyses.
7.  **`TO_VERIFY`:** Incomplete package execution, missing inputs, or validation script failures.

---

## 7. Reporting Guidelines and Interpretation Limits

1.  **Causal Language Guardrails:**
    *   The final report must use associative verbs: "associated with", "consistent with", "suggests a potential link".
    *   The final report must not use causal verbs: "drives", "causes", "mediates", "induces".
    *   No causal mediation analysis is permitted in this phase.
2.  **Literature Integration:**
    *   Literature search via the `citation-management` skill must be performed **only after** computational candidate mechanisms are defined and categorized.
    *   For each citation, the discussion must explicitly distinguish:
        *   Evidence generated in this project.
        *   Published PDAC evidence.
        *   Published evidence from other cancer types.
        *   Speculative microbial functional mechanisms.
    *   No literature citations can be used to rescue or support computationally non-significant or unstable findings.
    *   Do not cite or discuss any user-provided draft manuscript.
