# Phase 7A Tumor Microbiome & PDAC Transcriptional State Association Method Lock

This document locks the statistical framework, regression specifications, covariate hierarchies, sensitivity tests, sample-influence diagnostics, and evidence-grading rules for testing associations between the PDAC tumor microbiome and host transcriptional states. No association tests or downstream prioritizations were run prior to locking these methods.

## 1. Primary Scientific Question and Scope

The primary objective is to determine whether genus-level tumor microbiome composition is associated with the continuous basal–classical transcriptional state in PDAC. 

All analyses are restricted to the 62 mapped tumor samples from GSE172356 / PRJNA719915. No final microbiome–host association tests, differential abundance profiling, pathway enrichment, survival analysis, causal mediation, or target prioritization will be executed during Phase 7A.

---

## 2. Outcome Hierarchy (Task 2)

To prevent post hoc outcome selection or multiple-testing inflation, we lock the following outcome hierarchy:

*   **PRIMARY HOST OUTCOME:**
    *   `Moffitt50 basal–classical contrast` (calculated as the difference between the mean expression of the 25 basal-like and 25 classical genes, standardized to mean = 0, SD = 1).
*   **SECONDARY HOST OUTCOME:**
    *   `basal/classical coactivation score` (evaluating co-expression magnitude of both programs).
*   **SENSITIVITY HOST OUTCOMES:**
    *   `Moffitt49 no-LEMD1 contrast` (assesses sensitivity to the exclusion of the basal-like gene *LEMD1*).
    *   `singscore basal–classical contrast` (non-parametric rank-based single-sample contrast).
    *   `PurIST basal probability` (posterior probability of being basal-like).
    *   `Phase 4B assignment entropy` (representing transcriptomic classification uncertainty).
*   **DESCRIPTIVE OUTCOME ONLY:**
    *   `public Basal / Hybrid / Classical labels` (discrete subtyping groups from the original study). These discrete labels will not be treated as independently validated biological subtypes, and will not be used to optimize regression thresholds or selections.

---

## 3. Global Community-Level Tests (Task 1)

Before analyzing individual taxa, we will assess community-level associations using distance-based multivariate models.

*   **Distance Metric:** Primary Aitchison distance matrix (Euclidean distance on primary CLR matrix).
*   **Statistical Model:** Permutational Multivariate Analysis of Variance (PERMANOVA) using the `adonis2` function.
*   **Permutations:** 9,999 permutations (reproducible with random seed 2026).
*   **Primary Predictor:** Continuous Moffitt50 basal–classical contrast (standardized).
*   **Categorical Analyses:** For discrete public subtypes, a parallel PERMANOVA will be run. To ensure location effects are not confounded by variance differences, we will perform homogeneity of multivariate dispersions testing (**PERMDISP** / `betadisper`) and report location and dispersion results side-by-side.
*   **Sum of Squares Policy:** Sequential and marginal sums of squares will be clearly distinguished. For multi-variable models (if executed), marginal (Type III) sums of squares will be the standard.
*   **Effect Size:** Reported as $R^2$ (fraction of distance variance explained).
*   **Guardrails:** Preprocessing parameters (prev, pseudocounts, transformations) must not be selected based on the smallest PERMANOVA P-value. The public subtype PERMANOVA is descriptive only.

---

## 4. Genus-Level Primary Model (Task 2)

For each of the 122 primary retained genera, we prespecify:

*   **Primary Regression Model:**
    $$CLR\_genus \sim standardized\_Moffitt50\_contrast$$
*   **Estimation Method:** Ordinary Least Squares (OLS) regression.
*   **Standard Errors:** HC3 heteroscedasticity-robust standard errors to prevent bias from unequal variance in microbial abundances.
*   **Predictor Standardization:** The host Moffitt50 contrast score will be standardized (mean = 0, SD = 1).
*   **Coefficient Interpretation:** The OLS coefficient ($\beta$) will be interpreted as the relative change in Centered Log-Ratio (CLR) abundance per one-standard-deviation increase in the basal–classical contrast. The compositional nature of CLR coefficients must be stated explicitly in all reporting.
*   **Multiple-Testing Family:**
    *   A single family of 122 tests (one per retained genus) for association with the primary Moffitt50 contrast.
    *   Correction Method: Benjamini–Hochberg False Discovery Rate (BH FDR).
    *   Primary Significance Threshold: FDR-corrected $q < 0.05$.
    *   Reporting Requirements: For every test, we will report raw $P$-values, FDR-corrected $q$-values, regression coefficients ($\beta$), 95% confidence intervals, and standardized effect sizes (Cohen's $f^2$ or partial $R^2$).

---

## 5. Supporting Association Methods (Task 3)

To cross-validate findings from the parametric OLS model, we lock three supporting analyses:

1.  **Spearman Rank Correlation:** Assess monotonic, non-parametric relationships between CLR abundance and the standardized host score.
2.  **Permutation-Based Association Test:** Shuffling the host score 9,999 times to compute empirical, non-parametric $P$-values for each genus.
3.  **Bootstrap Confidence Intervals:** Resampling samples with replacement 2,000 times to calculate robust 95% confidence intervals for OLS coefficients and Spearman correlations.
4.  **MaAsLin2 Multivariable Modeling:** Run on CLR values with the following locked configuration:
    *   `normalization = "NONE"`
    *   `transform = "NONE"`
    *   *Constraint:* MaAsLin2 must be treated as supporting/sensitivity evidence rather than an independent replication of the same data matrix.

---

## 6. Model Hierarchy and Covariates (Task 4)

Given the sample size of $n=62$, we must avoid overfitting. We lock the following regression covariate hierarchy:

*   **Model 0 (Primary):**
    $$CLR\_genus \sim standardized\_host\_score$$
*   **Model 1 (Technical Sensitivity):**
    $$CLR\_genus \sim standardized\_host\_score + log10\_matrix\_total\_abundance\_proxy$$
    *Note:* The total sum of Bracken-normalized values is a technical library size proxy, not a direct absolute measurement of biological microbial load.
*   **Model 2 (Clinical Sensitivity):**
    $$CLR\_genus \sim standardized\_host\_score + age + sex + stage$$
    *   *Execution Constraint:* Model 2 will **ONLY** be run if:
        1.  At least 50 patients have complete covariate data.
        2.  Each categorical level (e.g., male/female, stages) contains at least 5 patients.
        3.  Collinearity diagnostics (VIF < 5) are acceptable.
    *   *Categorical Guardrail:* Clinical categories must not be collapsed post hoc after inspecting association results.
*   **Model 3P (Host TME Sensitivity: Tumor Purity):**
    $$CLR\_genus \sim standardized\_host\_score + inferred\_tumor\_purity$$
    *Execution Constraint:* Permitted after Phase 7A.5 validation because all 62 patients have complete ESTIMATE-derived purity values, maximum VIF = 1.40, and condition number = 1.82.
*   **Model 3I (Host TME Sensitivity: Immune Score):**
    $$CLR\_genus \sim standardized\_host\_score + immune\_score$$
    *Execution Constraint:* Permitted after Phase 7A.5 validation because all 62 patients have complete immune scores, maximum VIF = 1.24, and condition number = 1.61.
*   **Model 3S (Host TME Sensitivity: Stromal Score):**
    $$CLR\_genus \sim standardized\_host\_score + stromal\_score$$
    *Execution Constraint:* Permitted after Phase 7A.5 validation because all 62 patients have complete stromal scores, maximum VIF = 1.46, and condition number = 1.88.
*   **Combined TME Model Guardrail:**
    *   Tumor purity, immune score, stromal score, and ESTIMATE score must **not** be placed together in one regression model. Phase 7A.5 found severe collinearity for the combined TME screen (maximum VIF = Inf; condition number = 1.04e15; ESTIMATE score and inferred purity Spearman rho = -1.00).
    *   ESTIMATE-derived purity and cell-composition scores are inferred from the same host transcriptomic matrix as the Moffitt scores. They are sensitivity covariates, not independent experimental measurements. Adjustment may remove biological variation genuinely associated with the PDAC transcriptional state. Therefore, Model 0 remains primary and Models 3P/3I/3S assess robustness only.

---

## 7. Contamination Sensitivity (Task 5)

Primary genus-level testing retains all 122 genera. To ensure findings are not driven by background contamination, we lock five sensitivity checks:

1.  **Excluding HIGH_RISK_POTENTIAL_CONTAMINANT Taxa:** Rerun OLS and Spearman models excluding the 6 high-risk genera (`Elizabethkingia`, `Delftia`, `Brevundimonas`, `Comamonas`, `Caulobacter`, `Ralstonia`).
2.  **Excluding HIGH_RISK & MODERATE_RISK Taxa:** Rerun excluding both high-risk and moderate-risk environmental genera (`Paraburkholderia`, `Mesorhizobium`, `Novosphingobium`, `Dechloromonas`, `Sphingopyxis`, `Herbaspirillum`).
3.  **Leave-One-Genus-Out (LOGO) Checks:** For any top significant genus, recalculate sample CLR values and distance matrices without that genus, verifying that the association is not dependent on a single taxon dominating the geometric mean.
4.  **Matrix Total-Abundance Proxy Correlation:** Directly correlate each genus's CLR abundance with the total abundance proxy. If a significant taxon correlates strongly ($\rho > 0.5$, $P < 0.01$), flag it as potentially contamination-sensitive.
5.  **Flag-Specific Reporting:** Biologically plausible genera (e.g., `Pseudomonas`, `Acinetobacter`) that are also common contaminants will be flagged and discussed, but no genus will be labeled as "confirmed contamination" without experimental controls.

---

## 8. Preprocessing Sensitivity (Task 6)

For every candidate genus-level association, we will rerun OLS and Spearman models across the 9 pre-computed Phase 6C preprocessing sensitivity matrices:

*   10% prevalence filter matrix
*   30% prevalence filter matrix
*   Abundance > 10 detection matrix
*   Pseudocount 0.1 matrix
*   Pseudocount 1.0 matrix
*   Robust CLR (rCLR) matrix
*   Matrix excluding three technical extreme samples (`Basal-like1`, `Hybrid18`, `Hybrid23`)
*   Matrix excluding High-Risk contaminants
*   Matrix excluding High- plus Moderate-Risk contaminants

For each candidate association, we will report:
*   Sign consistency across all sensitivity runs.
*   The range of OLS coefficients ($\beta$).
*   The range of FDR-adjusted $q$-values.
*   The number and fraction of sensitivity analyses yielding the same direction.
*   Whether statistical significance depends on a single preprocessing choice.

---

## 9. Presence/Absence Sensitivity (Task 7)

For sufficiently prevalent but sparse genera, we define a binary model to test whether presence/absence (rather than CLR abundance) associates with host state.

*   **Locked Orientation:**
    $$genus\_presence\_absence \sim standardized\_host\_score$$
    using **Logistic Regression** (presence defined as Bracken abundance $> 0.0$).
*   **Sample Size Constraint:** Logistic regression will only be run for genera that are present in at least 10 samples **and** absent in at least 10 samples (i.e. prevalence between 16.1% and 83.9% of the 62-sample cohort). This avoids complete or quasi-complete separation.
*   **Rare Taxa Guardrail:** Do not apply presence/absence models to rare genera (prevalence < 16.1%) to prevent statistical separation artifacts.

---

## 10. Evidence Classification (Task 8)

Associations will be classified internally into five locked evidence categories:

*   **ROBUST_ASSOCIATION:**
    1. Primary OLS Benjamini–Hochberg $q < 0.05$.
    2. Bootstrap 95% confidence interval for OLS coefficient excludes zero.
    3. Same effect direction in Spearman rank correlation.
    4. Same direction in at least 75% (7 of 9) of locked preprocessing sensitivities.
    5. Association is not eliminated solely by removing one technical extreme sample.
    6. Contamination risk status explicitly evaluated and noted.
*   **SUGGESTIVE_ASSOCIATION:**
    *   Primary OLS BH $q$ is between 0.05 and 0.10, OR
    *   Primary OLS BH $q < 0.05$ but limited support (fails at least two sensitivity rules).
*   **METHOD_SENSITIVE:**
    *   Association shows a sign reversal or major rank change across reasonable preprocessing alternatives.
*   **CONTAMINATION_SENSITIVE:**
    *   Signal disappears after high-risk/moderate-risk taxon exclusion or is strongly correlated with the total-abundance proxy.
*   **NO_SUPPORTED_ASSOCIATION:**
    *   Insufficient statistical and robustness evidence (fails primary FDR and suggestive thresholds).
*   **TO_VERIFY:**
    *   Incomplete metadata or implementation uncertainty.

*Caution:* These are internal evidence categories for this project and do not represent external biological validation.

---

## 11. Secondary Outcomes (Task 9)

For coactivation, PurIST probability, Moffitt49, singscore, and entropy:
*   We will run separate hypothesis families of 122 tests.
*   Apply Benjamini–Hochberg FDR correction independently within each outcome's family.
*   Clearly report these as secondary or sensitivity analyses.
*   Do not pool or average P-values/q-values across outcomes as independent replication.
*   Require directionally coherent interpretation before drawing conclusions.

---

## 12. Descriptive Discrete Subtype Analyses (Task 10)

For the public Basal / Hybrid / Classical labels, we prespecify:
*   Descriptive community-level Aitchison PERMANOVA and PERMDISP.
*   Genus-level Kruskal–Wallis tests (122 tests, corrected using BH FDR).
*   Effect size reported as $\eta^2$ or $\epsilon^2$.
*   Post hoc pairwise Dunn's tests (with BH correction) executed **only** if the overall Kruskal–Wallis test is significant ($q < 0.05$).
*   These results must remain secondary and descriptive.

---

## 13. Sample-Level Influence (Task 11)

To ensure that findings are not driven by individual outlier patients, we will calculate:
*   Cook's distance for each OLS regression (flagging samples with $D_i > 4/n$).
*   DFBETAs (flagging samples with $|DFBETA_i| > 2/\sqrt{n}$).
*   Leave-one-sample-out sensitivity analysis.
*   Specific exclusion sensitivity for the three technical extreme samples: `Basal-like1`, `Hybrid18`, and `Hybrid23`.
*   *Reporting Rule:* We will explicitly report if any candidate association is driven by one or two patients. We will **not** automatically remove influential samples from the primary analysis.

---

## 14. Negative and Null Results (Task 12)

We enforce the publication-standard reporting of negative and null results:
*   If no global community association is found, it must be reported.
*   If no individual genus meets the primary FDR threshold, it must be stated.
*   If effect directions are inconsistent across compositional transformations, this must be documented.
*   Associations lost after contamination checks will be presented as null.
*   *Strict Rule:* Nominal P-values will not be promoted as "discoveries" or "trends" when the primary FDR threshold is not met.

---

## 15. Locked Computational Settings

*   **Permutations:** 9,999 for PERMANOVA and permutation tests.
*   **Bootstrap Iterations:** 2,000 for confidence intervals.
*   **Random Seed:** 2026.
*   **Testing Mode:** Two-sided tests for all hypothesis tests.
*   **Multiple Testing:** Benjamini–Hochberg FDR correction.

---

## 16. Overall Readiness Decision

Based on the locked specifications and available data:

*   **Decision:** `READY_WITH_HOST_TME_SENSITIVITY_MODELS`
*   **Justification:** 
    *   The primary filtered microbiome matrices (122 genera x 62 samples), continuous host scores, and public labels are finalized and validated.
    *   Clinical metadata (age, sex, stage) remain 100% missing in the clinical series; therefore, **Model 2 (clinical sensitivity) is not permitted** (as shown in the feasibility table `phase7a_model_matrix_feasibility.tsv`).
    *   Phase 7A.5 calculated and validated ESTIMATE-derived stromal score, immune score, ESTIMATE score, and inferred tumor purity for all 62 host-expression samples without using subtype or microbiome information.
    *   Downstream Phase 7 association tests may proceed with Model 0 (primary), Model 1 (technical total-abundance sensitivity), and the feasible host TME sensitivity models Model 3P, Model 3I, and Model 3S. Combined adjustment for purity, immune score, stromal score, and ESTIMATE score is blocked by collinearity.
