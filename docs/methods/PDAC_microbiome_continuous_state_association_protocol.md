# SOP: PDAC Tumor Microbiome & Continuous Transcriptional State Association Protocol

This protocol defines the step-by-step statistical procedures for evaluating associations between the tumor microbiome (at the genus level) and continuous transcriptional states in the 62-patient PDAC cohort (GSE172356 / PRJNA719915).

---

## Step 1: Input Data Preparation & Verification

Before executing any models, load and check the checksums and dimensions of the finalized inputs:

1.  **Microbiome Input Matrix:** 
    *   `03_processed/microbiome/PRJNA719915_genus_primary_CLR.tsv.gz` (122 genera x 62 samples).
    *   `03_processed/microbiome/PRJNA719915_primary_aitchison_distance.tsv.gz` (62 x 62 symmetric distance matrix).
2.  **Host Transcriptional Outcomes:** 
    *   `05_results/tables/phase5b_sample_continuous_scores.tsv` (contains continuous `basal_classical_contrast`, `coactivation_score`, `purist_basal_probability`, etc.).
3.  **Crosswalk Mapping:** 
    *   `01_metadata/microbiome_sample_crosswalk.tsv` (maps BioSample and SRA run accessions to host `patient_id` / tumor number).
4.  **Clinical Metadata:**
    *   `01_metadata/clinical_metadata.tsv` (confirm clinical fields completeness).
5.  **Host TME Sensitivity Covariates:**
    *   `01_metadata/host_tme_covariates.tsv` (ESTIMATE-derived stromal score, immune score, ESTIMATE score, and inferred tumor purity).
    *   Use only Phase 7A.5-permitted sensitivity covariates: inferred tumor purity, immune score, and stromal score in separate models.

*Verification Check:* Ensure all 62 patient rows are aligned exactly across host outcomes, clinical metadata, and microbiome columns.

---

## Step 2: Global Community-Level Association (PERMANOVA / PERMDISP)

Assess whether overall microbiome community composition is associated with the continuous basal–classical axis:

1.  **PERMANOVA Specification:**
    *   Function: `vegan::adonis2` (R) or `skbio.stats.distance.permanova` (Python).
    *   Formula: `Aitchison_Distance ~ standardized_Moffitt50_contrast`
    *   Permutations: 9,999.
    *   Random Seed: 2026.
    *   Report: $P$-value, pseudo-$F$ statistic, and $R^2$ (effect size).
2.  **Discrete Subtypes PERMANOVA/PERMDISP (Descriptive Only):**
    *   Formula: `Aitchison_Distance ~ public_subtype_label`
    *   Evaluate dispersion: `vegan::betadisper` (PERMDISP). Report whether within-group dispersion differences are statistically significant ($P < 0.05$).
    *   *Guardrail:* Report PERMANOVA location and PERMDISP dispersion results together. Do not interpret location differences as biological if dispersion is highly heterogeneous.

---

## Step 3: Genus-Level Primary Model (OLS with HC3 Robust SEs)

For each of the 122 primary retained genera, fit the primary linear model:

1.  **Standardize Predictor:** 
    *   Center and scale the Moffitt50 basal-classical contrast:
        $$x_i = \frac{Score_i - \bar{Score}}{SD(Score)}$$
2.  **Fit Model:**
    *   Formula: `CLR_genus ~ standardized_Moffitt50_contrast`
    *   Method: OLS.
    *   *Requirement:* Calculate standard errors, confidence intervals, and $P$-values using **HC3 heteroscedasticity-robust covariances** (e.g., using `statsmodels.regression.linear_model.OLS.fit(cov_type='HC3')` in Python).
3.  **Interpret Coefficients:**
    *   The coefficient $\beta$ represents the change in CLR abundance (log-ratio scale relative to geometric mean) per one-standard-deviation increase in the host score.
4.  **Multiple Testing Correction:**
    *   Gather all 122 raw $P$-values.
    *   Apply Benjamini–Hochberg FDR correction.
    *   Significance threshold is locked at $q < 0.05$.

---

## Step 4: Supporting Association Methods

To cross-validate the primary OLS model, execute the following supporting analyses:

1.  **Spearman Rank Correlation:**
    *   Calculate Spearman's $\rho$ and raw $P$-values.
    *   Construct 95% confidence intervals using 2,000 bootstrap resamples.
2.  **Permutation-Based Association Test:**
    *   Shuffle the host scores 9,999 times.
    *   Fit OLS for each shuffle and construct the empirical distribution of $t$-statistics.
    *   Calculate the empirical two-sided $P$-value.
3.  **MaAsLin2 Association:**
    *   Fit MaAsLin2 using the R/Python package.
    *   Settings: `normalization = "NONE"`, `transform = "NONE"`.
    *   Verify that the resulting coefficient sign and magnitude correspond to the primary OLS model.

---

## Step 5: Presence/Absence Logistic Regression

For genera that have moderate prevalence (preventing both all-zero OLS artifacts and logistic separation), fit a presence/absence model:

1.  **Taxon Selection:** 
    *   Only genera detected in at least 10 samples **and** absent in at least 10 samples (16.1% to 83.9% prevalence).
2.  **Binarization:**
    *   Convert Bracken abundance to binary: $y_i = 1$ if abundance $> 0.0$, else $y_i = 0$.
3.  **Model Specification:**
    *   Fit logistic regression:
        $$\ln\left(\frac{p_i}{1-p_i}\right) = \beta_0 + \beta_1 \cdot standardized\_Moffitt50\_contrast$$
    *   Check for separation: If the model fails to converge or exhibits extreme standard errors ($> 15$), exclude the genus and report it as mathematically unstable.

---

## Step 6: Contamination & Technical Sensitivity Checks

To determine if findings are robust to laboratory contamination or technical proxies:

1.  **Contaminant Exclusion Runs:**
    *   Rerun OLS models on `sensitivity/MICRO_SENS_NO_HIGH_RISK_centered_log_ratio.tsv.gz` (excluding 6 high-risk genera).
    *   Rerun OLS models on `sensitivity/MICRO_SENS_NO_CONTAMINANTS_centered_log_ratio.tsv.gz` (excluding 12 high- and moderate-risk genera).
2.  **LOGO Checks:**
    *   For any significant genus, recalculate the CLR values of the remaining 121 genera without that genus. Rerun the OLS model to confirm the association holds.
3.  **Matrix Total-Abundance Proxy Check:**
    *   Fit Model 1: `CLR_genus ~ standardized_Moffitt50_contrast + log10_matrix_total_abundance_proxy`.
    *   Confirm the coefficient $\beta_1$ for the host outcome remains stable (direction and significance) after accounting for the abundance proxy.
    *   Correlate each genus's CLR with the total-abundance proxy. Report correlation coefficients.
4.  **Host TME Sensitivity Models:**
    *   Fit Model 3P: `CLR_genus ~ standardized_Moffitt50_contrast + inferred_tumor_purity`.
    *   Fit Model 3I: `CLR_genus ~ standardized_Moffitt50_contrast + immune_score`.
    *   Fit Model 3S: `CLR_genus ~ standardized_Moffitt50_contrast + stromal_score`.
    *   Treat Models 3P/3I/3S as robustness checks only; Model 0 remains primary.
    *   Do not combine inferred tumor purity, immune score, stromal score, and ESTIMATE score in one model. Phase 7A.5 blocked the combined TME model because ESTIMATE score and inferred tumor purity are mathematically collinear and the combined screen had maximum VIF = Inf and condition number = 1.04e15.
    *   Interpretation guardrail: ESTIMATE-derived purity and cell-composition scores are inferred from the same host transcriptomic matrix as the Moffitt score and are not independent experimental measurements. Adjustment may remove transcriptional-state biology rather than only confounding.

---

## Step 7: Preprocessing Sensitivity Concordance

For every genus meeting the primary FDR threshold ($q < 0.05$), evaluate the association across all 9 pre-computed sensitivity matrices. Summarize the results in a concordance table reporting:

1.  **Sign Consistency:** Do all models agree on the direction of effect (+ or -)?
2.  **OLS Beta Range:** $[\beta_{min}, \beta_{max}]$ across all models.
3.  **FDR q Range:** $[q_{min}, q_{max}]$ across all models.
4.  **Concordance Fraction:** Number of sensitivity models with $P < 0.05$ divided by 9.
5.  **Transform Dependency:** Note if the association becomes non-significant ($P \ge 0.05$) under rCLR or when extreme samples are excluded.

---

## Step 8: Sample-Level Influence Diagnostics

Before finalizing findings, audit for individual patient leverage:

1.  **Cook's Distance:** Calculate Cook's $D_i$ for each sample in the OLS regression. Flag samples with $D_i > 4/n$ ($D_i > 0.0645$).
2.  **DFBETAs:** Calculate DFBETA for each sample. Flag samples with $|DFBETA_i| > 2/\sqrt{n}$ ($|DFBETA_i| > 0.254$).
3.  **Leave-One-Sample-Out:** Iteratively fit OLS leaving one patient out. Check if the raw $P$-value exceeds 0.05 in any run.
4.  **Extreme-Sample Exclusion Run:** Evaluate OLS coefficients on `sensitivity/MICRO_SENS_EXCLUDE_EXTREME_centered_log_ratio.tsv.gz` (excluding `Basal-like1`, `Hybrid18`, `Hybrid23`).
5.  *Guardrail:* If an association is driven by $\le 2$ patients, flag it as "leverage-sensitive" and do not grade it as robust. Do not delete influential samples to force significance.

---

## Step 9: Evidence Grading

Grade all candidate associations into these five predefined classes:

*   **ROBUST_ASSOCIATION:**
    *   Primary OLS BH $q < 0.05$.
    *   Bootstrap 95% CI for OLS coefficient excludes zero.
    *   Same direction in Spearman rank correlation.
    *   Same direction in at least 7 of 9 preprocessing sensitivity runs.
    *   Not driven solely by the removal of one technical sample.
    *   Evaluated for contamination status and not flagged as a high-risk contaminant.
*   **SUGGESTIVE_ASSOCIATION:**
    *   Primary OLS BH $q$ is between 0.05 and 0.10, OR
    *   Primary OLS BH $q < 0.05$ but fails robustness criteria (e.g. sign reversal in rCLR or driven by 1 sample).
*   **METHOD_SENSITIVE:**
    *   OLS sign changes or significance disappears completely under alternative prevalence/pseudocount choices.
*   **CONTAMINATION_SENSITIVE:**
    *   Signal is lost when high- or moderate-risk contaminant taxa are removed, or genus abundance is strongly correlated with the total-abundance proxy.
*   **NO_SUPPORTED_ASSOCIATION:**
    *   Fails to meet primary or suggestive thresholds.
*   **TO_VERIFY:**
    *   Uncertainty in metadata alignment or implementation.
