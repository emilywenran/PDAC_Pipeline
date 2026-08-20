# Risk Register

This document lists potential scientific, bioinformatic, and statistical risks, along with their likelihood, impact, mitigation plans, and monitoring protocols.

---

## 1. Risk Summary Matrix

| Risk ID | Description | Likelihood | Impact | Status |
| :--- | :--- | :--- | :--- | :--- |
| **R-01** | Small sample size ($n=62$) and high-dimensional predictors | High | High | Open |
| **R-02** | Circularity in predicting transcriptome-derived subtype labels | High | High | Open |
| **R-03** | Data leakage during cross-validation / resampling | Medium | High | Open |
| **R-04** | Instability and artificial partitioning of the "hybrid" subtype | High | Medium | Open |
| **R-05** | Low-biomass microbiome contamination ("kitome" and "screentome") | High | High | Open |
| **R-06** | Microbiome compositionality confounding standard statistics | High | High | Open |
| **R-07** | Tumor purity and cellular-composition confounding | High | High | Open |
| **R-08** | Technical batch effects in RNA-seq or shotgun metagenomic sequencing | Medium | Medium | Open |
| **R-09** | Missing clinical information (survival, stage, demographics) | High | Medium | Open |
| **R-10** | Absence of an identical external multi-omics validation cohort | High | High | Open |
| **R-11** | Overinterpretation of cross-sectional correlation as causation | High | High | Open |
| **R-12** | Multiple-testing burden inflating false discovery rates | High | High | Open |

---

## 2. Detailed Risk Descriptions & Mitigation Strategies

### R-01: Small Sample Size ($n=62$) and High-Dimensional Predictors
*   **Description:** With $n=62$ samples and thousands of transcriptomic and microbial features, models are highly prone to overfitting.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Restrict multivariate models to a maximum of 3-5 covariates to prevent overfitting (rule of thumb: $\ge 10$ samples per predictor).
    *   Prioritize dimension reduction (NMF, PCA) and signature-based scoring over individual gene-level modeling.
    *   Avoid deep learning or complex non-linear classifiers.
*   **Monitoring Strategy:** Track the degrees of freedom in all regression models. Require standard errors and confidence intervals to be reported alongside p-values.

### R-02: Circularity in Predicting Transcriptome-Derived Subtype Labels
*   **Description:** Subtype labels (basal-like vs. classical) are defined based on clustering of the transcriptome. Training a classifier on the same transcriptomic features to predict these labels creates a circular loop, yielding inflated classification metrics.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Do not train classifiers to predict host transcriptomic subtypes using the same host transcriptomic features.
    *   Define subtypes as continuous projection scores (Moffitt axis score) and test associations with non-transcriptomic variables (microbiome, clinical endpoints).
    *   If transcriptomic classification models are built, train on external datasets (e.g. TCGA) and project onto GSE172356 as an independent test.
*   **Monitoring Strategy:** Code audits of all classification workflows to ensure no model uses overlapping training features and target labels.

### R-03: Data Leakage during Cross-Validation
*   **Description:** Performing normalization, scaling, or feature selection on the entire dataset prior to splitting into cross-validation folds leaks information from the test set into the training set, leading to over-optimistic performance.
*   **Likelihood:** Medium | **Impact:** High
*   **Mitigation Strategy:** 
    *   Use pipeline frameworks (e.g., `tidymodels` in R, `scikit-learn Pipelines` in Python) that strictly split data before applying any transformations.
    *   Perform feature selection and parameter tuning entirely within the training folds.
*   **Monitoring Strategy:** Inspect the analysis scripts to verify that data splitting occurs as the very first step in any machine learning or resampling loop.

### R-04: Instability of the "Hybrid" Subtype
*   **Description:** The "hybrid" class may represent a statistical artifact of forcing a continuous distribution into three groups, leading to unstable classification and poor reproducibility.
*   **Likelihood:** High | **Impact:** Medium
*   **Mitigation Strategy:** 
    *   Evaluate stability using silhouette width analysis, Jaccard bootstrapping, and comparisons with random gene sets.
    *   Treat the basal-classical spectrum as a continuous gradient (e.g. using ssGSEA score subtraction) in primary analyses, keeping discrete labels as secondary descriptors.
*   **Monitoring Strategy:** Plot silhouette width distributions. Highlight samples with negative silhouette widths as "unstable."

### R-05: Low-Biomass Metagenomic Contamination
*   **Description:** Intratumoral microbiomes are characterized by low bacterial biomass. Reagent, kit, and laboratory contaminants can dominate sequencing reads and generate false biological associations.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Verify the availability of and download all negative controls (extraction, PCR, environmental) from PRJNA719915.
    *   Run the R package `decontam` (prevalence and frequency methods) to computationally identify and remove contaminant features before downstream modeling.
*   **Monitoring Strategy:** Report the number and percentage of reads/features identified as contaminants. Map identified contaminants against known database lists of common laboratory contaminants.

### R-06: Microbiome Compositionality
*   **Description:** Sequencing reads represent relative proportions, not absolute abundances. Standard statistical tests (e.g., Pearson correlation, t-tests, standard linear regression) on raw relative abundances yield spurious correlations and false positives.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Apply Centered Log-Ratio (CLR) transformation to taxons.
    *   Utilize composition-aware differential abundance and regression tools, specifically ALDEx2, ANCOM-BC, and MaAsLin2.
*   **Monitoring Strategy:** Review all code to ensure no direct regression or correlation is performed on raw relative abundances without compositional adjustment.

### R-07: Tumor Purity and Cellular-Composition Confounding
*   **Description:** Bulk RNA-seq is a mixture of cells. Basal-like signatures are associated with high inflammatory stroma. The abundance of specific microbes might correlate with the proportion of stromal or immune cells rather than the tumor subtype.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Estimate tumor purity (ESTIMATE) and cell fractions (ConsensusTME).
    *   Include tumor purity and immune scores as covariates in the MaAsLin2 association models to evaluate if microbial associations are independent of tumor composition.
*   **Monitoring Strategy:** Report and compare association results with and without adjusting for tumor purity.

### R-08: Technical Batch Effects
*   **Description:** Technical variation across sequencing runs, extraction batches, or center protocols can obscure biological variance.
*   **Likelihood:** Medium | **Impact:** Medium
*   **Mitigation Strategy:** 
    *   Audit metadata for sequencing batch indicators.
    *   Use PVCA to quantify batch-associated variance.
    *   If batch effects are significant, apply `ComBat` or add batch as a covariate.
*   **Monitoring Strategy:** Generate PCA and UMAP plots colored by technical batches; inspect for sample separation by batch.

### R-09: Missing Clinical Information
*   **Description:** GSE172356 clinical tables may contain missing values for staging, survival, or treatment, limiting survival modeling.
*   **Likelihood:** High | **Impact:** Medium
*   **Mitigation Strategy:** 
    *   Audit the clinical metadata in `01_metadata/` immediately.
    *   Label missing fields as `TO VERIFY`.
    *   Apply survival analysis (e.g., Cox proportional hazards) only to the subset of samples with verified survival times and censoring indicators.
*   **Monitoring Strategy:** Maintain a completeness table in the project logs tracking the percentage of missing values for all clinical variables.

### R-10: Absence of an Identical External Multi-Omics Cohort
*   **Description:** An identical external paired transcriptome–microbiome cohort is not assumed to exist, and any TCGA-derived microbiome signal must be treated as exploratory due to batch and contamination concerns.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Perform host validation (subtyping, continuous gradient, and pathways) in large cohorts (TCGA-PAAD, ICGC) where host transcriptomics is highly robust.
    *   Treat any TCGA-derived microbiome signal or external paired analysis as exploratory (Tier C/D) rather than direct replication, and carefully audit the data for batch effects and contamination.
    *   Compare identified taxa with published literature on the PDAC microbiome (e.g., Nejman et al., 2020) at the genus/species level, treating them as exploratory candidate findings.
*   **Monitoring Strategy:** Explicitly separate host validation results (highly robust) from microbiome validation results (strictly exploratory/literature-based) in all reports.

### R-11: Overinterpretation of Correlations
*   **Description:** Cross-sectional correlation between the microbiome and tumor subtypes can be easily overinterpreted as a causal mechanism (e.g., "bacteria drive the basal-like subtype").
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Strictly enforce the Manuscript Guardrails.
    *   Use observational, associative language (e.g., "is associated with," "correlates with," "is enriched in") instead of causal verbs ("causes," "drives," "promotes," "induces").
*   **Monitoring Strategy:** Conduct peer review of all manuscript drafts and reports to flag causal language.

### R-12: Multiple-Testing Burden
*   **Description:** Testing association across thousands of genes and hundreds of taxonomic features results in massive numbers of statistical tests, inflating false positives if uncorrected.
*   **Likelihood:** High | **Impact:** High
*   **Mitigation Strategy:** 
    *   Apply Benjamini-Hochberg False Discovery Rate (FDR) correction to all differential expression, differential abundance, and regression analyses.
    *   Pre-specify the primary taxa of interest (*Acinetobacter*, *Pseudomonas*, *Sphingopyxis*) to narrow the testing hypothesis space.
*   **Monitoring Strategy:** Report adjusted p-values ($q$-values) for all statistical outputs. Reject any association with $q \ge 0.05$.
