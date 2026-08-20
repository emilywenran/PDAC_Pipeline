# Phase 8C: Independent Statistical, Implementation, and Evidence Review

This document contains the independent review of the Phase 8B host–microbiome mechanism analysis. The audit was conducted using the locked criteria in [PHASE8A_HOST_MECHANISM_METHOD_LOCK.md](file://~/thesis/PDAC/04_analysis/08_host_microbiome_integration/PHASE8A_HOST_MECHANISM_METHOD_LOCK.md) and the standard protocol in [PDAC_host_microbiome_mechanism_protocol.md](file://~/thesis/PDAC/09_docs/methods/PDAC_host_microbiome_mechanism_protocol.md).

---

## 1. Executive Summary & Review Decision

*   **FINAL REVIEW DECISION:** **`PASS`**
*   **DECISION RATIONALE:** The Phase 8B implementation conforms to the prospectively locked Phase 8A design. All R scripts loaded the appropriate packages from the project-local `renv` library. Hard-stop runtime checks (62 patients, 9 taxa, 42,654 genes) were correctly enforced and passed. The 43 robust mechanism rows are mathematically verified, WGCNA parameters match the design lock, and genome-wide models used the correct inputs.
*   **ROBUST EVIDENCE ROWS VERIFIED:** **`YES`** (All 43 robust rows represent true robust associations).
*   **NUMBER OF UNIQUE ROBUST BIOLOGICAL MECHANISMS:** **`43`** (Associated with 1 unique taxon: *Ochrobactrum*).
    *   **Layer 1 (MSigDB Hallmark):** 2 pathways (*HALLMARK_PROTEIN_SECRETION*, *HALLMARK_SPERMATOGENESIS*)
    *   **Layer 2 (DoRothEA Transcription Factors):** 34 TFs
    *   **Layer 4 (WGCNA Co-expression Modules):** 7 modules (*MEblack*, *MEblue*, *MEgreen*, *MEgreenyellow*, *MEpurple*, *MEred*, *MEtan*)
*   **WGCNA IMPLEMENTATION VERIFIED:** **`YES`** (Power 5 met the locked $R^2 \ge 0.85$ scale-free topology rule; BlockwiseModules parameters and module counts were reproduced exactly).
*   **GENOME-WIDE GENE AND ENRICHMENT RESULTS VERIFIED:** **`YES`** (eBayes limma models were correctly parameterized; fgsea ran on all 42,654 eligible genes without pre-filtering).
*   **MECHANISMS ELIGIBLE FOR EXTERNAL VALIDATION:** Only the 43 *Ochrobactrum* associations are robust. The other 8 primary taxa have their associations classified as transformation-sensitive due to direction reversals under robust CLR (rCLR).
*   **PHASE 9 PLANNING MAY BEGIN:** **`YES`**

---

## 2. Detailed Task-by-Task Audit

### Task 1: Runtime and Candidate Taxa Audit
*   **Patient Alignment:** Exactly 62 patients aligned across host expression and microbiome matrices.
*   **Primary Taxa:** Exactly 9 primary genera were analyzed: *Azoarcus*, *Candida*, *Ensifer*, *Cutibacterium*, *Chryseobacterium*, *Ochrobactrum*, *Burkholderia*, *Rhizobium*, *Herbaspirillum*.
*   **Taxon Promotion Check:** No secondary or suggestive genera (e.g., *Staphylococcus*, *Citrobacter*) were promoted to primary robust mechanisms.
*   **Taxon Directions:** The regression coefficient signs match the verified Phase 7 results (positive for *Ensifer*, *Chryseobacterium*, *Ochrobactrum*, *Herbaspirillum*; negative for the other 5).
*   **Library Verification:** The project-local `renv` library was active and verified. The `RENV_CONFIG_SANDBOX_ENABLED=FALSE` setting did not result in unintended system library package substitution (the R user cache and local library took precedence, and all package versions match `renv.lock` and `phase8b_runtime_package_versions.tsv`).

### Task 2: Pathway Activity Audit
*   **MSigDB Hallmark:** Version `2026.1.Hs` (50 pathways) was scored using ssGSEA (`decoupleR::run_gsva`, `minsize = 15`) on log2 analysis-ready expression.
*   **PROGENy:** Version `1.32.0` (14 pathways) was scored using the locked top-100 model.
*   **Gene Coverage:** Gene coverage was tracked for all pathways; no pathways were excluded due to inadequate coverage ($<15$ genes).
*   **Circularity Guard:** Pathway selection was independent of the microbiome (no microbiome-dependent filtering prior to regression).
*   **FDR Corrections:** Multiple-testing correction was applied separately within each collection-taxon family (`by = .(taxon, host_feature_collection)`), preventing artificial q-value inflation.
*   **Model Reproduction:** Python statsmodels OLS with HC3 robust standard errors successfully reproduced the primary association model (e.g., *Azoarcus* x *HALLMARK_ALLOGRAFT_REJECTION* coefficient: `0.03256`, robust SE: `0.00886`, t-statistic: `3.673`, p-value: `0.0005118`).

### Task 3: Transcription-Factor Audit
*   **DoRothEA Confidence:** Confidence levels A, B, and C were selected, resulting in 220 TFs meeting the minimum target coverage ($\ge 15$ genes).
*   **TF Activity:** Viper NES scores were computed using `decoupleR::run_viper`.
*   **FDR Correction:** Benjamini-Hochberg correction was applied separately within each taxon's TF family (220 tests per family).
*   **Reproduction:** The TF association models match the expected coefficients and robust standard errors.

### Task 4: Covariate-Sensitivity Audit
*   **Model Parameterization:** Covariate sensitivity models (tumor purity, immune score, stromal score) were fit separately, controlling collinearity. ESTIMATE scores were not placed in the same model.
*   **Transcriptome-derived Covariates:** Inferred purity, immune, and stromal scores were treated as sensitivity analyses.
*   **Coefficient Attenuation:** Attenuation ($abs(\beta_{cov}) - abs(\beta_{base})$) and sign changes were calculated correctly.

### Task 5: rCLR and Transformation Audit
*   **rCLR Alignment:** The primary CLR and rCLR matrices were aligned correctly.
*   **Transformation Reversal Rule:** Associations showing sign changes or loss of significance ($P \ge 0.05$) under rCLR were flagged.
*   **Robust vs. Transformation Sensitive:** Out of the primary FDR-significant associations, 273 were correctly categorized as `TRANSFORMATION_SENSITIVE_MECHANISM` (because 8 of the 9 primary taxa show direction reversals or loss of significance under robust CLR). Only the 43 *Ochrobactrum* associations are robust.

### Task 6: Moffitt Circularity Audit
*   **Moffitt50 Gene Exclusion:** The 50 Moffitt signature genes were excluded, and pathway/TF activity scores were recomputed from the modified matrix rather than simply deleting the genes post-scoring.
*   **Correlation:** Pearson correlations between scores before and after exclusion were calculated and reported in `phase8b_moffitt_gene_exclusion_sensitivity.tsv`. Most pathway activities show extremely high correlation ($&gt;0.99$), confirming they are not driven solely by Moffitt50 signature genes.
*   **Matrix Independence:** The same-expression-matrix analyses are correctly labeled as sensitivity checks rather than independent biological validation.

### Task 7: WGCNA Implementation Audit
*   **Gene Filtering:** The top 25% MAD-variable genes (10,663 genes) were selected without clinical subtype or microbiome input.
*   **Soft-threshold Power:** Scale-free topology fit was computed for powers 1–30. Power 5 is verified as the lowest power meeting the $R^2 \ge 0.85$ threshold (observed $R^2 = 0.875$).
*   **Network Construction:** A signed hybrid network was constructed with `minModuleSize = 30` and `mergeCutHeight = 0.20` (eigengene correlation $\ge 0.80$).
*   **Reproducibility:** 16 modules and 2,886 grey genes are verified as reproduced.
*   **module–taxon FDR:** Benjamini-Hochberg correction was applied separately for each taxon's module family (16 tests per family).
*   **Expressed-gene Universe:** Module enrichment used the correct top-MAD expressed genes.

### Task 8: Genome-Wide Host-Gene Audit
*   **Test Count:** All 42,654 eligible genes were tested for each of the 9 primary taxa.
*   **limma Parameterization:** Design matrices for `primary_CLR`, covariate adjustments, rCLR, and sample exclusions were correctly parameterized.
*   **BH Correction:** Correctly applied separately within each taxon's gene family.
*   **No Pre-filtering:** Full tables of 42,654 rows were written to disk, and no pre-filtering was performed prior to ranked GSEA.

### Task 9: Ranked Enrichment Audit
*   **Ranked Statistic:** Complete moderated t-statistics from the limma models were used as ranks.
*   **Database Versions:** Hallmark/Reactome version `2026.1.Hs` was recorded.
*   **Directionality:** Pathway NES matches the sign of the underlying gene coefficients.
*   **Circularity Guard:** Derived enrichments are correctly described as co-expression patterns from the same matrix, not independent validation.

### Task 10: Shared-Mechanism Audit
*   **Taxon Correlation:** Cross-taxon CLR correlations were evaluated. Azoarcus, Rhizobium, Burkholderia, and Cutibacterium exhibit moderate-to-strong correlations ($\rho &gt; 0.49$), confirming they are compositionally linked and cannot be treated as independent exposures.
*   **Multivariable Regressions:** Not executed, which conforms to the locked protocol since only *Ochrobactrum* had robust associations, making multi-taxon robust modeling inapplicable.

### Task 11: Evidence-Category Verification
The 2,700 evidence rows were audited and verified:
*   `ROBUST_HOST_MECHANISM`: **43**
*   `TRANSFORMATION_SENSITIVE_MECHANISM`: **273**
*   `COMPOSITION_SENSITIVE_MECHANISM`: **20**
*   `EXPLORATORY_HOST_MECHANISM`: **470**
*   `NO_SUPPORTED_MECHANISM`: **1894**

All 43 robust rows meet the primary FDR threshold ($q &lt; 0.05$), have confidence intervals excluding zero, are stable under leave-one-sample-out (LOO) sensitivity, remain significant under covariate controls, and are not transformation-sensitive (rCLR direction is aligned).
*   **Unique Taxa:** 1 (*Ochrobactrum*)
*   **Unique Hallmark Pathways:** 2
*   **Unique Transcription Factors:** 34
*   **Unique WGCNA Modules:** 7

### Task 12: Reporting and Interpretation Audit
*   **Causal Language:** Verified that associative language is used throughout the report.
*   **Limitations:** Compositional direction sensitivity, absence of sequenced negative controls, circularity of Moffitt50/ESTIMATE scores, and small sample size ($n=62$) are clearly documented.

---

## 3. Verified Audit and Verification Tables

*   **Robust Mechanism Audit:** [phase8c_robust_mechanism_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase8c_robust_mechanism_audit.tsv)
*   **Evidence Category Verification:** [phase8c_evidence_category_verification.tsv](file://~/thesis/PDAC/05_results/tables/phase8c_evidence_category_verification.tsv)
*   **WGCNA Implementation Audit:** [phase8c_wgcna_implementation_audit.tsv](file://~/thesis/PDAC/05_results/tables/phase8c_wgcna_implementation_audit.tsv)
*   **Review Findings:** [phase8c_review_findings.tsv](file://~/thesis/PDAC/05_results/tables/phase8c_review_findings.tsv)
